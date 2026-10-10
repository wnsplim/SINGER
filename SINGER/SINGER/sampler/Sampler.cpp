//
//  Sampler.cpp
//  SINGER
//
//  Created by Yun Deng on 3/31/23.
//

#include "Sampler.hpp"
#include "TSP.hpp"

Sampler::Sampler() {}

Sampler::Sampler(double pop_size, double r, double m) {
    Ne = pop_size;
    mut_rate = m*pop_size;
    recomb_rate = r*pop_size;
}

void Sampler::set_precision(double c, double q) {
    bsp_c = c;
    tsp_q = q;
}

void Sampler::set_input_file_prefix(string f) {
    data.input = f;
}

void Sampler::set_output_file_prefix(string f) {
    output_prefix = f;
}

Node_ptr Sampler::new_sample(int i) {
    Node_ptr n = new_node(tip_times.empty() ? 0.0 : tip_times[i]);
    n->is_sample = true;
    n->set_index(i);
    return n;
}

void Sampler::read_tip_ages(string filename, double g) {
    data.ploidy = ploidy;
    tip_times = data.read_tip_ages(filename);
    tip_offset = *min_element(tip_times.begin(), tip_times.end())/(g*Ne);
    for (double &x : tip_times) {
        x = x/(g*Ne) - tip_offset;
    }
}

void Sampler::load_data(double start_pos, double end_pos) {
    data.ploidy = ploidy;
    data.read(start_pos, end_pos, true);
    vector<Node *> leaves;
    for (int i = 0; i < (int) data.derived_sites.size(); i++) {
        Node_ptr n = new_sample(i);
        for (double x : data.derived_sites[i]) {
            n->add_mutation(x);
        }
        sample_nodes.insert(n);
        leaves.push_back(n.get());
    }
    order_samples();
    sequence_length = end_pos - start_pos;
    set_missing(leaves);
}

void Sampler::set_missing(vector<Node *> &leaves) {
    for (int i = 0; i < (int) leaves.size(); i++) {
        leaves[i]->missing_sites = data.missing_sites[i];
        leaves[i]->masked_intervals = data.masked_intervals[i];
    }
}

void Sampler::order_samples() {
    ordered_sample_nodes = vector<Node_ptr>(sample_nodes.begin(), sample_nodes.end());
    random_engine.seed(random_seed);
    shuffle(ordered_sample_nodes.begin(), ordered_sample_nodes.end(), random_engine);
}

void Sampler::build_singleton_arg() {
    Node_ptr n = *ordered_sample_nodes.begin();
    arg = ARG(Ne, sequence_length);
    arg.time_offset = tip_offset;
    arg.any_missing = data.any_missing;
    arg.unassayed_sites = data.unassayed_sites;
    if (recomb_map.coordinates.empty() and mut_map.coordinates.empty()) {
        arg.discretize(min(max(1.0, round(rho_unit/recomb_rate)), 100.0));
    } else {
        arg.discretize(recomb_map, mut_map, start, rho_unit/Ne, recomb_rate/Ne);
    }
    arg.build_singleton_arg(n);
    arg.compute_rhos_thetas(recomb_rate, mut_rate, recomb_map, mut_map, start);
    for (const Node_ptr &s : sample_nodes) {
        for (auto &x : s->mutation_sites) {
            if (x.second == 1 and x.first < sequence_length and arg.thetas[arg.get_index(x.first)] == 0) {
                cerr << "Error: the VCF has a variant at position " << (long long) (start + x.first) << ", where the mutation rate is 0. " << endl;
                exit(1);
            }
        }
    }
    arg.masked = data.masked;
    arg.compute_assayed();
    vector<Node *> leaves(sample_nodes.size());
    for (const Node_ptr &s : sample_nodes) {
        leaves[s->index] = s.get();
    }
    set_entries(leaves);
}

void Sampler::set_entries(vector<Node *> &leaves) {
    arg.genotype_likelihoods.clear();
    for (const Genotype_entry &e : data.entries) {
        Genotype_likelihood g;
        g.pos = e.pos;
        g.leaves[0] = leaves[e.haplotypes[0]];
        g.leaves[1] = e.haplotypes[1] < 0 ? nullptr : leaves[e.haplotypes[1]];
        copy(e.L, e.L + 4, g.L);
        arg.genotype_likelihoods.push_back(g);
    }
    arg.index_genotype_likelihoods();
    for (double x : data.entry_sites) {
        arg.mutation_sites.insert(x);
        arg.mutation_branches[x];
    }
}

void Sampler::iterative_start() {
    start_log();
    build_singleton_arg();
    auto it = ordered_sample_nodes.begin();
    it++;
    Threader_smc threader = Threader_smc(bsp_c, tsp_q);
    while (it != ordered_sample_nodes.end()) {
        random_engine.seed(random_seed);
        threader.reset();
        threader.pe->penalty = penalty;
        threader.be->penalty = penalty;
        threader.pe->ancestral_prob = polar;
        threader.be->ancestral_prob = polar;
        arg.penalty = penalty;
        arg.ancestral_prob = polar;
        Node_ptr n = *it;
        threader.thread(arg, n);
        arg.check_incompatibility();
        cout << "Number of flippings: " << arg.count_flipping() << endl;
        it++;
        random_seed = random_engine();
        write_iterative_start();
        int placed = (int) arg.sample_nodes.size();
        if ((placed & (placed - 1)) == 0 or it == ordered_sample_nodes.end()) {
            sweep_samples();
        }
    }
    arg.resample_genotypes();
    cout << "orignal ARG length: " << arg.get_arg_length() << endl;
    rescale();
    cout << "rescaled ARG length: " << arg.get_arg_length() << endl;
    string node_file = output_prefix + "_start_nodes_" + to_string(sample_index) + ".txt";
    string branch_file= output_prefix + "_start_branches_" + to_string(sample_index) + ".txt";
    string recomb_file = output_prefix + "_start_recombs_" + to_string(sample_index) + ".txt";
    string mut_file = output_prefix + "_start_muts_" + to_string(sample_index) + ".txt";
    arg.write(node_file, branch_file, recomb_file, mut_file);
    string coord_file = output_prefix + "_coordinates.txt";
    arg.write_coordinates(coord_file);
}

void Sampler::sweep_samples() {
    Threader_smc threader = Threader_smc(bsp_c, tsp_q);
    vector<Node *> nodes(arg.sample_nodes.begin(), arg.sample_nodes.end());
    for (Node *n : nodes) {
        random_engine.seed(random_seed);
        threader.reset();
        threader.pe->penalty = penalty;
        threader.be->penalty = penalty;
        threader.pe->ancestral_prob = polar;
        threader.be->ancestral_prob = polar;
        threader.rethread_sample(arg, n);
        random_seed = random_engine();
    }
    arg.release_dead_nodes();
    write_iterative_start();
}

void Sampler::internal_sample(int num_iters, int spacing) {
    Threader_smc threader = Threader_smc(bsp_c, tsp_q);
    while (sample_index < num_iters) {
        cout << get_time() << " Iteration: " << to_string(sample_index) << endl;
        int moves = 0;
        cout << "Random seed: " << random_seed << endl;
        random_engine.seed(random_seed);
        while (moves < spacing) {
            threader.reset();
            threader.pe->penalty = penalty;
            threader.be->penalty = penalty;
            threader.pe->ancestral_prob = polar;
            threader.be->ancestral_prob = polar;
            arg.penalty = penalty;
            arg.ancestral_prob = polar;
            if (exact) {
                tuple<double, Branch, double> cut_point = arg.sample_uniform_cut();
                threader.exact_internal_rethread(arg, cut_point);
            } else {
                tuple<double, Branch, double> cut_point = arg.sample_internal_cut();
                threader.internal_rethread(arg, cut_point);
            }
            moves += 1;
            arg.clear_remove_info();
        }
        arg.resample_genotypes();
        arg.release_dead_nodes();
        rescale();
        random_seed = random_engine();
        write_sample();
        arg.check_incompatibility();
        cout << "Start: " << arg.start << " , End: " << arg.end << endl;
        string node_file = output_prefix + "_nodes_" + to_string(sample_index) + ".txt";
        string branch_file= output_prefix + "_branches_" + to_string(sample_index) + ".txt";
        string recomb_file = output_prefix + "_recombs_" + to_string(sample_index) + ".txt";
        string mut_file = output_prefix + "_muts_" + to_string(sample_index) + ".txt";
        sample_index += 1;
        arg.write(node_file, branch_file, recomb_file, mut_file);
        cout << "Number of trees: " << arg.recombinations.size() << endl;
        cout << "Number of flippings: " << arg.count_flipping() << endl;
    }
}

void Sampler::resume_internal_sample(int num_iters, int spacing) {
    string log_file = output_prefix + ".log";
    read_resume_point(log_file);
    vector<string> words = Data_reader::read_last_line(log_file);
    if (!seed_given) {
        random_seed = stoi(words[words.size() - 2]);
    }
    sample_index += 1;
    arg.check_incompatibility();
    cout << "Number of trees: " << arg.recombinations.size() << endl;
    cout << "Number of flippings: " << arg.count_flipping() << endl;
    internal_sample(num_iters, spacing);
}

void Sampler::debug_resume_internal_sample(int num_iters, int spacing) {
    retract_log(5);
    string log_file = output_prefix + ".log";
    vector<string> words = Data_reader::read_last_line(log_file);
    if (words.size() == 0 or words[2] == "initial_thread" or words[0] == "Time") {
        cout << "new seed: " << random_seed << endl;
        sample_index = 0;
        load_data(start, end);
        iterative_start();
        internal_sample(num_iters, spacing);
    } else {
        read_resume_point(log_file);
        sample_index += 1;
        cout << "new seed: " << random_seed << endl;
        internal_sample(num_iters, spacing);
    }
}

void Sampler::rescale() {
    if (scaling_rep == 0) {
        return;
    }
    Scaler scaler = Scaler();
    scaler.num_windows = scaling_bin;
    scaler.compute_deltas(arg, mut_map, start, mut_rate/Ne);
    for (int i = 0; i < scaling_rep; i++) {
        scaler.reset();
        scaler.rescale(arg, mut_rate);
    }
}

void Sampler::start_log() {
    string filename = output_prefix + ".log";
    ofstream file(filename, ios::out|ios::trunc);
    if (!file) {
        cerr << "Error opening the file: " << filename << endl;
        return;
    }
    file << "Time" << "\t"
    << "Iteration:" << "\t"
    << "Threading_type" << "\t"
    << "#Recombinations" << "\t"
    << "#Mutations_not_uniquely_mapped" << "\t"
    << "Last_updated_pos" << "\t"
    << "Random_seed" << "\t"
    << "Counter" << endl;
    file.close();
}

void Sampler::write_iterative_start() {
    string filename = output_prefix + ".log";
    ofstream file(filename, ios::out|ios::app);
    if (!file) {
        cerr << "Error opening the file: " << filename << endl;
        return;
    }
    file << get_time() << "\t"
    << arg.sample_nodes.size() << "\t"
    << "initial_thread" << "\t"
    << arg.recombinations.size() - 2 << "\t"
    << arg.num_unmapped() << "\t"
    << setprecision(numeric_limits<double>::max_digits10)
    << arg.end << "\t"
    << random_seed << "\t"
    << TSP::counter << endl;
}

void Sampler::write_sample() {
    string filename = output_prefix + ".log";
    ofstream file(filename, ios::out|ios::app);
    if (!file) {
        cerr << "Error opening the file: " << filename << endl;
        return;
    }
    file << get_time() << "\t"
    << sample_index << "\t"
    << "rethread" << "\t"
    << arg.recombinations.size() - 2 << "\t"
    << arg.num_unmapped() << "\t"
    << setprecision(numeric_limits<double>::max_digits10)
    << arg.end << "\t"
    << random_seed << "\t"
    << TSP::counter << endl;
}

void Sampler::load_resume_arg() {
    arg = ARG(Ne, sequence_length);
    string node_file, branch_file, recomb_file, mut_file, coord_file;
    node_file = output_prefix + "_nodes_" + to_string(sample_index) + ".txt";
    branch_file= output_prefix + "_branches_" + to_string(sample_index) + ".txt";
    recomb_file = output_prefix + "_recombs_" + to_string(sample_index) + ".txt";
    mut_file = output_prefix + "_muts_" + to_string(sample_index) + ".txt";
    coord_file = output_prefix + "_coordinates.txt";
    arg.num_samples = ploidy*(int) data.sample_names().size();
    arg.read(node_file, branch_file, recomb_file, mut_file);
    arg.read_coordinates(coord_file);
    vector<Node *> leaves(arg.sample_nodes.size());
    for (Node *n : arg.sample_nodes) {
        leaves[n->index] = n;
    }
    data.ploidy = ploidy;
    data.read(start, end, false);
    set_missing(leaves);
    arg.any_missing = data.any_missing;
    arg.unassayed_sites = data.unassayed_sites;
    arg.compute_rhos_thetas(recomb_rate, mut_rate, recomb_map, mut_map, start);
    arg.masked = data.masked;
    arg.compute_assayed();
    set_entries(leaves);
}

void Sampler::read_resume_point(string filename) {
    vector<string> words = Data_reader::read_last_line(filename);
    int log_length = (int) words.size();
    TSP::counter = stoi(words[log_length - 1]);
    sample_index = stoi(words[1]);
    load_resume_arg();
    arg.sequence_length = sequence_length;
    arg.end = stod(words[log_length - 3]);
    arg.end_tree = arg.get_tree_at(arg.end);
}

void Sampler::retract_log(int k) {
    const std::string file_path = output_prefix + ".log";
    std::ifstream in_file(file_path, std::ios::in | std::ios::ate);

    if (!in_file.is_open()) {
        std::cerr << "Unable to open log file: " << file_path << std::endl;
        return;
    }

    std::string first_line, second_line;
    in_file.seekg(0);
    std::getline(in_file, first_line);
    std::getline(in_file, second_line);

    in_file.seekg(0, std::ios::end);
    char c;
    int line_count = 0;
    long pos = in_file.tellg();

    while (pos > 0 && line_count < k) {
        in_file.seekg(--pos, std::ios::beg);
        in_file.get(c);
        if (c == '\n') {
            ++line_count;
        }
    }

    pos = (line_count < k) ? 0 : pos + 1;

    in_file.seekg(0, std::ios::beg);
    std::string content(pos, '\0');
    in_file.read(&content[0], pos);
    in_file.close();

    int remaining_lines = (int) std::count(content.begin(), content.end(), '\n') + (content.empty() ? 0 : 1);

    std::ofstream out_file(file_path, std::ios::out | std::ios::trunc);
    if (remaining_lines >= 2) {
        out_file << content;
    } else {
        out_file << first_line << '\n';
        if (!second_line.empty()) {
            out_file << second_line << '\n';
        }
    }
    out_file.close();
}

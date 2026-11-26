//
//  Sampler.cpp
//  SINGER
//
//  Created by Yun Deng on 3/31/23.
//

#include "Sampler.hpp"
#include "TSP.hpp"

Sampler::Sampler() {}

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
    approx_BSP::exact_calculator = *max_element(tip_times.begin(), tip_times.end()) > 0;
}

void Sampler::load_data(double start_pos, double end_pos) {
    data.ploidy = ploidy;
    data.read(start_pos, end_pos);
    vector<Node *> leaves;
    for (int i = 0; i < (int) data.derived_sites.size(); i++) {
        Node_ptr n = new_sample(i);
        sample_nodes.insert(n);
        leaves.push_back(n.get());
    }
    load_states(leaves);
    order_samples();
    sequence_length = end_pos - start_pos;
    set_missing(leaves);
}

void Sampler::load_states(vector<Node *> &leaves) {
    for (int i = 0; i < (int) leaves.size(); i++) {
        for (double x : data.derived_sites[i]) {
            leaves[i]->add_mutation(x);
        }
    }
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
    arg = ARG(Ne, sequence_length);
    arg.time_offset = tip_offset;
    arg.any_missing = data.any_missing;
    arg.unassayed_sites = data.unassayed_sites;
    if (recomb_map.coordinates.empty() and mut_map.coordinates.empty()) {
        arg.discretize(min(max(1.0, round(rho_unit/recomb_rate)), 100.0));
    } else {
        arg.discretize(recomb_map, mut_map, start, rho_unit/Ne, recomb_rate/Ne);
    }
    arg.build_singleton_arg(*ordered_sample_nodes.begin());
    arg.compute_rhos_thetas(recomb_rate, mut_rate, recomb_map, mut_map, start);
    for (const Node_ptr &s : sample_nodes) {
        for (auto &x : s->mutation_sites) {
            if (x.first >= 0 and x.first < sequence_length and arg.thetas[arg.get_index(x.first)] == 0) {
                cerr << "Error: the input has a variable site at position " << (long long) (start + x.first) << ", where the mutation rate is 0. " << endl;
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
    }
    arg.index_missing(leaves);
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
        arg.penalty = penalty;
        arg.ancestral_prob = polar;
        threader.thread(arg, *it);
        it++;
        random_seed = random_engine();
        arg.draw_mutations();
        arg.check_incompatibility();
        cout << "Number of flippings: " << arg.count_flipping() << endl;
        write_iterative_start();
        int placed = (int) arg.sample_nodes.size();
        if ((placed & (placed - 1)) == 0 or it == ordered_sample_nodes.end()) {
            sweep_samples();
        }
    }
    arg.resample_genotypes();
    cout << "orignal ARG length: " << arg.get_arg_length() << endl;
    arg.write_nodes(output_prefix + "_start_nodes_" + to_string(sample_index) + ".txt");
    arg.write_branches(output_prefix + "_start_branches_" + to_string(sample_index) + ".txt");
    arg.write_recombs(output_prefix + "_start_recombs_" + to_string(sample_index) + ".txt");
    if (!arg.genotype_likelihoods.empty()) {
        arg.write_cells(output_prefix + "_start_cells_" + to_string(sample_index) + ".txt");
    }
    arg.write_coordinates(output_prefix + "_coordinates.txt");
    write_mapping(output_prefix + "_start", to_string(sample_index), output_prefix + "_start_counts_" + to_string(sample_index) + ".txt", rescale(output_prefix + "_start_rescaled_nodes_" + to_string(sample_index) + ".txt"));
}

void Sampler::sweep_samples() {
    Threader_smc threader = Threader_smc(bsp_c, tsp_q);
    vector<Node *> nodes(arg.sample_nodes.begin(), arg.sample_nodes.end());
    for (Node *n : nodes) {
        random_engine.seed(random_seed);
        threader.reset();
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
            arg.penalty = penalty;
            arg.ancestral_prob = polar;
            if (exact) {
                threader.exact_internal_rethread(arg, arg.sample_uniform_cut());
            } else {
                threader.internal_rethread(arg, arg.sample_internal_cut());
            }
            moves += 1;
            arg.clear_remove_info();
        }
        arg.resample_genotypes();
        arg.release_dead_nodes();
        arg.write_nodes(output_prefix + "_nodes_" + to_string(sample_index) + ".txt");
        arg.write_branches(output_prefix + "_branches_" + to_string(sample_index) + ".txt");
        arg.write_recombs(output_prefix + "_recombs_" + to_string(sample_index) + ".txt");
        if (!arg.genotype_likelihoods.empty()) {
            arg.write_cells(output_prefix + "_cells_" + to_string(sample_index) + ".txt");
        }
        write_mapping(output_prefix, to_string(sample_index), output_prefix + "_counts_" + to_string(sample_index) + ".txt", rescale(output_prefix + "_rescaled_nodes_" + to_string(sample_index) + ".txt"));
        random_seed = random_engine();
        write_sample();
        arg.check_incompatibility();
        cout << "Start: " << arg.start << " , End: " << arg.end << endl;
        sample_index += 1;
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
    arg.penalty = penalty;
    arg.ancestral_prob = polar;
    arg.draw_mutations();
    arg.check_incompatibility();
    cout << "Number of trees: " << arg.recombinations.size() << endl;
    cout << "Number of flippings: " << arg.count_flipping() << endl;
    internal_sample(num_iters, spacing);
}


vector<double> Sampler::rescale(string filename) {
    vector<double> times;
    if (scaling_rep == 0) {
        return times;
    }
    Scaler scaler = Scaler();
    scaler.num_windows = scaling_bin;
    scaler.compute_deltas(arg, mut_map, start, mut_rate/Ne);
    for (int i = 0; i < scaling_rep; i++) {
        scaler.reset();
        scaler.rescale(arg, mut_rate);
    }
    times.resize(scaler.times.size());
    for (int i = 0; i < (int) times.size(); i++) {
        times[scaler.sorted_nodes[i]->index] = scaler.times[i];
    }
    ofstream file(filename);
    for (double t : times) {
        file << setprecision(numeric_limits<double>::max_digits10) << (t + arg.time_offset)*arg.Ne << "\n";
    }
    return times;
}

void Sampler::write_mapping(string base, string index, string counts_file, const vector<double> &rescaled) {
    function<double(Node *)> time_of = [&](Node *n) { return n == arg.root.get() ? n->time : rescaled[n->index]; };
    const function<double(Node *)> *tp = rescaled.empty() ? nullptr : &time_of;
    for (int r = 1; r <= mapping_rep; r++) {
        string file = base + "_muts_" + (mapping_rep > 1 ? to_string(r) + "_" : "") + index + ".txt";
        arg.draw_mutations(tp, arg.mutations);
        arg.write_mutations(file, arg.mutations);
    }
    vector<double> expected, informed_bases;
    arg.mapping_counts(tp, expected, informed_bases);
    arg.write_counts(counts_file, expected, informed_bases);
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
    arg.num_samples = ploidy*(int) data.sample_names().size();
    arg.read(output_prefix + "_nodes_" + to_string(sample_index) + ".txt", output_prefix + "_branches_" + to_string(sample_index) + ".txt", output_prefix + "_recombs_" + to_string(sample_index) + ".txt");
    arg.read_coordinates(output_prefix + "_coordinates.txt");
    vector<Node *> leaves(arg.sample_nodes.size());
    for (Node *n : arg.sample_nodes) {
        leaves[n->index] = n;
    }
    data.ploidy = ploidy;
    data.read(start, end);
    load_states(leaves);
    for (Node *n : leaves) {
        for (auto &x : n->mutation_sites) {
            arg.mutation_sites.insert(x.first);
        }
    }
    set_missing(leaves);
    arg.any_missing = data.any_missing;
    arg.unassayed_sites = data.unassayed_sites;
    arg.compute_rhos_thetas(recomb_rate, mut_rate, recomb_map, mut_map, start);
    arg.masked = data.masked;
    arg.compute_assayed();
    set_entries(leaves);
    if (!arg.genotype_likelihoods.empty()) {
        arg.read_cells(output_prefix + "_cells_" + to_string(sample_index) + ".txt");
    }
}

void Sampler::read_resume_point(string filename) {
    vector<string> words = Data_reader::read_last_line(filename);
    int log_length = (int) words.size();
    TSP::counter = stoi(words[log_length - 1]);
    sample_index = stoi(words[1]);
    load_resume_arg();
    arg.end = stod(words[log_length - 3]);
    arg.end_tree = arg.get_tree_at(arg.end);
}

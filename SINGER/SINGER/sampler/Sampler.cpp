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

Sampler::Sampler(double pop_size, Rate_map &rm, Rate_map &mm) {
    Ne = pop_size;
    mut_rate = mm.mean_rate()*Ne;
    recomb_rate = rm.mean_rate()*Ne;
    recomb_map = rm;
    mut_map = mm;
}

void Sampler::set_precision(double c, double q) {
    bsp_c = c;
    tsp_q = q;
}

void Sampler::set_input_file_prefix(string f) {
    input_prefix = f;
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

vector<string> Sampler::sample_names(string prefix) {
    ifstream file(prefix + ".vcf");
    string line;
    vector<string> names;
    while (getline(file, line)) {
        if (line.substr(0, 6) == "#CHROM") {
            istringstream iss(line);
            string field;
            int i = 0;
            while (iss >> field) {
                if (i >= 9) {
                    names.push_back(field);
                }
                i += 1;
            }
            return names;
        }
    }
    return names;
}

void Sampler::read_tip_ages(string filename, double g) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "tip ages file not found: " << filename << endl;
        exit(1);
    }
    vector<string> names = sample_names(input_prefix);
    int n = ploidy*(int) names.size();
    vector<pair<string, double>> rows;
    int ncol = 0;
    string line;
    while (getline(fin, line)) {
        istringstream iss(line);
        vector<string> parts;
        string tok;
        while (iss >> tok) {
            parts.push_back(tok);
        }
        if (parts.empty() or parts[0][0] == '#') {
            continue;
        }
        if (ncol == 0) {
            ncol = (int) parts.size();
        }
        if ((int) parts.size() != ncol or ncol > 2) {
            cerr << "tip ages file: expected 1 column (ages in VCF sample order) or 2 columns (name age)" << endl;
            exit(1);
        }
        rows.push_back({ncol == 2 ? parts[0] : "", stod(parts.back())});
    }
    if (rows.size() != names.size()) {
        cerr << "tip ages file has " << rows.size() << " rows for " << names.size() << " samples in the VCF" << endl;
        exit(1);
    }
    tip_times.assign(n, 0.0);
    if (ncol == 1) {
        for (int i = 0; i < n; i++) {
            tip_times[i] = rows[i/ploidy].second;
        }
    } else {
        map<string, double> by_name;
        for (auto &r : rows) {
            if (by_name.count(r.first) > 0) {
                cerr << "tip ages file: duplicate name " << r.first << endl;
                exit(1);
            }
            by_name[r.first] = r.second;
        }
        for (int i = 0; i < n; i++) {
            auto it = by_name.find(names[i/ploidy]);
            if (it == by_name.end()) {
                cerr << "tip ages file has no row for " << names[i/ploidy] << endl;
                exit(1);
            }
            tip_times[i] = it->second;
        }
    }
    for (double &x : tip_times) {
        x /= g*Ne;
    }
}

int Sampler::parse_genotype(const string &field, int expected_ploidy, int *calls) {
    size_t stop = field.find(':');
    if (stop == string::npos) {
        stop = field.size();
    }
    int n = 0;
    size_t i = 0;
    while (i < stop and n < expected_ploidy) {
        size_t j = i;
        while (j < stop and field[j] != '|' and field[j] != '/') {
            j++;
        }
        if (j == i + 1 and field[i] == '0') {
            calls[n] = 0;
        } else if (j == i + 1 and field[i] == '1') {
            calls[n] = 1;
        } else {
            calls[n] = -1;
        }
        n++;
        i = j + 1;
    }
    return n;
}

void Sampler::scan_missing(string prefix, double start_pos, double end_pos, vector<Node *> &leaves, int ploidy) {
    ifstream file(prefix + ".vcf");
    string line;
    long long prev_pos = -1;
    int calls[2];
    while (getline(file, line)) {
        if (line[0] == '#') {
            continue;
        }
        istringstream iss(line);
        string chrom, id, ref, alt, qual, filter, info, format, genotype;
        long long pos;
        iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format;
        if (pos < start_pos) {continue;}
        if (pos > end_pos) {break;}
        if (pos == prev_pos) {continue;}
        if (ref.size() > 1 or alt.size() > 1 or alt == "*" or alt == ".") {
            unassayed_site_list.push_back(pos - start_pos);
            continue;
        }
        streampos old_pos = file.tellg();
        string next_line;
        if (getline(file, next_line)) {
            istringstream next_iss(next_line);
            string next_chrom;
            long long next_pos;
            next_iss >> next_chrom >> next_pos;
            if (next_pos == pos and next_chrom == chrom) {
                unassayed_site_list.push_back(pos - start_pos);
                prev_pos = pos;
                continue;
            }
            file.seekg(old_pos);
        }
        int individual = 0;
        int site_missing = 0;
        while (iss >> genotype) {
            int n = parse_genotype(genotype, ploidy, calls);
            for (int k = 0; k < ploidy; k++) {
                if ((k < n ? calls[k] : -1) < 0) {
                    site_missing += 1;
                    leaves[ploidy*individual + k]->add_missing(pos - start_pos);
                }
            }
            individual += 1;
        }
        if (site_missing > 0) {
            any_missing = true;
        }
        if (site_missing == (int) leaves.size()) {
            unassayed_site_list.push_back(pos - start_pos);
        }
    }
    sort(unassayed_site_list.begin(), unassayed_site_list.end());
}

void Sampler::naive_read_vcf_haploid(string prefix, double start_pos, double end_pos) {
    string vcf_file = prefix + ".vcf";
    ifstream file(vcf_file);
    string line;
    int num_individuals = 0;
    long long prev_pos = -1;
    vector<Node_ptr> nodes = {};
    int valid_mutation = 0;
    int removed_mutation = 0;
    vector<double> genotypes = {};
    int calls[2];
    while (getline(file, line)) {
        if (line.substr(0, 6) == "#CHROM") {
            istringstream iss(line);
            vector<string> fields;
            string field;
            while (iss >> field) {
                fields.push_back(field);
            }
            num_individuals = (int) fields.size() - 9;
            nodes.resize(num_individuals);
            for (int i = 0; i < num_individuals; i++) {
                nodes[i] = new_sample(i);
                sample_nodes.insert(nodes[i]);
            }
            genotypes.resize(num_individuals);
            continue;
        } else if (line[0] == '#') {
            continue;        }
        istringstream iss(line);
        string chrom, id, ref, alt, qual, filter, info, format, genotype;
        long long pos;
        iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format;

        if (pos < start_pos) {continue;}
        if (pos > end_pos) {break;}
        if (pos == prev_pos) {continue;}        if (ref.size() > 1 or alt.size() > 1 or alt == "*" or alt == ".") {
            removed_mutation += 1;
            continue;
        }
        streampos old_pos = file.tellg();
        string next_line;
        if (getline(file, next_line)) {
            istringstream next_iss(next_line);
            string next_chrom;
            long long next_pos;
            next_iss >> next_chrom >> next_pos;
            if (next_pos == pos and next_chrom == chrom) {
                removed_mutation += 1;
                prev_pos = pos;
                continue;
            }
            file.seekg(old_pos);
        }
        int individual_index = 0;
        while (iss >> genotype) {
            int n = parse_genotype(genotype, 1, calls);
            genotypes[individual_index] = (n > 0 and calls[0] == 1) ? 1 : 0;
            individual_index += 1;
        }
        int genotype_sum = accumulate(genotypes.begin(), genotypes.end(), 0.0);
        if (genotype_sum >= 1 and genotype_sum < genotypes.size()) {
            valid_mutation += 1;
            for (int i = 0; i < genotypes.size(); i++) {
                if (genotypes[i] == 1) {
                    nodes[i]->add_mutation(pos - start_pos);
                }
            }
        }
    }
    ordered_sample_nodes = vector<Node_ptr>(sample_nodes.begin(), sample_nodes.end());
    shuffle(ordered_sample_nodes.begin(), ordered_sample_nodes.end(), random_engine);
    sequence_length = end_pos - start_pos;
    cout << "valid mutations: " << valid_mutation << endl;
    cout << "removed mutations: " << removed_mutation << endl;
    vector<Node *> leaves;
    for (const Node_ptr &n : sample_nodes) {
        leaves.push_back(n.get());
    }
    scan_missing(prefix, start_pos, end_pos, leaves, 1);
}

void Sampler::naive_read_vcf(string prefix, double start_pos, double end_pos) {
    string vcf_file = prefix + ".vcf";
    ifstream file(vcf_file);
    string line;
    int num_individuals = 0;
    long long prev_pos = -1;
    vector<Node_ptr> nodes = {};
    int valid_mutation = 0;
    int removed_mutation = 0;
    vector<double> genotypes = {};
    int calls[2];
    while (getline(file, line)) {
        if (line.substr(0, 6) == "#CHROM") {
            istringstream iss(line);
            vector<string> fields;
            string field;
            while (iss >> field) {
                fields.push_back(field);
            }
            num_individuals = (int) fields.size() - 9;
            nodes.resize(2*num_individuals);
            for (int i = 0; i < 2*num_individuals; i++) {
                nodes[i] = new_sample(i);
                sample_nodes.insert(nodes[i]);
            }
            genotypes.resize(2*num_individuals);
            continue;
        } else if (line[0] == '#') {
            continue;        }
        istringstream iss(line);
        string chrom, id, ref, alt, qual, filter, info, format, genotype;
        long long pos;
        iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format;

        if (pos < start_pos) {continue;}
        if (pos > end_pos) {break;}
        if (pos == prev_pos) {continue;}        if (ref.size() > 1 or alt.size() > 1 or alt == "*" or alt == ".") {
            removed_mutation += 1;
            continue;
        }
        streampos old_pos = file.tellg();
        string next_line;
        if (getline(file, next_line)) {
            istringstream next_iss(next_line);
            string next_chrom;
            long long next_pos;
            next_iss >> next_chrom >> next_pos;
            if (next_pos == pos and next_chrom == chrom) {
                removed_mutation += 1;
                prev_pos = pos;
                continue;
            }
            file.seekg(old_pos);
        }
        int individual_index = 0;
        while (iss >> genotype) {
            int n = parse_genotype(genotype, 2, calls);
            for (int k = 0; k < 2; k++) {
                genotypes[2*individual_index + k] = (k < n and calls[k] == 1) ? 1 : 0;
            }
            individual_index += 1;
        }
        int genotype_sum = accumulate(genotypes.begin(), genotypes.end(), 0.0);
        if (genotype_sum >= 1 and genotype_sum < genotypes.size()) {
            valid_mutation += 1;
            for (int i = 0; i < genotypes.size(); i++) {
                if (genotypes[i] == 1) {
                    nodes[i]->add_mutation(pos - start_pos);
                }
            }
        }
    }
    ordered_sample_nodes = vector<Node_ptr>(sample_nodes.begin(), sample_nodes.end());
    shuffle(ordered_sample_nodes.begin(), ordered_sample_nodes.end(), random_engine);
    sequence_length = end_pos - start_pos;
    cout << "valid mutations: " << valid_mutation << endl;
    cout << "removed mutations: " << removed_mutation << endl;
}

void Sampler::guide_read_vcf(string prefix, double start, double end) {
    random_engine.seed(random_seed);
    string vcf_file = prefix + ".vcf";
    string index_file = prefix + ".index";
    ifstream idx_stream(index_file);
    if (!idx_stream.is_open()) {
        cerr << "Index file not found: " + index_file << endl;
        exit(1);
    }
    string line;
    long byte_offset = -1;
    while (getline(idx_stream, line)) {
        istringstream iss(line);
        double segment_start;
        long offset;
        iss >> segment_start >> offset;
        if (segment_start == start) {
            byte_offset = offset;
            break;
        }
    }
    if (byte_offset == -1) {
        cerr << "Start position not found in index file: " + index_file << endl;
        exit(1);
    }
    ifstream vcf_stream(vcf_file, ios::binary);
    if (!vcf_stream.is_open()) {
        cerr << "VCF file not found: " + vcf_file << endl;
    }
    vcf_stream.seekg(byte_offset, ios::beg);
    long long prev_pos = -1;
    vector<Node_ptr> nodes = {};
    int valid_mutation = 0;
    int removed_mutation = 0;
    vector<double> genotypes = {};
    int calls[2];
    while (getline(vcf_stream, line)) {
        istringstream iss(line);
        string chrom, id, ref, alt, qual, filter, info, format, genotype;
        long long pos;
        iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format;
        if (pos == prev_pos) {continue;}        if (pos >= end) {break;}        if (ref.size() > 1 or alt.size() > 1 or alt == "*" or alt == ".") {
            removed_mutation += 1;
            continue;
        }        streampos old_pos = vcf_stream.tellg();
        string next_line;
        if (getline(vcf_stream, next_line)) {
            istringstream next_iss(next_line);
            string next_chrom;
            long long next_pos;
            next_iss >> next_chrom >> next_pos;
            if (next_pos == pos and next_chrom == chrom) {
                removed_mutation += 1;
                prev_pos = pos;
                continue;
            }
            vcf_stream.seekg(old_pos);
        }
        int individual_index = 0;
        while (iss >> genotype) {
            if (genotypes.size() < 2*individual_index + 2) {
                genotypes.resize(2*individual_index + 2);
            }
            int n = parse_genotype(genotype, 2, calls);
            for (int k = 0; k < 2; k++) {
                genotypes[2*individual_index + k] = (k < n and calls[k] == 1) ? 1 : 0;
            }
            individual_index += 1;
        }
        if (nodes.size() == 0) {
            nodes.resize(genotypes.size());
            for (int i = 0; i < nodes.size(); i++) {
                nodes[i] = new_sample(i);
                sample_nodes.insert(nodes[i]);
            }
        } else {
            assert(nodes.size() == genotypes.size());
        }
        int genotype_sum = accumulate(genotypes.begin(), genotypes.end(), 0.0);
        if (genotype_sum >= 1 and genotype_sum < genotypes.size()) {
            valid_mutation += 1;
            for (int i = 0; i < genotypes.size(); i++) {
                if (genotypes[i] == 1) {
                    nodes[i]->add_mutation(pos - start);
                }
            }
        }
    }
    if (valid_mutation < 3) {
        cerr << "there are too few variants in this region, algorithm not run" << endl;
    }
    ordered_sample_nodes = vector<Node_ptr>(sample_nodes.begin(), sample_nodes.end());
    sequence_length = end - start;
    cout << "valid mutations: " << valid_mutation << endl;
    cout << "removed mutations: " << removed_mutation << endl;
}

void Sampler::load_vcf(string prefix, double start, double end) {
    string index_file = prefix + ".index";
    ifstream idx_stream(index_file);
    if (idx_stream.is_open()) {
        guide_read_vcf(prefix, start, end);
    } else {
        naive_read_vcf(prefix, start, end);
    }
    vector<Node *> leaves(sample_nodes.size());
    for (const Node_ptr &n : sample_nodes) {
        leaves[n->index] = n.get();
    }
    scan_missing(prefix, start, end, leaves, 2);
}

void Sampler::build_singleton_arg() {
    double bin_size = max(1.0, rho_unit/recomb_rate);
    bin_size = min(bin_size, 100.0);
    Node_ptr n = *ordered_sample_nodes.begin();
    arg = ARG(Ne, sequence_length);
    arg.any_missing = any_missing;
    arg.unassayed_sites = unassayed_site_list;
    arg.discretize(bin_size);
    arg.build_singleton_arg(n);
    if (mut_rate > 0 and recomb_rate > 0) {
        arg.compute_rhos_thetas(recomb_rate, mut_rate);
    } else {
        arg.compute_rhos_thetas(recomb_map, mut_map);
    }
    arg.discount_unassayed();
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
    }
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
    sample_index += 1;
    arg.check_incompatibility();
    cout << "Number of trees: " << arg.recombinations.size() << endl;
    cout << "Number of flippings: " << arg.count_flipping() << endl;
    internal_sample(num_iters, spacing);
}

void Sampler::debug_resume_internal_sample(int num_iters, int spacing) {
    retract_log(5);
    string log_file = output_prefix + ".log";
    vector<string> words = read_last_line(log_file);
    if (words.size() == 0 or words[2] == "initial_thread" or words[0] == "Time") {
        cout << "new seed: " << random_seed << endl;
        sample_index = 0;
        load_vcf(input_prefix, start, end);
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
    scaler.compute_deltas(arg);
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
    arg.num_samples = ploidy*(int) sample_names(input_prefix).size();
    arg.read(node_file, branch_file, recomb_file, mut_file);
    arg.read_coordinates(coord_file);
    vector<Node *> leaves(arg.sample_nodes.size());
    for (Node *n : arg.sample_nodes) {
        leaves[n->index] = n;
    }
    scan_missing(input_prefix, start, end, leaves, 2);
    arg.any_missing = any_missing;
    arg.unassayed_sites = unassayed_site_list;
    if (mut_rate > 0 and recomb_rate > 0) {
        arg.compute_rhos_thetas(recomb_rate, mut_rate);
    } else {
        arg.compute_rhos_thetas(recomb_map, mut_map);
    }
    arg.discount_unassayed();
}

vector<string> Sampler::read_last_line(string filename) {
    ifstream file(filename, ios::in);
    vector<string> words;
    if (!file.is_open()) {
        cerr << "Error opening the file: " << filename << endl;
        exit(1);
    }
    file.seekg(-1, ios_base::end);
    bool found = false;
    int newlines_encountered = 0;
    char ch;
    while (file.get(ch)) {
        if (ch == '\n') {
            newlines_encountered++;
            if (newlines_encountered == 2) {
                found = true;
                break;
            }
        }

        file.seekg(-2, ios_base::cur);
    }

    if (!found) {
        return words;
    }

    string last_line_with_content;
    getline(file, last_line_with_content);

    stringstream ss(last_line_with_content);
    string word;
    while (ss >> word) {
        words.push_back(word);
    }
    return words;
}

void Sampler::read_resume_point(string filename) {
    vector<string> words = read_last_line(filename);
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

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
        write_log((int) arg.sample_nodes.size(), "initial_thread");
        int placed = (int) arg.sample_nodes.size();
        if ((placed & (placed - 1)) == 0 or it == ordered_sample_nodes.end()) {
            sweep_samples();
        }
    }
    arg.resample_genotypes();
    cout << "orignal ARG length: " << arg.get_arg_length() << endl;
    arg.write_coordinates(output_prefix + "_coordinates.txt");
    write_output(output_prefix + "_start_" + to_string(sample_index) + ".bin");
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
    write_log((int) arg.sample_nodes.size(), "sweep");
}

template <class T> static void put(ofstream &f, T x) {
    f.write(reinterpret_cast<const char *>(&x), sizeof(T));
}

static double quantile(const vector<double> &v, double q) {
    double h = q*(v.size() - 1);
    int lo = (int) floor(h);
    int hi = min(lo + 1, (int) v.size() - 1);
    return v[lo] + (h - lo)*(v[hi] - v[lo]);
}

static void summarize(const Tree &tree, int n, double Ne, double offset, double *stats, vector<int> &history) {
    vector<Node *> inner;
    vector<pair<Node *, Node *>> edges;
    double tbl = 0;
    for (auto &p : tree.parents) {
        if (p.second->index != -1) {
            tbl += (p.second->time - p.first->time)*Ne;
            edges.push_back({p.second, p.first});
        }
        if (!p.first->is_sample) {
            inner.push_back(p.first);
        }
    }
    sort(inner.begin(), inner.end(), compare_node());
    sort(edges.begin(), edges.end());
    vector<double> ages(inner.size());
    for (int i = 0; i < (int) inner.size(); i++) {
        ages[i] = (inner[i]->time + offset)*Ne;
    }
    stats[0] = tbl;
    stats[1] = quantile(ages, 0.25);
    stats[2] = quantile(ages, 0.5);
    stats[3] = quantile(ages, 0.75);
    vector<Node *> active(n);
    for (auto &p : tree.parents) {
        if (p.first->is_sample) {
            active[p.first->index] = p.first;
        }
    }
    history.clear();
    for (Node *u : inner) {
        auto e = lower_bound(edges.begin(), edges.end(), make_pair(u, (Node *) nullptr));
        int a = (int) (find(active.begin(), active.end(), e->second) - active.begin());
        int b = (int) (find(active.begin(), active.end(), (e + 1)->second) - active.begin());
        if (a > b) {
            swap(a, b);
        }
        history.push_back(b*(b - 1)/2 + a);
        active.erase(active.begin() + b);
        active.erase(active.begin() + a);
        active.push_back(u);
    }
}

void Sampler::refresh_summary() {
    int n = (int) arg.sample_nodes.size();
    Tree tree;
    int at = -2;
    map<double, Recombination>::iterator it;
    for (int j = 0; j < (int) log_positions.size(); j++) {
        if (!log_dirty[j]) {
            continue;
        }
        if (at == j - 1) {
            while (it->first <= log_positions[j]) {
                tree.forward_update(it->second);
                it++;
            }
        } else {
            tree = arg.get_tree_at(log_positions[j]);
            it = arg.recombinations.upper_bound(log_positions[j]);
        }
        at = j;
        summarize(tree, n, arg.Ne, arg.time_offset, &log_stats[4*j], log_histories[j]);
        log_dirty[j] = 0;
        log_changed[j] = 1;
    }
}

static void append_number(string &line, double v) {
    char buffer[32];
    auto r = to_chars(buffer, buffer + 32, v);
    line.append(buffer, r.ptr);
}

void Sampler::write_summary_row() {
    string line = to_string(moves_done) + "\t" + to_string(arg.recombinations.size() - 2);
    for (double v : log_stats) {
        line += '\t';
        append_number(line, v);
    }
    summary_file << line << '\n';
}

void Sampler::write_trees_row() {
    put(trees_file, (int64_t) moves_done);
    put(trees_file, (int32_t) count(log_changed.begin(), log_changed.end(), 1));
    for (int j = 0; j < (int) log_positions.size(); j++) {
        if (log_changed[j]) {
            put(trees_file, (int32_t) j);
            for (int c : log_histories[j]) {
                if (code_width == 2) {
                    put(trees_file, (uint16_t) c);
                } else {
                    put(trees_file, (uint32_t) c);
                }
            }
            log_changed[j] = 0;
        }
    }
}

void Sampler::start_summary(bool resume) {
    if (thin_log == 0) {
        return;
    }
    int n = (int) arg.sample_nodes.size();
    code_width = n*(n - 1)/2 - 1 < 65536 ? 2 : 4;
    log_positions.clear();
    for (int j = 0; j < thin_log_pos; j++) {
        log_positions.push_back((j + 0.5)*arg.sequence_length/thin_log_pos);
    }
    log_stats.assign(4*thin_log_pos, 0);
    log_histories.assign(thin_log_pos, {});
    log_dirty.assign(thin_log_pos, 1);
    log_changed.assign(thin_log_pos, 0);
    refresh_summary();
    string summary_name = output_prefix + "_summary.log", trees_name = output_prefix + "_summary_trees.bin";
    if (resume) {
        ifstream text(summary_name);
        string line;
        long long offset = 0;
        while (getline(text, line) and (line[0] == '#' or line[0] == 's' or stoll(line) <= moves_done)) {
            offset += line.size() + 1;
        }
        text.close();
        filesystem::resize_file(summary_name, offset);
        ifstream in(trees_name, ios::binary);
        offset = 4*sizeof(int64_t) + (1 + thin_log_pos)*sizeof(double);
        long long entry = sizeof(int32_t) + (n - 1)*code_width;
        in.seekg(offset);
        int64_t iteration;
        while (in.read(reinterpret_cast<char *>(&iteration), sizeof(int64_t)) and iteration <= moves_done) {
            int32_t m;
            in.read(reinterpret_cast<char *>(&m), sizeof(int32_t));
            offset += sizeof(int64_t) + sizeof(int32_t) + m*entry;
            in.seekg(offset);
        }
        in.close();
        filesystem::resize_file(trees_name, offset);
        summary_file.open(summary_name, ios::app);
        trees_file.open(trees_name, ios::binary | ios::app);
        fill(log_changed.begin(), log_changed.end(), 0);
    } else {
        summary_file.open(summary_name, ios::trunc);
        string line = "# positions:";
        for (double x : log_positions) {
            line += ' ';
            append_number(line, x);
        }
        line += "\nstate\trecombinations";
        for (int j = 0; j < thin_log_pos; j++) {
            for (string s : {"tbl_", "q25_", "q50_", "q75_"}) {
                line += "\t" + s + to_string(j);
            }
        }
        summary_file << line << '\n';
        write_summary_row();
        summary_file.flush();
        trees_file.open(trees_name, ios::binary | ios::trunc);
        for (int64_t h : {(int64_t) thin_log_pos, (int64_t) n, (int64_t) thin_log, (int64_t) code_width}) {
            put(trees_file, h);
        }
        put(trees_file, arg.sequence_length);
        for (double x : log_positions) {
            put(trees_file, x);
        }
        write_trees_row();
        trees_file.flush();
    }
}

void Sampler::run_moves(Threader_smc &threader, int num_moves) {
    for (int moves = 0; moves < num_moves; moves++) {
        threader.reset();
        arg.penalty = penalty;
        arg.ancestral_prob = polar;
        threader.exact_internal_rethread(arg, arg.sample_uniform_cut());
        arg.clear_remove_info();
        moves_done += 1;
        if (thin_log > 0) {
            if (threader.accepted) {
                for (auto j = lower_bound(log_positions.begin(), log_positions.end(), arg.start); j != log_positions.end() and *j <= arg.end; j++) {
                    log_dirty[j - log_positions.begin()] = 1;
                }
            }
            if (moves_done % thin_log == 0) {
                refresh_summary();
                write_summary_row();
                write_trees_row();
            }
        }
    }
    arg.release_dead_nodes();
}

void Sampler::burn_in(int num_moves, int spacing) {
    if (num_moves == 0) {
        return;
    }
    Threader_smc threader = Threader_smc(bsp_c, tsp_q);
    random_engine.seed(random_seed);
    for (int done = 0; done < num_moves;) {
        int k = min(spacing, num_moves - done);
        run_moves(threader, k);
        done += k;
        write_log(moves_done, "burn_in");
    }
    random_seed = random_engine();
}

void Sampler::internal_sample(int num_iters, int spacing) {
    Threader_smc threader = Threader_smc(bsp_c, tsp_q);
    while (sample_index < num_iters) {
        int seed = random_seed;
        random_engine.seed(random_seed);
        run_moves(threader, spacing);
        write_output(output_prefix + "_" + to_string(sample_index) + ".bin");
        if (thin_log > 0) {
            summary_file.flush();
            trees_file.flush();
        }
        random_seed = random_engine();
        write_log(moves_done, "main_run");
        cout << get_time() << " Sample: " << sample_index << ", Iteration: " << moves_done << ", Random seed: " << seed << ", Start: " << arg.start << ", End: " << arg.end << ", Number of recombinations: " << arg.recombinations.size() - 2 << endl;
        sample_index += 1;
    }
}

void Sampler::resume_internal_sample(int num_iters, int spacing, int burn) {
    vector<string> words;
    ifstream log(output_prefix + ".log");
    string line;
    while (getline(log, line)) {
        istringstream iss(line);
        vector<string> row;
        string w;
        while (iss >> w) {
            row.push_back(w);
        }
        if (row.size() > 2 and (row[2] == "main_run" or row[2] == "sweep")) {
            words = row;
        }
    }
    if (words.empty()) {
        cerr << "Error: " << output_prefix << ".log has no sweep or main_run row to resume from; start the run again without -resume. " << endl;
        exit(1);
    }
    bool sampled = words[2] == "main_run";
    moves_done = sampled ? stoll(words[1]) : 0;
    int done_samples = (int) ((moves_done - burn)/spacing);
    read_resume_point(words, output_prefix + (sampled ? "_" + to_string(done_samples - 1) : "_start_0") + ".bin");
    if (!seed_given) {
        random_seed = stoi(words[words.size() - 2]);
    }
    arg.penalty = penalty;
    arg.ancestral_prob = polar;
    cout << "Number of recombinations: " << arg.recombinations.size() - 2 << endl;
    start_summary(true);
    if (sampled) {
        sample_index = done_samples;
    } else {
        sample_index = 0;
        burn_in(burn, spacing);
    }
    internal_sample(num_iters, spacing);
}


vector<double> Sampler::rescale() {
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
    return times;
}

void Sampler::write_output(string filename) {
    arg.number_nodes();
    vector<Cell_row> cells = arg.cell_rows();
    vector<double> rescaled = rescale();
    function<double(Node *)> time_of = [&](Node *n) { return n == arg.root.get() ? n->time : rescaled[n->index]; };
    const function<double(Node *)> *tp = rescaled.empty() ? nullptr : &time_of;
    vector<map<double, vector<Mutation_row>>> draws;
    for (int r = 1; r <= mapping_rep; r++) {
        arg.draw_mutations(tp, arg.mutations);
        draws.push_back(arg.mutations);
    }
    vector<double> expected, informed_bases;
    arg.mapping_counts(tp, expected, informed_bases);
    arg.write_sample(filename, cells, rescaled, draws, expected, informed_bases);
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
    << "Last_updated_pos" << "\t"
    << "Random_seed" << "\t"
    << "Counter" << endl;
    file.close();
}

void Sampler::write_log(long long iteration, string type) {
    string filename = output_prefix + ".log";
    ofstream file(filename, ios::out|ios::app);
    if (!file) {
        cerr << "Error opening the file: " << filename << endl;
        return;
    }
    file << get_time() << "\t"
    << iteration << "\t"
    << type << "\t"
    << arg.recombinations.size() - 2 << "\t"
    << setprecision(numeric_limits<double>::max_digits10)
    << arg.end << "\t"
    << random_seed << "\t"
    << TSP::counter << endl;
}

void Sampler::load_resume_arg(string filename) {
    arg = ARG(Ne, sequence_length);
    arg.num_samples = ploidy*(int) data.sample_names().size();
    arg.read_sample(filename);
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
        arg.read_cells(filename);
    }
}

void Sampler::read_resume_point(const vector<string> &words, string filename) {
    int log_length = (int) words.size();
    TSP::counter = stoi(words[log_length - 1]);
    load_resume_arg(filename);
}

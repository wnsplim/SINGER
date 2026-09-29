//
//  ARG.cpp
//  SINGER
//
//  Created by Yun Deng on 4/14/22.
//

#include "ARG.hpp"

ARG::ARG() {}

ARG::ARG(double N, double l) {
    root->set_index(-1);
    Ne = N;
    sequence_length = l;
    end = l;
    Recombination r = Recombination({}, {});
    r.set_pos(0.0);
    recombinations[0] = r;
    r = Recombination({}, {});
    r.set_pos(INT_MAX);
    recombinations[INT_MAX] = r;
    mutation_sites.insert(INT_MAX);
    mutation_branches[INT_MAX] = {};
}

ARG::~ARG() {
}

void ARG::discretize(double s) {
    bin_size = s;
    auto recomb_it = recombinations.upper_bound(0);
    double curr_pos = 0;
    while (curr_pos < sequence_length) {
        coordinates.push_back(curr_pos);
        if (recomb_it->first < curr_pos + s) {
            curr_pos = recomb_it->first;
            recomb_it++;
        } else {
            curr_pos = min(curr_pos + s, sequence_length);
        }
    }
    coordinates.push_back(sequence_length);
    bin_num = (int) coordinates.size() - 1;
}

int ARG::get_index(double x) {
    auto it = upper_bound(coordinates.begin(), coordinates.end(), x);
    --it;
    int index = (int) distance(coordinates.begin(), it);
    return index;
}

void ARG::compute_rhos_thetas(double r, double m) {
    int n = (int) coordinates.size() - 1;
    for (int i = 0; i < n; i++) {
        rhos.push_back(r*(coordinates[i+1] - coordinates[i]));
        thetas.push_back(m*(coordinates[i+1] - coordinates[i]));
    }
}

void ARG::compute_rhos_thetas(Rate_map &rm, Rate_map &mm) {
    int n = (int) coordinates.size() - 1;
    for (int i = 0; i < n; i++) {
        rhos.push_back(rm.segment_distance(coordinates[i], coordinates[i+1])*Ne);
        thetas.push_back(mm.segment_distance(coordinates[i], coordinates[i+1])*Ne);
    }
}

void ARG::build_singleton_arg(const Node_ptr &n) {
    add_sample(n);
    Branch branch = Branch(n, root);
    Recombination r = Recombination({}, {branch});
    r.set_pos(0.0);
    recombinations[0] = r;
    for (double x : mutation_sites) {
        mutation_branches[x] = {branch};
    }
}

void ARG::add_sample(const Node_ptr &n) {
    node_owner.push_back(n);
    n->is_sample = true;
    sample_nodes.insert(n.get());
    for (auto &x : n->mutation_sites) {
        mutation_sites.insert(x.first);
    }
    removed_branches.clear();
    removed_branches[0] = Branch(n, root);
    removed_branches[sequence_length] = Branch();
    anchors_changed(0, INT_MAX);
    start_tree = get_tree_at(0);
    cut_pos = 0;
    start = 0;
    end = sequence_length;
}

void ARG::add_node(Node *n) {
    if (n != root.get() and n != nullptr) {
        node_set.insert(n);
    }
}

void ARG::add_new_node(double t, bool sample) {
    if (!isinf(t)) {
        Node_ptr n = new_node(t);
        n->index = (int) node_set.size();
        n->is_sample = sample;
        node_owner.push_back(n);
        node_set.insert(n.get());
        node_list.push_back(n.get());
        if (sample) {
            sample_nodes.insert(n.get());
        }
    }
}
 
void ARG::resort_after_time_change() {
    for (auto &x : recombinations) {
        Recombination &r = x.second;
        r.deleted_branches = set<Branch>(r.deleted_branches.begin(), r.deleted_branches.end());
        r.inserted_branches = set<Branch>(r.inserted_branches.begin(), r.inserted_branches.end());
    }
    for (auto &x : mutation_branches) {
        x.second = set<Branch>(x.second.begin(), x.second.end());
    }
    node_set.clear();
    create_node_set();
    anchors_changed(0, INT_MAX);
    start_tree = get_tree_at(start);
    end_tree = get_tree_at(end);
}

Tree ARG::get_tree_at(double x) {
    if (anchor_step == 0) {
        anchor_step = max(8, (int) sqrt((double) recombinations.size()));
        anchor_dirty.push_back({0, (double) INT_MAX});
    }
    if (!anchor_dirty.empty()) {
        rebuild_anchors();
    }
    Tree tree;
    auto recomb_it = recombinations.begin();
    auto a = anchors.upper_bound(x);
    if (a != anchors.begin()) {
        --a;
        tree = a->second;
        recomb_it = recombinations.upper_bound(a->first);
    }
    while (recomb_it->first <= x) {
        Recombination &r = recomb_it->second;
        tree.forward_update(r);
        recomb_it++;
    }
    return tree;
}

void ARG::anchors_changed(double lo, double hi) {
    anchors.erase(anchors.lower_bound(lo), anchors.upper_bound(hi));
    anchor_dirty.push_back({lo, hi});
}

void ARG::rebuild_anchors() {
    for (auto &d : anchor_dirty) {
        Tree tree;
        auto it = recombinations.begin();
        auto a = anchors.lower_bound(d.first);
        if (a != anchors.begin()) {
            --a;
            tree = a->second;
            it = recombinations.upper_bound(a->first);
        }
        int count = 0;
        while (it != recombinations.end() and it->first <= d.second) {
            tree.forward_update(it->second);
            count += 1;
            if (count % anchor_step == 0 and it->first >= d.first) {
                anchors[it->first] = tree;
            }
            ++it;
        }
    }
    anchor_dirty.clear();
}

void ARG::window_copy_into(ARG &c, double lo_pos) {
    c.Ne = Ne;
    c.root = root;
    c.cut_node = cut_node;
    c.cut_node_owner = cut_node_owner;
    c.cut_time = cut_time;
    c.anchors.clear();
    c.anchor_dirty.clear();
    c.anchor_step = anchor_step;
    c.mutation_sites.clear();
    for (auto it = mutation_sites.lower_bound(lo_pos); it != mutation_sites.end() and *it <= end; ++it) {
        c.mutation_sites.insert(c.mutation_sites.end(), *it);
    }
    c.mutation_sites.insert(c.mutation_sites.end(), *mutation_sites.rbegin());
    c.mutation_branches.clear();
    c.recombinations.clear();
    c.recombinations.insert(*recombinations.begin());
    for (auto it = recombinations.lower_bound(lo_pos); it != recombinations.end() and it->first <= end; ++it) {
        c.recombinations.insert(c.recombinations.end(), *it);
    }
    c.recombinations.insert(*recombinations.rbegin());
    c.bin_num = bin_num;
    c.penalty = penalty;
    c.ancestral_prob = ancestral_prob;
    c.any_missing = any_missing;
    c.unassayed_sites = unassayed_sites;
    c.assayed = assayed;
    c.sequence_length = sequence_length;
    c.bin_size = bin_size;
    if (c.coordinates.size() != coordinates.size()) {
        c.coordinates = coordinates;
        c.rhos = rhos;
        c.thetas = thetas;
    }
    c.joining_branches = joining_branches;
    c.removed_branches = removed_branches;
    c.start = start;
    c.end = end;
    c.cut_pos = cut_pos;
    c.cut_tree = cut_tree;
    c.start_tree = start_tree;
    c.end_tree = end_tree;
}

Node *ARG::get_query_node_at(double x) {
    auto query_it = removed_branches.upper_bound(x);
    query_it--;
    return query_it->second.lower_node;
}

void ARG::remove(tuple<double, Branch, double> cut_point) {
    double pos;
    Branch center_branch;
    double t;
    tie(pos, center_branch, t) = cut_point;
    cut_time = t;
    cut_node_owner = new_node(cut_time);
    cut_node = cut_node_owner.get();
    cut_node->set_index(-2);
    Tree forward_tree = cut_tree;
    Tree backward_tree = cut_tree;
    auto f_it = recombinations.upper_bound(pos);
    auto b_it = recombinations.upper_bound(pos);
    Branch prev_joining_branch;
    Branch next_joining_branch;
    Branch prev_removed_branch = center_branch;
    Branch next_removed_branch = center_branch;
    while (next_removed_branch != Branch()) {
        Recombination &r = f_it->second;
        prev_joining_branch = forward_tree.find_joining_branch(prev_removed_branch);
        forward_tree.forward_update(r);
        next_removed_branch = r.trace_forward(t, prev_removed_branch);
        if (next_removed_branch.upper_node == root.get()) {
            next_removed_branch = Branch();
        }
        next_joining_branch = forward_tree.find_joining_branch(next_removed_branch);
        r.remove(prev_removed_branch, next_removed_branch, prev_joining_branch, next_joining_branch, cut_node);
        removed_branches[min(r.pos, sequence_length)] = next_removed_branch;
        joining_branches[min(r.pos, sequence_length)] = next_joining_branch;
        f_it++;
        prev_removed_branch = next_removed_branch;
    }
    next_removed_branch = center_branch;
    prev_removed_branch = center_branch;
    while (prev_removed_branch != Branch()) {
        b_it--;
        Recombination &r = b_it->second;
        removed_branches[r.pos] = prev_removed_branch;
        next_joining_branch = backward_tree.find_joining_branch(next_removed_branch);
        joining_branches[r.pos] = next_joining_branch;
        backward_tree.backward_update(r);
        prev_removed_branch = r.trace_backward(t, next_removed_branch);
        if (prev_removed_branch.upper_node == root.get()) {
            prev_removed_branch = Branch();
        }
        if (prev_removed_branch == Branch()) {
            backward_tree.forward_update(r);
        }
        prev_joining_branch = backward_tree.find_joining_branch(prev_removed_branch);
        r.remove(prev_removed_branch, next_removed_branch, prev_joining_branch, next_joining_branch, cut_node);
        next_removed_branch = prev_removed_branch;
    }
    start = removed_branches.begin()->first;
    end = removed_branches.rbegin()->first;
    remove_empty_recombinations();
    anchors_changed(start, end);
    remap_mutations();
    cut_tree.remove(center_branch, cut_node);
    backward_tree.remove(removed_branches.begin()->second, cut_node);
    start_tree = move(backward_tree);
    end_tree = move(forward_tree);
}

void ARG::add(map<double, Branch> &new_joining_branches, map<double, Branch> &added_branches) {
    auto join_it = new_joining_branches.begin();
    auto add_it = added_branches.begin();
    auto recomb_it = recombinations.lower_bound(start);
    Branch prev_joining_branch = Branch();
    Branch next_joining_branch = Branch();
    Branch prev_added_branch = Branch();
    Branch next_added_branch = Branch();
    while (add_it != added_branches.end() and add_it->first < sequence_length) {
        assert(add_it->first <= join_it->first);
        if (recomb_it->first == add_it->first) {
            Recombination &r = recomb_it->second;
            recomb_it++;
            if (join_it->first == add_it->first) {
                next_joining_branch = join_it->second;
                join_it++;
            }
            next_added_branch = add_it->second;
            add_it++;
            r.add(prev_added_branch, next_added_branch, prev_joining_branch, next_joining_branch, cut_node);
            prev_joining_branch = next_joining_branch;
            prev_added_branch = next_added_branch;
        } else {
            if (join_it->first == add_it->first) {
                next_joining_branch = join_it->second;
                join_it++;
            }
            next_added_branch = add_it->second;
            new_recombination(add_it->first, prev_added_branch, prev_joining_branch, next_added_branch, next_joining_branch);
            add_it++;
            prev_joining_branch = next_joining_branch;
            prev_added_branch = next_added_branch;
        }
    }
    remove_empty_recombinations();
    anchors_changed(start, end);
    impute(new_joining_branches, added_branches);
    start_tree.add(added_branches.begin()->second, new_joining_branches.begin()->second, cut_node);
}

Branch ARG::lineage_branch_before(map<double, Branch> &lineage, double x) {
    auto it = lineage.lower_bound(x);
    if (it == lineage.begin()) {
        return Branch();
    }
    return prev(it)->second;
}

void ARG::smc_sample_recombinations(map<double, Branch> &lineage) {
    RSP_smc rsp = RSP_smc();
    Tree tree = start_tree;
    auto it = recombinations.upper_bound(start);
    while (it->first < end) {
        Recombination &r = it->second;
        if (r.pos != 0 and r.pos < sequence_length) {
            rsp.sample_recombination(r, cut_time, tree, lineage_branch_before(lineage, r.pos));
            assert(r.start_time > 0);
        }
        tree.forward_update(r);
        it++;
    }
}

void ARG::approx_sample_recombinations() {
    RSP_smc rsp = RSP_smc();
    auto it = recombinations.lower_bound(start);
    while (it->first <= end and it->first < sequence_length) {
        Recombination &r = it->second;
        if (r.pos > 0 and r.pos < sequence_length) {
            rsp.approx_sample_recombination(r, cut_time);
            assert(r.start_time > 0);
            assert(r.start_time <= r.inserted_node->time);
            assert(r.start_time <= r.deleted_node->time);
        }
        it++;
    }
}

void ARG::adjust_recombinations() {
    RSP_smc rsp = RSP_smc();
    auto it = recombinations.upper_bound(0);
    while (it->first < sequence_length) {
        Recombination &r = it->second;
        if (r.pos != 0 and r.pos < sequence_length) {
            rsp.adjust(r, 0);
            assert(r.start_time > 0);
            assert(r.start_time <= r.inserted_node->time);
            assert(r.start_time <= r.deleted_node->time);
        }
        it++;
    }
}

int ARG::count_flipping() {
    int count = 0;
    for (auto &x : mutation_branches) {
        set<Branch> &branches = x.second;
        if (branches.size() > 1 and branches.rbegin()->upper_node == root.get()) {
            count += 1;
        }
    }
    return count;
}

void ARG::read_coordinates(string filename) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    double x;
    while (fin >> x) {
        coordinates.push_back(x);
    }
    bin_num = (int) coordinates.size() - 1;
    return;
}

void ARG::write_coordinates(string filename) {
    ofstream file;
    file.open(filename);
    file << std::setprecision(std::numeric_limits<double>::max_digits10) << std::fixed;
    for (auto x : coordinates) {
        file << x << "\n";
    }
    return;
}

void ARG::write(string node_file, string branch_file, string recomb_file, string mutation_file) {
    write_nodes(node_file);
    write_branches(branch_file);
    write_recombs(recomb_file);
    write_mutations(mutation_file);
}

void ARG::read(string node_file, string branch_file, string recomb_file, string mut_file) {
    read_nodes(node_file);
    read_branches(branch_file);
    read_recombs(recomb_file);
    read_muts(mut_file);
}

void ARG::impute(map<double, Branch> &new_joining_branches, map<double, Branch> &added_branches) {
    double start = added_branches.begin()->first;
    double end = added_branches.rbegin()->first;
    auto mut_it = mutation_sites.lower_bound(start);
    auto join_it = new_joining_branches.begin();
    auto add_it = added_branches.begin();
    double m = 0;
    Branch joining_branch = Branch();
    Branch added_branch = Branch();
    Branch table_joining = Branch();
    Branch table_added = Branch();
    double table_theta = -1;
    double p[32];
    int bi = get_index(start);
    while (add_it->first < end) {
        added_branch = add_it->second;
        if (join_it->first == add_it->first) {
            joining_branch = join_it->second;
            join_it++;
        }
        add_it++;
        while (*mut_it < add_it->first) {
            m = *mut_it;
            while (coordinates[bi + 1] <= m) {
                bi++;
            }
            double unit_theta = thetas[bi]/(coordinates[bi + 1] - coordinates[bi]);
            if (unit_theta != table_theta or joining_branch != table_joining or added_branch != table_added) {
                joining_state_table(joining_branch, added_branch, unit_theta, p);
                table_theta = unit_theta;
                table_joining = joining_branch;
                table_added = added_branch;
            }
            map_mutation(m, joining_branch, added_branch, p);
            mut_it++;
        }
    }
}

void ARG::joining_state_table(const Branch &joining_branch, const Branch &added_branch, double unit_theta, double *p) {
    double tm = added_branch.upper_node->time;
    double ll = tm - joining_branch.lower_node->time;
    double lu = joining_branch.upper_node->time - tm;
    double l0 = tm - added_branch.lower_node->time;
    double u0 = ll*unit_theta;
    double u1 = isinf(lu) ? 1.0 : lu*unit_theta;
    double u2 = l0*unit_theta;
    double w0 = isinf(lu) ? ancestral_prob : 1.0;
    double w1 = isinf(lu) ? 1 - ancestral_prob : 1.0;
    double pen[4] = {1.0, penalty, penalty*penalty, penalty*penalty*penalty};
    for (int c = 0; c < 8; c++) {
        int sl = c & 1, su = (c >> 1) & 1, s0 = (c >> 2) & 1;
        int base = abs(sl - su);
        int k0 = sl + su + s0;
        double t0 = w0*(sl ? u0 : 1.0)*(su ? u1 : 1.0)*(s0 ? u2 : 1.0)*pen[k0 - base];
        double t1 = w1*(sl ? 1.0 : u0)*(su ? 1.0 : u1)*(s0 ? 1.0 : u2)*pen[3 - k0 - base];
        p[c] = t1/(t0 + t1);
    }
    if (!any_missing) {
        return;
    }
    double p0 = u0*penalty, p1 = u1*penalty, p2 = u2*penalty;
    for (int mask = 1; mask < 4; mask++) {
        bool ml = mask & 1, m0 = mask & 2;
        for (int c = 0; c < 8; c++) {
            int sl = c & 1, su = (c >> 1) & 1, s0 = (c >> 2) & 1;
            double t0 = w0*(ml or sl == 0 ? 1.0 : p0)*(su == 0 ? 1.0 : p1)*(m0 or s0 == 0 ? 1.0 : p2);
            double t1 = w1*(ml or sl == 1 ? 1.0 : p0)*(su == 1 ? 1.0 : p1)*(m0 or s0 == 1 ? 1.0 : p2);
            p[8*mask + c] = t1/(t0 + t1);
        }
    }
}

void ARG::compute_assayed() {
    int n = (int) coordinates.size() - 1;
    assayed.resize(n);
    for (int i = 0; i < n; i++) {
        assayed[i] = coordinates[i+1] - coordinates[i];
    }
    for (double x : unassayed_sites) {
        if (x >= 0 and x < sequence_length) {
            assayed[get_index(x)] -= 1;
        }
    }
    for (auto &m : masked) {
        double lo = max(0.0, m.first);
        double hi = min(sequence_length, m.second);
        for (int i = get_index(lo); lo < hi and i < n and coordinates[i] < hi; i++) {
            assayed[i] -= min(hi, coordinates[i+1]) - max(lo, coordinates[i]);
        }
    }
    for (double &a : assayed) {
        a = max(0.0, a);
    }
}

void ARG::map_mutation(double x, Branch joining_branch, Branch added_branch, const double *joining_state_prob) {
    double sl, su, s0, sm;
    Branch new_branch;
    Node *lower_node = joining_branch.lower_node;
    Node *query_node = added_branch.lower_node;
    bool ml = any_missing and lower_node->is_missing(x);
    bool m0 = any_missing and query_node->is_missing(x);
    sl = lower_node->get_state(x);
    su = joining_branch.upper_node->get_state(x);
    s0 = query_node->get_state(x);
    int c = ((int) sl) | (((int) su) << 1) | (((int) s0) << 2) | (ml ? 8 : 0) | (m0 ? 16 : 0);
    sm = uniform_random() < joining_state_prob[c] ? 1 : 0;
    added_branch.upper_node->write_state(x, sm);
    if (sl != su) {
        mutation_branches[x].erase(joining_branch);
    }
    if (sm != sl and !ml) {
        new_branch = Branch(lower_node, added_branch.upper_node);
        mutation_branches[x].insert(new_branch);
    }
    if (sm != su) {
        new_branch = Branch(added_branch.upper_node, joining_branch.upper_node);
        mutation_branches[x].insert(new_branch);
    }
    if (sm != s0 and !m0) {
        mutation_branches[x].insert(added_branch);
    }
    for (const Branch &b : mutation_branches[x]) {
        assert(b.lower_node->get_state(x) != b.upper_node->get_state(x));
    }
}

void ARG::remap_mutations() {
    double x = joining_branches.begin()->first;
    double y = joining_branches.rbegin()->first;
    auto mut_it = mutation_branches.lower_bound(x);
    auto join_it = joining_branches.begin();
    auto remove_it = removed_branches.begin();
    Node *joining_node = nullptr;
    Branch joining_branch = Branch();
    Branch removed_branch = Branch();
    Branch lower_branch = Branch();
    Branch upper_branch = Branch();
    while (mut_it->first < y) {
        while (join_it->first <= mut_it->first) {
            joining_branch = join_it->second;
            assert(joining_branch != Branch());
            join_it++;
        }
        while (remove_it->first <= mut_it->first) {
            removed_branch = remove_it->second;
            joining_node = removed_branch.upper_node;
            assert(joining_node != nullptr);
            remove_it++;
        }
        assert(removed_branch != Branch());
        lower_branch = Branch(joining_branch.lower_node, joining_node);
        upper_branch = Branch(joining_node, joining_branch.upper_node);
        mut_it->second.erase(removed_branch);
        mut_it->second.erase(lower_branch);
        mut_it->second.erase(upper_branch);
        Node *lower_node = lower_branch.lower_node;
        double sl = lower_node->get_state(mut_it->first);
        double su = upper_branch.upper_node->get_state(mut_it->first);
        bool ml = any_missing and lower_node->is_missing(mut_it->first);
        if (sl != su and !ml) {
            mut_it->second.insert(joining_branch);
        }
        for (const Branch &b : mut_it->second) {
            double m = mut_it->first;
            assert(b.lower_node->get_state(m) != b.upper_node->get_state(m));
        }
        mut_it++;
    }
}

int ARG::num_unmapped() {
    int count = 0;
    for (auto &x : mutation_branches) {
        set<Branch> &branches = x.second;
        if (branches.size() > 1) {
            if (branches.rbegin()->upper_node == root.get()) {
                if (branches.size() > 2) {
                    count += 1;
                }
            } else {
                if (branches.size() > 1) {
                    count += 1;
                }
            }
        }
    }
    return count;
}

void ARG::check_incompatibility() {
    int count = 0;
    for (auto &x : mutation_branches) {
        set<Branch> &branches = x.second;
        if (branches.size() > 1) {
            if (branches.rbegin()->upper_node == root.get()) {
                if (branches.size() > 2) {
                    count += 1;
                }
            } else {
                if (branches.size() > 1) {
                    count += 1;
                }
            }
        }
    }
    cout << "Number of incompatibilities: " << count << endl;
}

void ARG::clear_remove_info() {
    removed_branches.clear();
    joining_branches.clear();
    cut_node = nullptr;
}

double ARG::site_weight(const Flat_tree &tree, int bin, double pos, Node *summed) {
    double theta = thetas[bin]/(coordinates[bin + 1] - coordinates[bin]);
    Node *root_node = root.get();
    double w[2] = {1.0, 1.0};
    for (int k = 0; k < (int) tree.parents.size(); k++) {
        Node *l = tree.parents[k].first;
        Node *u = tree.parents[k].second;
        if (u == root_node) {
            for (int s = 0; s < 2; s++) {
                double sm = (l == summed) ? s : l->get_state(pos);
                w[s] *= (sm == 0) ? ancestral_prob : (1 - ancestral_prob)*penalty;
            }
            continue;
        }
        double length = penalty*theta*tree.lengths[k];
        if (any_missing and l->missing_sites.size() > 0 and l->is_missing(pos)) {
            w[0] *= 1 + length;
            w[1] *= 1 + length;
            continue;
        }
        if (l != summed and u != summed) {
            if (l->get_state(pos) != u->get_state(pos)) {
                w[0] *= length;
                w[1] *= length;
            }
            continue;
        }
        for (int s = 0; s < 2; s++) {
            double sl = (l == summed) ? s : l->get_state(pos);
            double su = (u == summed) ? s : u->get_state(pos);
            if (sl != su) {
                w[s] *= length;
            }
        }
    }
    return (summed == nullptr) ? w[0] : w[0] + w[1];
}

double ARG::mutation_log_likelihood(map<double, Branch> &lineage, double x, double y, double pos_p, const Tree &tree_p) {
    Flat_tree tree;
    tree.assign(tree_p);
    for (auto it = recombinations.lower_bound(pos_p); it != recombinations.end() and it->first <= x; ++it) {
        tree.forward_update(it->second);
    }
    auto recomb_it = recombinations.upper_bound(x);
    auto mut_it = mutation_sites.lower_bound(x);
    auto lin_it = lineage.begin();
    double ll = 0;
    double crho = -1;
    double p = 1.0;
    int since_full = 0;
    int i_end = get_index(y);
    for (int i = get_index(x); i < i_end; i++) {
        double w = coordinates[i + 1] - coordinates[i];
        double c = penalty*thetas[i]/w;
        if (coordinates[i] == recomb_it->first) {
            Recombination &r = recomb_it->second;
            tree.forward_update(r);
            recomb_it++;
            if (c == crho and since_full < 256) {
                for (const Branch &b : r.deleted_branches) {
                    if (b.upper_node != root.get()) p /= 1 + crho*(b.upper_node->time - b.lower_node->time);
                }
                for (const Branch &b : r.inserted_branches) {
                    if (b.upper_node != root.get()) p *= 1 + crho*(b.upper_node->time - b.lower_node->time);
                }
                since_full += 1;
            } else {
                crho = -1;
            }
        }
        if (c != crho) {
            crho = c;
            p = 1.0;
            Node *root_node = root.get();
            for (int k = 0; k < (int) tree.parents.size(); k++) {
                if (tree.parents[k].second != root_node) p *= 1 + crho*tree.lengths[k];
            }
            since_full = 0;
        }
        ll -= assayed[i]*log1p((p - 1)/penalty);
        double q = coordinates[i] + 0.5*w;
        while (lin_it != lineage.end() and lin_it->first < q) {
            ++lin_it;
        }
        Node *summed = (lin_it == lineage.begin()) ? Branch().upper_node : prev(lin_it)->second.upper_node;
        while (*mut_it < coordinates[i + 1]) {
            ll += log(site_weight(tree, i, *mut_it, summed));
            mut_it++;
        }
    }
    return ll;
}

double ARG::corrected_smc_prior(map<double, Branch> &lineage, int lo, int hi, double p, const Tree &tree_p) {
    Flat_tree tree;
    tree.assign(tree_p);
    for (auto it = recombinations.lower_bound(p); it != recombinations.end() and it->first <= coordinates[lo]; ++it) {
        tree.forward_update(it->second);
    }
    double log_likelihood = 0.0;
    if (lo == 0) {
        Tree full_tree = tree_p;
        for (auto it = recombinations.lower_bound(p); it != recombinations.end() and it->first <= coordinates[lo]; ++it) {
            full_tree.forward_update(it->second);
        }
        log_likelihood = full_tree.prior_likelihood();
    }
    RSP_smc rsp = RSP_smc();
    rsp.set_tree(tree);
    double visible_length = tree.length() - rsp.unchanged_recomb_length(tree);
    auto recomb_it = recombinations.upper_bound(coordinates[lo]);
    for (int i = lo; i < hi; i++) {
        double rho = rhos[i];
        if (coordinates[i+1] == recomb_it->first) {
            Recombination &rec = recomb_it->second;
            recomb_it++;
            log_likelihood += log(rho);
            log_likelihood += (rec.start_time > 0) ? rsp.log_start_density(rec) : rsp.log_start_marginal(rec, cut_time, lineage_branch_before(lineage, rec.pos));
            tree.forward_update(rec);
            rsp.set_tree(tree);
            visible_length = tree.length() - rsp.unchanged_recomb_length(tree);
        } else {
            if (rho*visible_length >= 1) {
                cerr << "bin " << i << " is too wide for one record per bin: rho*visible length = " << rho*visible_length << endl;
                exit(1);
            }
            log_likelihood += log1p(-rho*visible_length);
        }
    }
    return log_likelihood;
}

set<double> ARG::get_check_points() {
    double start_pos = removed_branches.begin()->first;
    double end_pos = removed_branches.rbegin()->first;
    auto recomb_it = recombinations.lower_bound(start_pos);
    map<Node *, double> deleted_nodes = {};
    vector<tuple<Node *, double, double>> node_span = {};
    while (recomb_it->first <= end_pos) {
        Recombination &r = recomb_it->second;
        deleted_nodes[r.deleted_node] = r.pos;
        Node *inserted_node = r.inserted_node;
        if (deleted_nodes.count(inserted_node) > 0 and inserted_node != root.get()) {
            node_span.push_back({inserted_node, deleted_nodes[inserted_node], r.pos});
            deleted_nodes.erase(inserted_node);
        }
        recomb_it++;
    }
    Node *n;
    double x;
    double y;
    set<double> check_points = {};
    for (auto ns : node_span) {
        tie(n, x, y) = ns;
        if (!check_disjoint_nodes(x, y)) {
            check_points.insert(y);
        }
    }
    return check_points;
}

bool ARG::check_disjoint_nodes(double x, double y) {
    auto recomb_it = recombinations.lower_bound(x);
    double t = recomb_it->second.deleted_node->time;
    Branch b = recomb_it->second.merging_branch;
    while (recomb_it->first < y) {
        b = recomb_it->second.trace_forward(t, b);
        if (b == Branch()) {
            return false;
        }
        recomb_it++;
    }
    if (b != recomb_it->second.target_branch) {
        return false;
    }
    return true;
}

void ARG::new_recombination(double pos, Branch prev_added_branch, Branch prev_joining_branch, Branch next_added_branch, Branch next_joining_branch) {
    set<Branch> deleted_branches;
    set<Branch> inserted_branches;
    deleted_branches.insert(prev_added_branch);
    deleted_branches.insert(Branch(prev_joining_branch.lower_node, prev_added_branch.upper_node));
    deleted_branches.insert(Branch(prev_added_branch.upper_node, prev_joining_branch.upper_node));
    deleted_branches.insert(next_joining_branch);
    inserted_branches.insert(next_added_branch);
    inserted_branches.insert(Branch(next_joining_branch.lower_node, next_added_branch.upper_node));
    inserted_branches.insert(Branch(next_added_branch.upper_node, next_joining_branch.upper_node));
    inserted_branches.insert(prev_joining_branch);
    Recombination r = Recombination(deleted_branches, inserted_branches);
    r.set_pos(pos);
    recombinations[pos] = r;
    return;
}

void ARG::remove_empty_recombinations() {
    auto recomb_it = recombinations.lower_bound(start);
    while (recomb_it->first <= end) {
        Recombination &r = recomb_it->second;
        if (r.deleted_branches.size() == 0 and r.inserted_branches.size() == 0 and r.pos < sequence_length) {
            recomb_it = recombinations.erase(recomb_it);
        } else {
            ++recomb_it;
        }
    }
}

void ARG::release_dead_nodes() {
    unordered_set<Node *> live = {cut_node};
    auto keep_branch = [&](const Branch &b) {
        live.insert(b.lower_node);
        live.insert(b.upper_node);
    };
    auto keep_tree = [&](const Tree &t) {
        for (auto &x : t.parents) {
            live.insert(x.first);
            live.insert(x.second);
        }
    };
    for (auto &x : recombinations) {
        const Recombination &r = x.second;
        for (const Branch *b : {&r.source_branch, &r.target_branch, &r.source_sister_branch, &r.source_parent_branch,
                                &r.recombined_branch, &r.merging_branch, &r.lower_transfer_branch, &r.upper_transfer_branch}) {
            keep_branch(*b);
        }
        for (const Branch &b : r.deleted_branches) {
            keep_branch(b);
        }
        for (const Branch &b : r.inserted_branches) {
            keep_branch(b);
        }
        live.insert(r.deleted_node);
        live.insert(r.inserted_node);
    }
    for (Node *n : sample_nodes) {
        live.insert(n);
    }
    for (auto &x : mutation_branches) {
        for (const Branch &b : x.second) {
            keep_branch(b);
        }
    }
    for (auto *m : {&joining_branches, &removed_branches}) {
        for (auto &x : *m) {
            keep_branch(x.second);
        }
    }
    for (auto &x : anchors) {
        keep_tree(x.second);
    }
    keep_tree(cut_tree);
    keep_tree(start_tree);
    keep_tree(end_tree);
    auto first_dead = partition(node_owner.begin(), node_owner.end(), [&](const Node_ptr &p) { return live.count(p.get()) > 0; });
    dead_nodes.assign(make_move_iterator(first_dead), make_move_iterator(node_owner.end()));
    node_owner.erase(first_dead, node_owner.end());
}

void ARG::create_node_set() {
    node_set.clear();
    for (auto &x : recombinations) {
        for (const Branch &b : x.second.inserted_branches) {
            add_node(b.upper_node);
        }
    }
    for (Node *n : sample_nodes) {
        add_node(n);
    }
    assert(cut_node == nullptr or node_set.count(cut_node) > 0);
}

void ARG::write_nodes(string filename) {
    anchors_changed(0, INT_MAX);
    node_set.clear();
    create_node_set();
    vector<Node *> order;
    for (Node *n : node_set) {
        if (n->is_sample) {
            order.push_back(n);
        }
    }
    sort(order.begin(), order.end(), [](Node *a, Node *b) { return a->index < b->index; });
    for (Node *n : node_set) {
        if (!n->is_sample) {
            order.push_back(n);
        }
    }
    ofstream file;
    file.open(filename);
    int index = 0;
    for (Node *n : order) {
        if (!n->is_sample) {
            n->set_index(index);
        }
        file << std::setprecision(std::numeric_limits<double>::max_digits10) << (n->time + time_offset)*Ne << "\n";
        index += 1;
    }
    file.close();
    node_set.clear();
}

void ARG::write_branches(string filename) {
    map<Branch, double> branch_map;
    vector<tuple<double, double, double, double>> branch_info;
    double pos;
    for (auto x : recombinations) {
        if (x.first < sequence_length) {
            pos = x.first;
            Recombination &r = x.second;
            for (Branch b : r.inserted_branches) {
                branch_map[b] = pos;
            }
            for (Branch b : r.deleted_branches) {
                int k1 = b.upper_node->index;
                int k2 = b.lower_node->index;
                branch_info.push_back({k1, k2, branch_map.at(b), pos});
                branch_map.erase(b);
            }
        }
    }
    for (auto x : branch_map) {
        Branch b = x.first;
        int k1 = b.upper_node->index;
        int k2 = b.lower_node->index;
        branch_info.push_back({k1, k2, x.second, sequence_length});
    }
    sort(branch_info.begin(), branch_info.end(), compare_edge);
    ofstream file;
    file.open(filename);
    file << std::setprecision(std::numeric_limits<double>::max_digits10) << std::fixed;
    for (int i = 0; i < branch_info.size(); i++) {
        auto [k1, k2, x, l] = branch_info[i];
        file << x << " " << l << " " << k1 << " " << k2 << "\n";
    }
    file.close();
}

void ARG::write_recombs(string filename) {
    ofstream file;
    file.open(filename);
    file << std::setprecision(std::numeric_limits<double>::max_digits10) << std::fixed;
    for (auto x : recombinations) {
        Recombination &r = x.second;
        if (x.first > 0 and x.first < sequence_length) {
            file << r.pos << " " << r.source_branch.lower_node->index << " " << r.source_branch.upper_node->index << " " << Ne*(r.start_time + time_offset) << "\n";
        }
    }
    file.close();
}

void ARG::write_mutations(string filename) {
    ofstream file;
    file.open(filename);
    file.precision(numeric_limits<double>::max_digits10);
    for (auto &x : mutation_branches) {
        double m = x.first;
        for (auto &y : x.second) {
            if (m < sequence_length and m > 0) {
                file << m << " " << y.lower_node->index << " " << y.upper_node->index << " " << y.lower_node->get_state(m) << "\n";
            }
        }
    }
}

void ARG::read_nodes(string filename) {
    root->set_index(-1);
    node_set.clear();
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    vector<double> times;
    double x;
    while (fin >> x) {
        times.push_back(x);
    }
    double youngest = *min_element(times.begin(), times.begin() + num_samples);
    time_offset = youngest/Ne;
    for (int i = 0; i < (int) times.size(); i++) {
        add_new_node((times[i] - youngest)/Ne, i < num_samples);
    }
}

void ARG::read_branches(string filename) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    vector<Node *> nodes = node_list;
    double x;
    double y;
    double p;
    double c;
    double left;
    double right;
    Node *parent_node;
    Node *child_node;
    Branch b;
    map<double, set<Branch>> deleted_branches = {{0, {}}};
    map<double, set<Branch>> inserted_branches = {};
    while (fin >> x >> y >> p >> c) {
        left = x;
        right = y;
        if (p < 0) {
            parent_node = root.get();
        } else {
            parent_node = nodes[int(p)];
        }
        child_node = nodes[int(c)];
        b = Branch(child_node, parent_node);
        if (deleted_branches.count(right) > 0) {
            deleted_branches[right].insert(b);
        } else {
            deleted_branches[right] = {b};
        }
        if (inserted_branches.count(left) > 0) {
            inserted_branches[left].insert(b);
        } else {
            inserted_branches[left] = {b};
        }
    }
    deleted_branches.erase(sequence_length);
    for (auto x : deleted_branches) {
        double pos = x.first;
        set<Branch> db = deleted_branches.at(pos);
        set<Branch> ib = inserted_branches.at(pos);
        Recombination r = Recombination(db, ib);
        r.set_pos(pos);
        recombinations[pos] = r;
    }
}

void ARG::read_recombs(string filename) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    create_node_set();
    vector<Node *> nodes = node_list;
    map<double, Branch> source_branches = {};
    map<double, double> start_times = {};
    double pos;
    int n1;
    int n2;
    double t;
    Node *ln;
    Node *un;
    Branch b;
    while (fin >> pos >> n1 >> n2 >> t) {
        ln = nodes[n1];
        if (n2 == -1) {
            un = root.get();
        } else {
            un = nodes[n2];
        }
        b = Branch(ln, un);
        source_branches[pos] = b;
        start_times[pos] = t;
    }
    for (auto &x : recombinations) {
        pos = x.first;
        if (pos > 0 and pos < sequence_length) {
            t = start_times.at(pos);
            b = source_branches.at(pos);
            if (pos > 0 and pos < sequence_length) {
                x.second.start_time = t/Ne - time_offset;
                x.second.source_branch = b;
                x.second.find_nodes();
                x.second.find_target_branch();
                x.second.find_recomb_info();
            }
        }
        if (pos > 0 and pos < sequence_length) {
            assert(x.second.start_time <= x.second.deleted_node->time);
            assert(x.second.start_time <= x.second.inserted_node->time);
        }
    }
    node_set.clear();
}

void ARG::read_muts(string filename) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    create_node_set();
    vector<Node *> nodes = node_list;
    double pos;
    int n1;
    int n2;
    double s;
    Node *ln;
    Node *un;
    Branch b;
    while (fin >> pos >> n1 >> n2 >> s) {
        if (pos <= sequence_length) {
            mutation_sites.insert(pos);
            ln = nodes[n1];
            if (n2 == -1) {
                un = root.get();
            } else {
                un = nodes[n2];
            }
            if (s == 1) {
                ln->write_state(pos, 1);
                un->write_state(pos, 0);
            } else {
                un->write_state(pos, 1);
                ln->write_state(pos, 0);
            }
            b = Branch(ln, un);
            mutation_branches[pos].insert(b);
        }
    }
    Tree tree = Tree();
    auto m_it = mutation_branches.begin();
    auto r_it = recombinations.begin();
    while (r_it->first < sequence_length) {
        Recombination &r = r_it->second;
        tree.forward_update(r);
        while (m_it->first < next(r_it)->first) {
            tree.impute_states(m_it->first, m_it->second);
            m_it++;
        }
        r_it++;
    }
}

double ARG::get_arg_length() {
    Tree tree = get_tree_at(0);
    auto recomb_it = recombinations.upper_bound(0);
    double l = 0, span = 0;
    double tree_length = tree.length();
    double prev_pos = 0;
    double next_pos = 0;
    while (next_pos <= sequence_length) {
        next_pos = recomb_it->first;
        span = min(sequence_length, next_pos) - prev_pos;
        l += tree_length*span;
        Recombination &r = recomb_it->second;
        recomb_it++;
        tree.forward_update(r);
        tree_length = tree.length();
        prev_pos = next_pos;
    }
    return l;
}

tuple<double, Branch, double> ARG::sample_internal_cut() {
    if (end >= sequence_length - 0.1) {
        cut_pos = 0;
        cut_tree = get_tree_at(0);
    } else {
        cut_tree = move(end_tree);
        cut_pos = end;
    }
    Branch b;
    double t;
    tie(b, t) = cut_tree.sample_cut_point();
    while (t == b.lower_node->time or t == b.upper_node->time) {
        tie(b, t) = cut_tree.sample_cut_point();
    }
    return {cut_pos, b, t};
}

tuple<double, Branch, double> ARG::sample_uniform_cut() {
    cut_pos = uniform_random()*sequence_length;
    cut_tree = get_tree_at(cut_pos);
    Branch b;
    double t;
    tie(b, t) = cut_tree.sample_uniform_cut_point();
    while (t == b.lower_node->time or t == b.upper_node->time) {
        tie(b, t) = cut_tree.sample_uniform_cut_point();
    }
    return {cut_pos, b, t};
}

bool compare_edge(const tuple<int, int, double, double>& edge1, const tuple<int, int, double, double>& edge2) {
    if (get<0>(edge1) < get<0>(edge2)) {
        return true;
    } else if (get<0>(edge1) > get<0>(edge2)) {
        return false;
    }
    if (get<1>(edge1) < get<1>(edge2)) {
        return true;
    } else if (get<1>(edge1) > get<1>(edge2)) {
        return false;
    }
    return get<2>(edge1) < get<2>(edge2);
}

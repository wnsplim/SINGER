//
//  Threader_smc.cpp
//  SINGER
//
//  Created by Yun Deng on 4/4/23.
//

#include "Threader_smc.hpp"

bool Threader_smc::no_data = false;

Threader_smc::Threader_smc(double c, double q) {
    cutoff = c;
    gap = q;
}

Threader_smc::~Threader_smc() {
}

void Threader_smc::reset() {
    bsp.reset();
    tsp.reset();
    new_joining_branches.clear();
    added_branches.clear();
    cut_time = 0;
    start = 0;
    end = 0;
    start_index = 0;
    end_index = 0;
}

static bool held_between(ARG &a, Node *n, double x, double y) {
    auto it = a.recombinations.find(x);
    if (it == a.recombinations.end() or it->second.inserted_node != n) {
        return false;
    }
    for (++it; it != a.recombinations.end() and it->first < y; ++it) {
        if (it->second.deleted_node == n) {
            return false;
        }
    }
    return true;
}

static bool has_gap_node(const map<double, Branch> &branches, ARG &a) {
    map<Node *, double> left_at;
    Node *current = nullptr;
    for (auto &x : branches) {
        Node *n = x.second.upper_node;
        if (n == current) {
            continue;
        }
        if (current != nullptr) {
            left_at[current] = x.first;
        }
        auto it = left_at.find(n);
        if (it != left_at.end() and !held_between(a, n, it->second, x.first)) {
            return true;
        }
        current = n;
    }
    return false;
}

void Threader_smc::thread(ARG &a, Node_ptr n) {
    cout << "Iteration: " << a.sample_nodes.size() << endl;
    cut_time = n->time;
    a.cut_time = cut_time;
    a.add_sample(n);
    get_boundary(a);
    cout << get_time() << " : begin BSP" << endl;
    run_BSP(a);
    cout << "BSP avg num of states: " << bsp.avg_num_states() << endl;
    cout << get_time() << " : begin sampling branches" << endl;
    sample_joining_branches(a);
    cout << get_time() << " : begin TSP" << endl;
    run_TSP(a);
    cout << get_time() << " : begin sampling points" << endl;
    sample_joining_points(a);
    cout << get_time() << " : begin adding" << endl;
    a.add(new_joining_branches, added_branches);
    cout << get_time() << " : begin sampling recombination" << endl;
    a.approx_sample_recombinations();
    a.clear_remove_info();
    cout << get_time() << " : finish" << endl;
    cout << a.recombinations.size() << endl;
}

void Threader_smc::internal_rethread(ARG &a, tuple<double, Branch, double> cut_point) {
    cut_time = get<2>(cut_point);
    a.remove(cut_point);
    get_boundary(a);
    set_check_points(a);
    run_BSP(a);
    sample_joining_branches(a);
    run_TSP(a);
    sample_joining_points(a);
    double ar = has_gap_node(added_branches, a) ? 0 : acceptance_ratio(a);
    double q = random();
    if (q < ar) {
        a.add(new_joining_branches, added_branches);
    } else {
        a.add(a.joining_branches, a.removed_branches);
    }
    a.approx_sample_recombinations();
    a.clear_remove_info();
}

void Threader_smc::exact_internal_rethread(ARG &a, tuple<double, Branch, double> cut_point) {
    cut_time = get<2>(cut_point);
#ifndef NDEBUG
    full_length = a.get_tree_at(get<0>(cut_point)).length();
#endif
    a.remove(cut_point);
    get_boundary(a);
    set_check_points(a);
    double ar = 0;
    if (!has_bridges(a)) {
        run_BSP(a);
        sample_joining_branches(a);
        run_TSP(a);
        sample_joining_points(a);
        ar = bridges_kept(a) and !has_gap_node(added_branches, a) ? exact_acceptance_ratio(a) : 0;
    }
    double q = random();
    if (q < ar) {
        a.add(new_joining_branches, added_branches);
        a.smc_sample_recombinations(added_branches);
    } else {
        a.add(a.joining_branches, a.removed_branches);
        a.smc_sample_recombinations(a.removed_branches);
    }
    a.clear_remove_info();
}

static bool joined_at(map<double, Branch> &m, double x, double y, Node *n) {
    auto it = m.upper_bound(x);
    --it;
    for (; it != m.end() and it->first < y; ++it) {
        if (it->second.upper_node != n) {
            return false;
        }
    }
    return true;
}

bool Threader_smc::has_bridges(ARG &a) {
    map<Node *, double> deleted_at;
    for (auto it = a.recombinations.upper_bound(a.start); it != a.recombinations.end() and it->first <= a.end; ++it) {
        Recombination &r = it->second;
        auto d = deleted_at.find(r.inserted_node);
        if (d != deleted_at.end()) {
            if (joined_at(a.removed_branches, d->second, it->first, r.inserted_node)) {
                return true;
            }
            deleted_at.erase(d);
        }
        deleted_at[r.deleted_node] = it->first;
    }
    return false;
}

bool Threader_smc::bridges_kept(ARG &a) {
    map<Node *, double> deleted_at;
    for (auto it = a.recombinations.upper_bound(a.start); it != a.recombinations.end() and it->first <= a.end; ++it) {
        Recombination &r = it->second;
        auto d = deleted_at.find(r.inserted_node);
        if (d != deleted_at.end()) {
            if (joined_at(a.removed_branches, d->second, it->first, r.inserted_node) and !joined_at(added_branches, d->second, it->first, r.inserted_node)) {
                return false;
            }
            deleted_at.erase(d);
        }
        deleted_at[r.deleted_node] = it->first;
    }
    return true;
}

void Threader_smc::get_boundary(ARG &a) {
    start = a.start;
    end = a.end;
    start_index = a.get_index(start);
    end_index = a.get_index(end);
}

void Threader_smc::set_check_points(ARG &a) {
    set<double> check_points = a.get_check_points();
    bsp.set_check_points(check_points);
    tsp.set_check_points(check_points);
}

void Threader_smc::run_BSP(ARG &a) {
    bsp.reserve_memory(end_index - start_index);
    bsp.set_cutoff(cutoff);
    bsp.set_emission(pe);
    bsp.start(a.start_tree, cut_time);
    auto recomb_it = a.recombinations.upper_bound(start);
    auto mut_it = a.mutation_sites.lower_bound(start);
    auto query_it = a.removed_branches.begin();
    vector<double> mutations;
    vector<double> mut_set = {};
    set<Branch> deletions = {};
    set<Branch> insertions = {};
    Node *query_node = nullptr;
    bool varying = varying_rate(a, start_index, end_index);
    Tree tree;
    if (varying) {
        tree = a.start_tree;
    }
    double curr_rho = a.thetas[start_index]/(a.coordinates[start_index + 1] - a.coordinates[start_index]);
    double p_tree = branch_product(a.start_tree, pe->penalty, curr_rho, a.cut_node);
    pe->any_missing = a.any_missing;
    for (int i = start_index; i < end_index; i++) {
        if (a.coordinates[i] == query_it->first) {
            query_node = query_it->second.lower_node;
            query_it++;
        }
        if (a.coordinates[i] == recomb_it->first) {
            Recombination &r = recomb_it->second;
            recomb_it++;
            bsp.transfer(r);
            if (varying) {
                tree.forward_update(r);
                curr_rho = -1;
            } else {
                p_tree = update_branch_product(p_tree, r, pe->penalty*curr_rho, a.cut_node);
            }
        } else if (a.coordinates[i] != start) {
            bsp.forward(a.rhos[i - 1]);
        }
        double w = a.coordinates[i + 1] - a.coordinates[i];
        double rho = a.thetas[i]/w;
        if (varying and rho != curr_rho) {
            p_tree = branch_product(tree, pe->penalty, rho, a.cut_node);
            curr_rho = rho;
        }
        pe->set_tree_product(pe->penalty, no_data ? 0.0 : rho, p_tree, called_width(query_node, a, i));
        mut_set.clear();
        while (*mut_it < a.coordinates[i + 1]) {
            mut_set.push_back(*mut_it);
            mut_it++;
        }
        if (mut_set.size() > 0 and !no_data) {
            bsp.mut_emit(a.thetas[i], w, mut_set, query_node);
        } else {
            bsp.null_emit(a.thetas[i], query_node);
        }
    }
    if (bsp.check_points.count(end) > 0) {
        Recombination &r = a.recombinations[end];
        bsp.sanity_check(r);
    }
}

void Threader_smc::run_TSP(ARG &a) {
    run_TSP(a, new_joining_branches);
}

void Threader_smc::run_TSP(ARG &a, map<double, Branch> &jb) {
    tsp.reserve_memory(end_index - start_index);
    tsp.set_gap(gap);
    tsp.set_emission(be);
    tsp.cc = make_shared<coalescent_calculator>(cut_time);
    tsp.cc->start(a.start_tree);
    tsp.cc1 = make_shared<coalescent_calculator>(cut_time);
    tsp.cc1->extra = 1;
    tsp.cc1->start(a.start_tree);
    Branch start_branch = jb.begin()->second;
    tsp.start(start_branch, cut_time);
    auto recomb_it = a.recombinations.upper_bound(start);
    auto join_it = jb.upper_bound(start);
    auto mut_it = a.mutation_sites.lower_bound(start);
    auto query_it = a.removed_branches.lower_bound(start);
    Branch prev_branch = start_branch;
    Branch next_branch = start_branch;
    Node *query_node = nullptr;
    vector<double> mut_set = {};
    bool varying = varying_rate(a, start_index, end_index);
    Tree tree = a.start_tree;
    RSP_smc rsp;
    auto set_sister_mass = [&](Branch &b) {
        double lo = b.lower_node->time;
        if (lo < cut_time) {
            rsp.set_tree(tree);
            tsp.sister_mass = rsp.sister_mass(lo, cut_time);
        } else {
            tsp.sister_mass = 0;
        }
    };
    set_sister_mass(start_branch);
    double curr_rho = a.thetas[start_index]/(a.coordinates[start_index + 1] - a.coordinates[start_index]);
    double p_tree = branch_product(a.start_tree, be->penalty, curr_rho, a.cut_node);
    be->any_missing = a.any_missing;
    for (int i = start_index; i < end_index; i++) {
        if (a.coordinates[i] == query_it->first) {
            query_node = query_it->second.lower_node;
            query_it++;
        }
        if (a.coordinates[i] == join_it->first) {
            next_branch = join_it->second;
            join_it++;
        }
        if (a.coordinates[i] == recomb_it->first) {
            Recombination &r = recomb_it->second;
            recomb_it++;
            tsp.transfer(r, prev_branch, next_branch);
            tree.forward_update(r);
            if (!varying) {
                p_tree = update_branch_product(p_tree, r, be->penalty*curr_rho, a.cut_node);
            } else {
                curr_rho = -1;
            }
            prev_branch = next_branch;
            set_sister_mass(next_branch);
        } else if (prev_branch != next_branch) {
            tsp.recombine(prev_branch, next_branch);
            prev_branch = next_branch;
            set_sister_mass(next_branch);
        } else if (a.coordinates[i] != start) {
            tsp.forward(a.rhos[i - 1]);
        }
        double w = a.coordinates[i+1] - a.coordinates[i];
        double rho = a.thetas[i]/w;
        if (varying and rho != curr_rho) {
            p_tree = branch_product(tree, be->penalty, rho, a.cut_node);
            curr_rho = rho;
        }
        be->set_tree_product(be->penalty, no_data ? 0.0 : rho, p_tree, called_width(query_node, a, i));
        mut_set.clear();
        while (*mut_it < a.coordinates[i+1]) {
            mut_set.push_back(*mut_it);
            mut_it++;
        }
        if (mut_set.size() > 0 and !no_data) {
            tsp.mut_emit(a.thetas[i], w, mut_set, query_node);
        } else {
            tsp.null_emit(a.thetas[i], query_node);
        }
    }
    if (tsp.check_points.count(end) > 0) {
        Recombination &r = a.recombinations[end];
        tsp.sanity_check(r);
    }
}

void Threader_smc::sample_joining_branches(ARG &a) {
    new_joining_branches = bsp.sample_joining_branches(start_index, a.coordinates);
}

void Threader_smc::sample_joining_points(ARG &a) {
    map<double, Node *> added_nodes = tsp.sample_joining_nodes(start_index, a.coordinates);
    a.node_owner.insert(a.node_owner.end(), tsp.node_owner.begin(), tsp.node_owner.end());
    tsp.node_owner.clear();
    auto add_it = added_nodes.begin();
    auto end_it = added_nodes.end();
    Node *query_node = nullptr;
    Node *added_node = nullptr;
    double x;
    while (add_it != end_it) {
        x = add_it->first;
        added_node = add_it->second;
        query_node = a.get_query_node_at(x);
        added_branches[x] = Branch(query_node, added_node);
        add_it++;
    }
}

static double pruned_height(ARG &a) {
    for (auto &x : a.cut_tree.parents) {
        if (x.second == a.root.get()) {
            return x.first->time;
        }
    }
    return 0;
}

double Threader_smc::acceptance_ratio(ARG &a) {
    double cut_height = pruned_height(a);
    double old_height = cut_height;
    double new_height = cut_height;
    auto old_join_it = a.joining_branches.upper_bound(a.cut_pos);
    old_join_it--;
    auto new_join_it = new_joining_branches.upper_bound(a.cut_pos);
    new_join_it--;
    auto old_add_it = a.removed_branches.upper_bound(a.cut_pos);
    old_add_it--;
    auto new_add_it = added_branches.upper_bound(a.cut_pos);
    new_add_it--;
    if (old_join_it->second.upper_node == a.root.get()) {
        old_height = old_add_it->second.upper_node->time;
    }
    if (new_join_it->second.upper_node == a.root.get()) {
        new_height = new_add_it->second.upper_node->time;
    }
    return old_height/new_height;
}

double Threader_smc::cut_ratio(ARG &a) {
    double cut_height = pruned_height(a);
    double pruned_length = a.cut_tree.length();
    auto old_join_it = a.joining_branches.upper_bound(a.cut_pos);
    old_join_it--;
    auto new_join_it = new_joining_branches.upper_bound(a.cut_pos);
    new_join_it--;
    auto old_add_it = a.removed_branches.upper_bound(a.cut_pos);
    old_add_it--;
    auto new_add_it = added_branches.upper_bound(a.cut_pos);
    new_add_it--;
    double old_join_time = old_add_it->second.upper_node->time;
    double new_join_time = new_add_it->second.upper_node->time;
    double old_length = pruned_length + old_join_time - cut_time;
    double new_length = pruned_length + new_join_time - cut_time;
    if (old_join_it->second.upper_node == a.root.get()) {
        old_length += old_join_time - cut_height;
    }
    if (new_join_it->second.upper_node == a.root.get()) {
        new_length += new_join_time - cut_height;
    }
    assert(fabs(old_length - full_length) <= 1e-9*full_length);
    return old_length/new_length;
}

double Threader_smc::exact_acceptance_ratio(ARG &a) {
    double jac = cut_ratio(a);
    double q_new = bsp.branch_log_q(new_joining_branches, start_index, a.coordinates);
    double q_old = bsp.branch_log_q(a.joining_branches, start_index, a.coordinates);
    double h_new = tsp.sel_log_q + tsp.time_log_q;
    if (a.joining_branches != new_joining_branches) {
        set<double> check_points = tsp.check_points;
        tsp.reset();
        tsp.check_points = check_points;
        run_TSP(a, a.joining_branches);
    }
    double h_old = tsp.eval_joining_nodes(a.joining_branches, a.removed_branches, start_index, a.coordinates);
    int lo = max(start_index - 1, 0);
    int hi = min(end_index, a.bin_num - 1);
    pos_lo = a.coordinates[lo];
    tree_lo = a.start_tree;
    {
        auto it = a.recombinations.upper_bound(start);
        while (it != a.recombinations.begin()) {
            --it;
            if (it->first < pos_lo) break;
            tree_lo.backward_update(it->second);
        }
    }
    a.window_copy_into(new_arg, pos_lo);
    new_arg.add(new_joining_branches, added_branches);
    double pi_new = new_arg.corrected_smc_prior(added_branches, lo, hi, pos_lo, tree_lo) + (no_data ? 0.0 : new_arg.mutation_log_likelihood(added_branches, start, end, pos_lo, tree_lo));
    a.window_copy_into(old_arg, pos_lo);
    old_arg.add(a.joining_branches, a.removed_branches);
    double pi_old = old_arg.corrected_smc_prior(a.removed_branches, lo, hi, pos_lo, tree_lo) + (no_data ? 0.0 : old_arg.mutation_log_likelihood(a.removed_branches, start, end, pos_lo, tree_lo));
    double log_a = (pi_new - pi_old) + (q_old - q_new) + (h_old - h_new) + log(jac);
    return isfinite(q_new) and isfinite(h_new) ? exp(log_a) : 0;
}

double Threader_smc::random() {
    return uniform_random();
}

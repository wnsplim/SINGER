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
    resume_index = -1;
    site_slot.clear();
    emit_trees = 0;
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
    set_check_points(a);
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

void Threader_smc::rethread_sample(ARG &a, Node *n) {
    cut_time = n->time + 1e-10;
    a.cut_pos = 0;
    a.cut_tree = a.get_tree_at(0);
    a.remove({0.0, Branch(n, a.cut_tree.parents[n]), cut_time});
    get_boundary(a);
    set_check_points(a);
    run_BSP(a);
    sample_joining_branches(a);
    run_TSP(a);
    sample_joining_points(a);
    a.add(new_joining_branches, added_branches);
    a.approx_sample_recombinations();
    a.clear_remove_info();
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
    Node *query = a.removed_branches.begin()->second.lower_node;
    a.set_collapsed(query, a.removed_branches.begin()->first, a.removed_branches.rbegin()->first);
    auto [entry_lo, entry_hi] = a.entries_in(query, start, end);
    bool entries = entry_lo != entry_hi;
    a.journal_on = entries;
    double ar = 0;
    if (!joins_deleted_node(a, a.removed_branches)) {
        run_BSP(a);
        sample_joining_branches(a);
        if (!entries) {
            find_resume_index(a);
        }
        run_TSP(a);
        sample_joining_points(a);
        if (keeps_deleted_node(a) and !has_gap_node(added_branches, a) and !joins_deleted_node(a, added_branches) and !a.empties_cut_span_end_recombination(new_joining_branches, added_branches)) {
            ar = entries ? redrawn_acceptance_ratio(a, cut_point) : exact_acceptance_ratio(a);
        }
    }
    a.journal_on = false;
    double q = random();
    if (q < ar and entries) {
        a.add(new_joining_branches, added_branches);
        a.recombinations.erase(a.recombinations.lower_bound(start), a.recombinations.upper_bound(end));
        a.recombinations.insert(redrawn_records.begin(), redrawn_records.end());
    } else if (q < ar) {
        a.add(new_joining_branches, added_branches);
        a.smc_sample_recombinations(added_branches);
    } else {
        a.undo_journal();
        a.add(a.joining_branches, a.removed_branches);
        a.smc_sample_recombinations(a.removed_branches);
    }
    a.state_journal.clear();
    query->likelihood_sites.clear();
    a.clear_remove_info();
}

double Threader_smc::log_pi_ratio(ARG &a) {
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
    double pi_new = new_arg.corrected_smc_prior(added_branches, lo, hi, pos_lo, tree_lo);
    a.window_copy_into(old_arg, pos_lo);
    old_arg.add(a.joining_branches, a.removed_branches);
    double pi_old = old_arg.corrected_smc_prior(a.removed_branches, lo, hi, pos_lo, tree_lo);
    if (no_data) {
        return pi_new - pi_old;
    }
    return pi_new - pi_old + data_log_ratio(a);
}

double Threader_smc::data_log_ratio(ARG &a) {
    Flat_tree &flat = ratio_flat;
    Tree_order &t = ratio_order;
    vector<array<double, 2>> &ref = ratio_ref, &pat = ratio_pattern;
    flat.assign(a.start_tree);
    double prior[2] = {a.ancestral_prob, 1 - a.ancestral_prob};
    double pen = a.penalty;
    const ARG *owner = a.entry_owner == nullptr ? &a : a.entry_owner;
    auto [qi, qe] = a.entries_in(a.removed_branches.begin()->second.lower_node, start, end);
    const vector<Genotype_likelihood> &gl = a.entries();
    auto recomb_it = a.recombinations.upper_bound(start);
    auto mut_it = a.mutation_sites.lower_bound(start);
    auto oj = a.joining_branches.begin(), on = a.removed_branches.begin();
    auto nj = new_joining_branches.begin(), nn = added_branches.begin();
    Branch bo, bn, ro, rn;
    vector<char> main;
    vector<int> path;
    int samples = (int) a.sample_nodes.size();
    vector<char> missing(samples, 0), target(samples, 0), seen(samples, 0);
    vector<int> missing_list, touched;
    int observed_main = 0, observed_sub = 0;
    double unit = -1, p_common = 1;
    bool order_stale = true, ref_valid = false, pattern_valid = false, pattern_value_valid = false, z_valid = false;
    double ref_value = 0, pattern_value = 0, z_diff = 0;

    auto above = [&](const vector<array<double, 2>> &msg, int d) {
        path.clear();
        for (int k = d; k != t.top; k = t.up[k]) {
            path.push_back(k);
        }
        array<double, 2> A = {prior[0], prior[1]};
        int v = t.top;
        for (int j = (int) path.size() - 1; j >= 0; j--) {
            int c = path[j];
            array<double, 2> B = A;
            for (int sib : t.kids[v]) {
                if (sib != c) {
                    double len = unit*flat.lengths[sib];
                    B[0] *= msg[sib][0] + len*msg[sib][1];
                    B[1] *= msg[sib][1] + len*msg[sib][0];
                }
            }
            if (j == 0) {
                return B;
            }
            double len = unit*flat.lengths[c];
            A = {B[0] + len*B[1], B[1] + len*B[0]};
            v = c;
        }
        return A;
    };
    auto join = [&](const Branch &b, double time, Node *query, const array<double, 2> &nu, const array<double, 2> &m_d, const array<double, 2> &up) {
        double l2 = unit*(time - b.lower_node->time), l3 = unit*(time - query->time);
        double m0 = (m_d[0] + l2*m_d[1])*(nu[0] + l3*nu[1]), m1 = (m_d[1] + l2*m_d[0])*(nu[1] + l3*nu[0]);
        if (b.upper_node == a.root.get()) {
            return prior[0]*m0 + prior[1]*m1;
        }
        double l1 = unit*(b.upper_node->time - time);
        return up[0]*(m0 + l1*m1) + up[1]*(m1 + l1*m0);
    };
    auto log_z = [&](const Branch &b, double time, Node *query) {
        double l2 = unit*(time - b.lower_node->time), l3 = unit*(time - query->time);
        double p = p_common*(1 + l2)*(1 + l3);
        if (b.upper_node != a.root.get()) {
            double l1 = unit*(b.upper_node->time - time);
            p *= (1 + l1)/(1 + l1 + l2);
        }
        return log1p((p - 1)/pen);
    };
    auto pattern_log_ratio = [&](const vector<array<double, 2>> &msg) {
        Node *query = ro.lower_node;
        int kq = query->flat_index, kd_o = bo.lower_node->flat_index, kd_n = bn.lower_node->flat_index;
        double e1 = observed_main == 0 and observed_sub == 0 ? 1.0 : 0.0;
        double c = (1 - pen)*(prior[0] + prior[1]*e1);
        array<double, 2> up_o = bo.upper_node == a.root.get() ? array<double, 2>{0, 0} : above(msg, kd_o);
        array<double, 2> up_n = bn.upper_node == a.root.get() ? array<double, 2>{0, 0} : above(msg, kd_n);
        double w_o = join(bo, ro.upper_node->time, query, msg[kq], msg[kd_o], up_o);
        double w_n = join(bn, rn.upper_node->time, query, msg[kq], msg[kd_n], up_n);
        return log((w_n - c)/(w_o - c));
    };
    auto apply_pattern = [&]() {
        if (!pattern_valid) {
            pat = ref;
            for (int s : missing_list) {
                int k = t.leaf[s];
                pat[k] = {1.0, 1.0};
                a.repeel_path(flat, t, unit, k, pat);
            }
            pattern_valid = true;
            pattern_value_valid = false;
        }
        for (int s : touched) {
            if (target[s] != missing[s]) {
                missing[s] = target[s];
                int k = t.leaf[s];
                pat[k] = missing[s] ? array<double, 2>{1.0, 1.0} : array<double, 2>{1.0, 0.0};
                (main[k] ? observed_main : observed_sub) += missing[s] ? -1 : 1;
                a.repeel_path(flat, t, unit, k, pat);
                pattern_value_valid = false;
            }
        }
        if (!pattern_value_valid) {
            pattern_value = pattern_log_ratio(pat);
            pattern_value_valid = true;
        }
    };
    auto entry_log_ratio = [&](const Genotype_likelihood &g, double pos) {
        Node *query = ro.lower_node;
        Peel &p = ratio_peel;
        vector<array<double, 2>> &msg = ratio_cell;
        a.peel_site(flat, t, unit, pos, p);
        int nc = g.leaves[1] == nullptr ? 2 : 4;
        int kl[2] = {t.index_of(flat, g.leaves[0]), g.leaves[1] == nullptr ? -1 : t.index_of(flat, g.leaves[1])};
        int kq = query->flat_index, kd_o = bo.lower_node->flat_index, kd_n = bn.lower_node->flat_index;
        double w_o = 0, w_n = 0;
        for (int c = 0; c < nc; c++) {
            if (g.L[c] == 0) {
                continue;
            }
            msg = p.msg;
            int count[2] = {p.count[0], p.count[1]};
            for (int j = 0; j < 2 and kl[j] >= 0; j++) {
                if (!p.missing[kl[j]]) {
                    int h = j == 0 ? (c & 1) : (c >> 1);
                    count[(int) g.leaves[j]->get_state(pos)] -= 1;
                    count[h] += 1;
                    msg[kl[j]] = {h == 0 ? 1.0 : 0.0, h == 1 ? 1.0 : 0.0};
                    a.repeel_path(flat, t, unit, kl[j], msg);
                }
            }
            double empty = (1 - pen)*(prior[0]*(count[1] == 0) + prior[1]*(count[0] == 0));
            array<double, 2> up_o = bo.upper_node == a.root.get() ? array<double, 2>{0, 0} : above(msg, kd_o);
            array<double, 2> up_n = bn.upper_node == a.root.get() ? array<double, 2>{0, 0} : above(msg, kd_n);
            w_o += g.L[c]*(join(bo, ro.upper_node->time, query, msg[kq], msg[kd_o], up_o) - empty);
            w_n += g.L[c]*(join(bn, rn.upper_node->time, query, msg[kq], msg[kd_n], up_n) - empty);
        }
        return log(w_n/w_o);
    };

    double ll = 0;
    for (int i = start_index; i < end_index; i++) {
        double x = a.coordinates[i];
        double w = a.coordinates[i + 1] - x;
        bool moved = false;
        while (oj != a.joining_branches.end() and oj->first <= x) {
            bo = oj->second;
            ++oj;
            moved = true;
        }
        while (on != a.removed_branches.end() and on->first <= x) {
            ro = on->second;
            ++on;
            moved = true;
        }
        while (nj != new_joining_branches.end() and nj->first <= x) {
            bn = nj->second;
            ++nj;
            moved = true;
        }
        while (nn != added_branches.end() and nn->first <= x) {
            rn = nn->second;
            ++nn;
            moved = true;
        }
        if (recomb_it->first == x) {
            flat.forward_update(recomb_it->second);
            ++recomb_it;
            order_stale = true;
        }
        double c = pen*a.thetas[i]/w;
        if (order_stale or c != unit) {
            if (order_stale) {
                t.assign(flat, a.root.get());
                main.assign(flat.parents.size(), 0);
                main[t.top] = 1;
                for (int k = (int) flat.parents.size() - 1; k >= 0; k--) {
                    if (k != t.top and t.up[k] >= 0 and main[t.up[k]]) {
                        main[k] = 1;
                    }
                }
                order_stale = false;
            }
            unit = c;
            a.reference_peel(flat, t, unit, ref);
            p_common = 1;
            observed_main = observed_sub = 0;
            for (int k = 0; k < (int) flat.parents.size(); k++) {
                Node *l = flat.parents[k].first;
                if (l->is_sample) {
                    (main[k] ? observed_main : observed_sub) += 1;
                }
                if (k != t.top and flat.parents[k].second != a.cut_node) {
                    p_common *= 1 + unit*flat.lengths[k];
                }
            }
            for (int s : missing_list) {
                (main[t.leaf[s]] ? observed_main : observed_sub) -= 1;
            }
            ref_valid = pattern_valid = z_valid = false;
        }
        if (moved) {
            ref_valid = pattern_value_valid = z_valid = false;
        }
        if (!z_valid) {
            z_diff = log_z(bn, rn.upper_node->time, rn.lower_node) - log_z(bo, ro.upper_node->time, ro.lower_node);
            z_valid = true;
        }
        ll -= a.assayed[i]*z_diff;
        if (!ref_valid) {
            int all_main = observed_main, all_sub = observed_sub;
            for (int s : missing_list) {
                (main[t.leaf[s]] ? all_main : all_sub) += 1;
            }
            swap(observed_main, all_main);
            swap(observed_sub, all_sub);
            ref_value = pattern_log_ratio(ref);
            swap(observed_main, all_main);
            swap(observed_sub, all_sub);
            ref_valid = true;
        }
        double rest = a.assayed[i];
        if (a.any_missing) {
            const Bin_gaps &g = owner->gaps_of(i);
            if (g.uniform and g.full == missing_list) {
                touched.clear();
                apply_pattern();
                ll += g.count*pattern_value;
                rest -= g.count;
            } else if (!g.sample.empty()) {
                touched = missing_list;
                for (int s : touched) {
                    target[s] = 0;
                    seen[s] = 1;
                }
                for (size_t j = 0; j < g.sample.size(); j++) {
                    if (g.before[j] > 0) {
                        apply_pattern();
                        ll += g.before[j]*pattern_value;
                        rest -= g.before[j];
                    }
                    int s = g.sample[j];
                    target[s] = g.d[j] < 0;
                    if (!seen[s]) {
                        seen[s] = 1;
                        touched.push_back(s);
                    }
                }
                if (g.tail > 0) {
                    apply_pattern();
                    ll += g.tail*pattern_value;
                    rest -= g.tail;
                }
                missing_list.clear();
                for (int s : touched) {
                    seen[s] = 0;
                    if (missing[s]) {
                        missing_list.push_back(s);
                    }
                }
                sort(missing_list.begin(), missing_list.end());
            }
        }
        while (*mut_it < a.coordinates[i + 1]) {
            while (qi != qe and gl[*qi].pos < *mut_it) {
                qi++;
            }
            if (qi != qe and gl[*qi].pos == *mut_it) {
                ll += entry_log_ratio(gl[*qi], *mut_it);
            } else {
                const Site_messages &s = site_pool[site_slot.at(*mut_it)];
                Node *query = ro.lower_node;
                int kd_o = bo.lower_node->flat_index, kd_n = bn.lower_node->flat_index;
                double w_o = join(bo, ro.upper_node->time, query, s.nu, s.in[kd_o], s.above[kd_o]);
                double w_n = join(bn, rn.upper_node->time, query, s.nu, s.in[kd_n], s.above[kd_n]);
                ll += log((w_n - s.empty_term)/(w_o - s.empty_term));
            }
            rest -= 1;
            ++mut_it;
        }
        ll += rest*ref_value;
    }
    return ll;
}

// the proposal reads the moved individual's genotype cells, which the ratio sums out, so they are redrawn on the proposed ARG before the decision and the reverse densities read the redrawn cells
double Threader_smc::redrawn_acceptance_ratio(ARG &a, tuple<double, Branch, double> cut_point) {
    double jac = cut_ratio(a);
    double q_new = bsp.branch_log_q(new_joining_branches, start_index, a.coordinates);
    double h_new = tsp.sel_log_q + tsp.time_log_q;
    if (!(isfinite(q_new) and isfinite(h_new))) {
        return 0;
    }
    double log_pi = log_pi_ratio(a);
    map<double, Branch> old_jb = a.joining_branches, old_rb = a.removed_branches;
    Node *query = old_rb.begin()->second.lower_node;
    a.add(new_joining_branches, added_branches);
    a.smc_sample_recombinations(added_branches);
    a.resample_entries(query, start, end);
    redrawn_records.clear();
    redrawn_records.insert(a.recombinations.lower_bound(start), a.recombinations.upper_bound(end));
    a.clear_remove_info();
    a.cut_tree = a.get_tree_at(get<0>(cut_point));
    a.remove({get<0>(cut_point), prev(added_branches.upper_bound(get<0>(cut_point)))->second, cut_time});
    a.joining_branches = old_jb;
    a.removed_branches = old_rb;
    bsp.reset();
    tsp.reset();
    site_slot.clear();
    emit_trees = 0;
    set_check_points(a);
    a.set_collapsed(query, old_rb.begin()->first, old_rb.rbegin()->first);
    run_BSP(a);
    double q_old = bsp.branch_log_q(old_jb, start_index, a.coordinates);
    run_TSP(a, old_jb);
    double h_old = tsp.eval_joining_nodes(old_jb, old_rb, start_index, a.coordinates);
    return exp(log_pi + (q_old - q_new) + (h_old - h_new) + log(jac));
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

bool Threader_smc::joins_deleted_node(ARG &a, map<double, Branch> &lineage) {
    map<Node *, double> deleted_at;
    for (auto it = a.recombinations.upper_bound(a.start); it != a.recombinations.end() and it->first <= a.end; ++it) {
        Recombination &r = it->second;
        auto d = deleted_at.find(r.inserted_node);
        if (d != deleted_at.end()) {
            if (joined_at(lineage, d->second, it->first, r.inserted_node)) {
                return true;
            }
            deleted_at.erase(d);
        }
        deleted_at[r.deleted_node] = it->first;
    }
    return false;
}

bool Threader_smc::keeps_deleted_node(ARG &a) {
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

void Threader_smc::set_emission_sites(ARG &a, Emission &e, bool &stale, int bin, vector<pair<double, int>> &mut_set, Node *query_node, const Tree *tree) {
    double unit = a.penalty*a.thetas[bin]/(a.coordinates[bin + 1] - a.coordinates[bin]);
    for (auto &m : mut_set) {
        if (site_slot.count(m.first) == 0) {
            if (stale) {
                if (tree != nullptr) {
                    emit_flat.assign(*tree);
                }
                emit_order.assign(emit_flat, a.root.get());
                if (emit_trees == (int) emit_nodes.size()) {
                    emit_nodes.emplace_back();
                }
                vector<Node *> &nodes = emit_nodes[emit_trees++];
                nodes.clear();
                for (auto &x : emit_flat.parents) {
                    nodes.push_back(x.first);
                }
                stale = false;
            }
            int slot = (int) site_slot.size();
            if (slot == (int) site_pool.size()) {
                site_pool.emplace_back();
            }
            a.emission_site(emit_flat, emit_order, unit, m.first, query_node, emit_peel, site_pool[slot]);
            site_pool[slot].tree = emit_trees - 1;
            site_slot[m.first] = slot;
        }
    }
    e.sites.clear();
    for (auto &m : mut_set) {
        e.sites.push_back(&site_pool[site_slot[m.first]]);
    }
}

void Threader_smc::run_BSP(ARG &a) {
    bsp.reserve_memory(end_index - start_index);
    bsp.set_cutoff(cutoff);
    bsp.set_emission(pe);
    bsp.start(a.start_tree, cut_time);
    auto recomb_it = a.recombinations.upper_bound(start);
    auto mut_it = a.site_at(start);
    auto query_it = a.removed_branches.begin();
    vector<double> mutations;
    vector<pair<double, int>> mut_set = {};
    set<Branch> deletions = {};
    set<Branch> insertions = {};
    Node *query_node = nullptr;
    bool varying = varying_rate(a, start_index, end_index);
    Tree tree;
    if (varying) {
        tree = a.start_tree;
    }
    emit_flat.assign(a.start_tree);
    bool stale = true;
    double curr_rho = a.thetas[start_index]/(a.coordinates[start_index + 1] - a.coordinates[start_index]);
    double p_tree = branch_product(a.start_tree, a.penalty, curr_rho, a.cut_node);
    for (int i = start_index; i < end_index; i++) {
        if (a.coordinates[i] == query_it->first) {
            query_node = query_it->second.lower_node;
            query_it++;
        }
        if (a.coordinates[i] == recomb_it->first) {
            Recombination &r = recomb_it->second;
            recomb_it++;
            bsp.transfer(r);
            emit_flat.forward_update(r);
            stale = true;
            if (varying) {
                tree.forward_update(r);
                curr_rho = -1;
            } else {
                p_tree = update_branch_product(p_tree, r, a.penalty*curr_rho, a.cut_node);
            }
        } else if (a.coordinates[i] != start) {
            bsp.forward(a.rhos[i - 1]);
        }
        double w = a.coordinates[i + 1] - a.coordinates[i];
        double rho = a.thetas[i]/w;
        if (varying and rho != curr_rho) {
            p_tree = branch_product(tree, a.penalty, rho, a.cut_node);
            curr_rho = rho;
        }
        pe->set_tree_product(a.penalty, no_data ? 0.0 : rho, p_tree, called_width(query_node, a, i));
        mut_set.clear();
        while (mut_it->first < a.coordinates[i + 1]) {
            mut_set.push_back(*mut_it);
            mut_it++;
        }
        if (mut_set.size() > 0 and !no_data) {
            set_emission_sites(a, *pe, stale, i, mut_set, query_node, nullptr);
            bsp.mut_emit(query_node);
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

void Threader_smc::find_resume_index(ARG &a) {
    auto new_it = new_joining_branches.upper_bound(start);
    auto old_it = a.joining_branches.upper_bound(start);
    Branch new_branch = new_joining_branches.begin()->second;
    Branch old_branch = a.joining_branches.begin()->second;
    for (int i = start_index; i < end_index; i++) {
        auto old_at = old_it;
        if (a.coordinates[i] == new_it->first) {
            new_branch = new_it->second;
            new_it++;
        }
        if (a.coordinates[i] == old_it->first) {
            old_branch = old_it->second;
            old_it++;
        }
        if (new_branch != old_branch) {
            if (i > start_index) {
                resume_index = i;
                resume_join_it = old_at;
            }
            return;
        }
    }
}

void Threader_smc::run_TSP(ARG &a, map<double, Branch> &jb, bool resume) {
    tsp.reserve_memory(end_index - start_index);
    tsp.set_gap(gap);
    tsp.set_emission(be);
    map<double, Recombination>::iterator recomb_it;
    map<double, Branch>::iterator join_it;
    vector<pair<double, int>>::iterator mut_it;
    map<double, Branch>::iterator query_it;
    Branch prev_branch;
    Branch next_branch;
    Node *query_node = nullptr;
    vector<pair<double, int>> mut_set = {};
    bool varying;
    Tree &tree = tsp_tree;
    RSP_smc &rsp = tsp_rsp;
    auto set_sister_mass = [&](Branch &b) {
        double lo = b.lower_node->time;
        if (lo < cut_time) {
            rsp.set_tree(tree);
            tsp.sister_mass = rsp.sister_mass(lo, cut_time);
        } else {
            tsp.sister_mass = 0;
        }
    };
    double curr_rho;
    double p_tree;
    int first = start_index;
    if (resume) {
        tsp.restore_state(saved_tsp);
        recomb_it = saved_recomb_it;
        join_it = resume_join_it;
        mut_it = saved_mut_it;
        query_it = saved_query_it;
        prev_branch = saved_prev_branch;
        next_branch = saved_next_branch;
        query_node = saved_query_node;
        varying = varying_rate(a, start_index, end_index);
        tree = saved_tree;
        curr_rho = saved_curr_rho;
        p_tree = saved_p_tree;
        first = resume_index;
    } else {
        tsp.cc = make_shared<coalescent_calculator>(cut_time);
        tsp.cc->start(a.start_tree);
        tsp.cc1 = make_shared<coalescent_calculator>(cut_time);
        tsp.cc1->extra = 1;
        tsp.cc1->start(a.start_tree);
        Branch start_branch = jb.begin()->second;
        tsp.start(start_branch, cut_time);
        recomb_it = a.recombinations.upper_bound(start);
        join_it = jb.upper_bound(start);
        mut_it = a.site_at(start);
        query_it = a.removed_branches.lower_bound(start);
        prev_branch = start_branch;
        next_branch = start_branch;
        varying = varying_rate(a, start_index, end_index);
        tree = a.start_tree;
        set_sister_mass(start_branch);
        curr_rho = a.thetas[start_index]/(a.coordinates[start_index + 1] - a.coordinates[start_index]);
        p_tree = branch_product(a.start_tree, a.penalty, curr_rho, a.cut_node);
    }
    bool stale = true;
    for (int i = first; i < end_index; i++) {
        if (i == resume_index and !resume) {
            tsp.save_state(saved_tsp);
            saved_recomb_it = recomb_it;
            saved_mut_it = mut_it;
            saved_query_it = query_it;
            saved_prev_branch = prev_branch;
            saved_next_branch = next_branch;
            saved_query_node = query_node;
            saved_tree = tree;
            saved_curr_rho = curr_rho;
            saved_p_tree = p_tree;
        }
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
            stale = true;
            if (!varying) {
                p_tree = update_branch_product(p_tree, r, a.penalty*curr_rho, a.cut_node);
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
            p_tree = branch_product(tree, a.penalty, rho, a.cut_node);
            curr_rho = rho;
        }
        be->set_tree_product(a.penalty, no_data ? 0.0 : rho, p_tree, called_width(query_node, a, i));
        mut_set.clear();
        while (mut_it->first < a.coordinates[i+1]) {
            mut_set.push_back(*mut_it);
            mut_it++;
        }
        if (mut_set.size() > 0 and !no_data) {
            set_emission_sites(a, *be, stale, i, mut_set, query_node, &tree);
            if (stale) {
                Node *n = next_branch.lower_node;
                const vector<Node *> &nodes = emit_nodes[be->sites[0]->tree];
                n->flat_index = (int) (lower_bound(nodes.begin(), nodes.end(), n, compare_node()) - nodes.begin());
            }
            tsp.mut_emit(query_node);
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
        if (resume_index >= 0) {
            run_TSP(a, a.joining_branches, true);
        } else {
            set<double> check_points = tsp.check_points;
            tsp.reset();
            tsp.check_points = check_points;
            run_TSP(a, a.joining_branches);
        }
    }
    double h_old = tsp.eval_joining_nodes(a.joining_branches, a.removed_branches, start_index, a.coordinates);
    double log_a = log_pi_ratio(a) +(q_old - q_new) + (h_old - h_new) + log(jac);
    return isfinite(q_new) and isfinite(h_new) ? exp(log_a) : 0;
}

double Threader_smc::random() {
    return uniform_random();
}

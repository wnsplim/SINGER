//
//  Threader_smc.cpp
//  SINGER
//
//  Created by Yun Deng on 4/4/23.
//

#include "Threader_smc.hpp"

bool Threader_smc::no_data = false;

static void check_drawn_density(double q_new, double h_new) {
    if (!(isfinite(q_new) and isfinite(h_new))) {
        cerr << "Error: the proposal density of the drawn path is not finite (branch " << q_new << ", time " << h_new << ")." << endl;
        exit(1);
    }
}

static void check_ratio(double ar) {
    if (isnan(ar)) {
        cerr << "Error: the acceptance ratio of a proposal is not a number." << endl;
        exit(1);
    }
}

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
    if (keeps_deleted_node(a) and !has_gap_node(added_branches, a)) {
        a.add(new_joining_branches, added_branches);
    } else {
        a.add(a.joining_branches, a.removed_branches);
    }
    a.approx_sample_recombinations();
    a.clear_remove_info();
}

void Threader_smc::internal_rethread(ARG &a, tuple<double, Branch, double> cut_point) {
    cut_time = get<2>(cut_point);
    a.remove(cut_point);
    if (a.augment) {
        a.detach_histories();
    }
    get_boundary(a);
    set_check_points(a);
    run_BSP(a);
    sample_joining_branches(a);
    run_TSP(a);
    sample_joining_points(a);
    double ar = has_gap_node(added_branches, a) ? 0 : acceptance_ratio(a);
    check_ratio(ar);
    double q = random();
    if (q < ar) {
        if (a.augment) {
            a.attach_histories(new_joining_branches, added_branches);
            a.keep_histories();
        }
        a.add(new_joining_branches, added_branches);
    } else {
        a.add(a.joining_branches, a.removed_branches);
        if (a.augment) {
            a.restore_histories();
        }
    }
    a.approx_sample_recombinations();
    a.clear_remove_info();
}

void Threader_smc::exact_internal_rethread(ARG &a, tuple<double, Branch, double> cut_point) {
    cut_time = get<2>(cut_point);
    a.remove(cut_point);
    if (a.augment) {
        a.detach_histories();
    }
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
        if (!has_gap_node(added_branches, a) and !joins_deleted_node(a, added_branches) and !a.empties_cut_span_end_recombination(new_joining_branches, added_branches)) {
            ar = entries ? redrawn_acceptance_ratio(a, cut_point) : exact_acceptance_ratio(a);
        }
    }
    check_ratio(ar);
    a.journal_on = false;
    double q = random();
    if (q < ar and entries) {
        a.add(new_joining_branches, added_branches);
        a.recombinations.erase(a.recombinations.lower_bound(start), a.recombinations.upper_bound(end));
        a.recombinations.insert(redrawn_records.begin(), redrawn_records.end());
    } else if (q < ar) {
        if (a.augment) {
            a.attach_histories(new_joining_branches, added_branches);
            a.keep_histories();
        }
        a.add(new_joining_branches, added_branches);
        a.smc_sample_recombinations(added_branches);
    } else {
        a.undo_journal();
        a.add(a.joining_branches, a.removed_branches);
        a.smc_sample_recombinations(a.removed_branches);
        if (a.augment) {
            a.restore_histories();
        }
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
    if (a.augment) {
        return pi_new - pi_old + history_log_ratio(a);
    }
    return pi_new - pi_old + data_log_ratio(a);
}

double Threader_smc::history_log_ratio(ARG &a) {
    Flat_tree &flat = ratio_flat;
    Tree_order &t = ratio_order;
    flat.assign(a.start_tree);
    vector<int> &below = ratio_below;
    vector<Pattern_group> &groups = ratio_groups;
    auto recomb_it = a.recombinations.upper_bound(start);
    auto mut_it = a.mutation_sites.lower_bound(start);
    auto oj = a.joining_branches.begin(), on = a.removed_branches.begin();
    auto nj = new_joining_branches.begin(), nn = added_branches.begin();
    Branch bo, bn, ro, rn;
    bool stale = true, unknown = a.unknown_base();
    int S = a.model.states;
    long double ll = 0;
    double w4[max_states];
    auto node_log_value = [&](const Branch &b, double time, double t_d, int xu, int xl, unsigned mask, double u) {
        return logl(a.node_weights(b, time, t_d, xu, xl, mask, u, w4)) - logl(a.model.transition_entry(u*(b.upper_node->time - b.lower_node->time), xu, xl));
    };
    for (int i = start_index; i < end_index; i++) {
        double x = a.coordinates[i];
        while (oj != a.joining_branches.end() and oj->first <= x) {
            bo = oj->second;
            ++oj;
        }
        while (on != a.removed_branches.end() and on->first <= x) {
            ro = on->second;
            ++on;
        }
        while (nj != new_joining_branches.end() and nj->first <= x) {
            bn = nj->second;
            ++nj;
        }
        while (nn != added_branches.end() and nn->first <= x) {
            rn = nn->second;
            ++nn;
        }
        if (recomb_it->first == x) {
            flat.forward_update(recomb_it->second);
            ++recomb_it;
            stale = true;
        }
        if (stale) {
            t.assign(flat, a.root.get());
            stale = false;
        }
        double unit = a.penalty*a.thetas[i]/(a.coordinates[i + 1] - x);
        Node *query = ro.lower_node;
        double to = ro.upper_node->time, tn = rn.upper_node->time;
        a.pattern_groups(i, groups);
        for (Pattern_group &g : groups) {
            a.column_below(flat, t, x, &g.forced, below);
            bool data = query->is_sample ? !g.forced[query->index] : below[query->flat_index] > 0;
            if (!data or below[t.top] == 0) {
                continue;
            }
            Join_case co = a.join_case(flat, t, below, bo, to, query, true, unit), cn = a.join_case(flat, t, below, bn, tn, query, true, unit);
            const Joint_terms &J = a.joint_terms(cn, co, bn, tn, bo, to, unit*a.data_length(flat, t, below), unit);
            for (int b = unknown ? -1 : 0; b < (unknown ? 0 : S); b++) {
                double n = b < 0 ? g.unknown : g.count[b];
                if (n > 0) {
                    ll += n*(a.joint_log_mix(J, J.F, b) - a.joint_log_mix(J, J.R, b));
                }
            }
        }
        auto column = [&](double pos, int base) {
            a.column_below(flat, t, pos, nullptr, below);
            unsigned mask = a.query_mask(query, below, pos, base);
            if (mask == 0) {
                return;
            }
            double u = below[t.top] > 0 ? unit*a.column_rate(pos) : unit;
            Join_case co = a.join_case(flat, t, below, bo, to, query, true, u), cn = a.join_case(flat, t, below, bn, tn, query, true, u);
            if (co.anchor < 0) {
                ll += logl(a.join_value(cn, -1, mask)) - logl(a.join_value(co, -1, mask));
                return;
            }
            const vector<Mutation_row> &rows = a.mutations.at(pos);
            if (base >= 0 and a.only_replaced_changes(rows, base, cn, co, bn, bo, pos)) {
                const Joint_terms &J = a.joint_terms(a.join_case(flat, t, below, bn, tn, query, true, unit), a.join_case(flat, t, below, bo, to, query, true, unit), bn, tn, bo, to, unit*a.data_length(flat, t, below), unit);
                int b = unknown ? -1 : base;
                ll += a.joint_log_mix(J, J.F, b) - a.joint_log_mix(J, J.R, b);
                return;
            }
            auto log_value = [&](const Join_case &c, const Branch &br, double time) {
                if (c.kind == E1) {
                    int xu = a.node_state(br.upper_node, pos);
                    return node_log_value(br, time, query->time, xu, a.branch_state(rows, br, br.lower_node->time, pos), mask, u);
                }
                return logl(a.join_value(c, a.join_state(c, flat, rows, br, time, pos), mask));
            };
            ll += log_value(cn, bn, tn) - log_value(co, bo, to);
        };
        while (*mut_it < a.coordinates[i + 1]) {
            column(*mut_it, -1);
            ++mut_it;
        }
        for (auto &e : a.exception_sites[i]) {
            column(e.first, e.second);
        }
    }
    return ll;
}

double Threader_smc::data_log_ratio(ARG &a) {
    Flat_tree &flat = ratio_flat;
    Tree_order &t = ratio_order;
    vector<vector<Msg>> &ref = ratio_ref, &pat = ratio_pattern;
    flat.assign(a.start_tree);
    const Site_model &model = a.model;
    int S = model.states, nbase = S == 2 ? 1 : S, K = (int) model.rates.size();
    bool product = model.branch_product;
    ref.resize(nbase*K);
    pat.resize(nbase*K);
    Msg prior = {};
    a.root_prior(prior.data());
    double pen = a.penalty;
    double T[max_states*max_states];
    Msg zero = {};
    auto [qi, qe] = a.entries_in(a.removed_branches.begin()->second.lower_node, start, end);
    const vector<Genotype_likelihood> &gl = a.genotype_likelihoods;
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
    array<double, max_states> ref_value = {}, pattern_value = {};
    double ref_value_u = 0, pattern_value_u = 0;
    bool unknown = a.unknown_base();
    double z_diff = 0;

    auto above = [&](const vector<Msg> &msg, int d, double u) {
        path.clear();
        for (int k = d; k != t.top; k = t.up[k]) {
            path.push_back(k);
        }
        Msg A = prior;
        int v = t.top;
        for (int j = (int) path.size() - 1; j >= 0; j--) {
            int c = path[j];
            Msg B = A;
            for (int sib : t.children[v]) {
                if (sib != c) {
                    model.transition(u*flat.lengths[sib], T);
                    push_up(T, S, product, msg[sib], B);
                }
            }
            if (j == 0) {
                return B;
            }
            model.transition(u*flat.lengths[c], T);
            transmit_up(T, S, product, B, A);
            v = c;
        }
        return A;
    };
    auto join = [&](const Branch &b, double time, Node *query, const Msg &nu, const Msg &m_d, const Msg &up, double u) {
        Msg d, n, m, v;
        model.transition(u*(time - b.lower_node->time), T);
        transmit(T, S, product, m_d, d);
        model.transition(u*(time - query->time), T);
        transmit(T, S, product, nu, n);
        for (int q = 0; q < S; q++) {
            m[q] = d[q]*n[q];
        }
        double w = 0;
        if (b.upper_node == a.root.get()) {
            for (int q = 0; q < S; q++) {
                w += prior[q]*m[q];
            }
            return w;
        }
        model.transition(u*(b.upper_node->time - time), T);
        transmit(T, S, product, m, v);
        for (int q = 0; q < S; q++) {
            w += up[q]*v[q];
        }
        return w;
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
    auto pattern_w = [&](const vector<vector<Msg>> &set, int b, double &w_o, double &w_n) {
        Node *query = ro.lower_node;
        int kq = query->flat_index, kd_o = bo.lower_node->flat_index, kd_n = bn.lower_node->flat_index;
        w_o = w_n = 0;
        for (int k = 0; k < K; k++) {
            const vector<Msg> &msg = set[b*K + k];
            double u = unit*model.rates[k];
            Msg up_o = bo.upper_node == a.root.get() ? zero : above(msg, kd_o, u);
            Msg up_n = bn.upper_node == a.root.get() ? zero : above(msg, kd_n, u);
            w_o += join(bo, ro.upper_node->time, query, msg[kq], msg[kd_o], up_o, u);
            w_n += join(bn, rn.upper_node->time, query, msg[kq], msg[kd_n], up_n, u);
        }
    };
    auto pattern_log_ratio = [&](const vector<vector<Msg>> &set, int b) {
        double e1 = observed_main == 0 and observed_sub == 0 ? 1.0 : 0.0;
        double c = product ? (1 - pen)*(prior[0] + prior[1]*e1) : 0.0;
        double w_o, w_n;
        pattern_w(set, b, w_o, w_n);
        return log((w_n - c)/(w_o - c));
    };
    auto pattern_log_ratio_unknown = [&](const vector<vector<Msg>> &set) {
        double sum_o = 0, sum_n = 0;
        for (int b = 0; b < S; b++) {
            double w_o, w_n;
            pattern_w(set, b, w_o, w_n);
            sum_o += w_o;
            sum_n += w_n;
        }
        return log(sum_n/sum_o);
    };
    auto apply_pattern = [&]() {
        if (!pattern_valid) {
            for (int i = 0; i < nbase*K; i++) {
                pat[i] = ref[i];
                for (int s : missing_list) {
                    int k = t.leaf[s];
                    pat[i][k] = ones;
                    a.repeel_path(flat, t, unit*model.rates[i % K], k, pat[i]);
                }
            }
            pattern_valid = true;
            pattern_value_valid = false;
        }
        for (int s : touched) {
            if (target[s] != missing[s]) {
                missing[s] = target[s];
                int k = t.leaf[s];
                for (int i = 0; i < nbase*K; i++) {
                    pat[i][k] = missing[s] ? ones : a.leaf_message(product ? 0 : 1 << (i/K));
                    a.repeel_path(flat, t, unit*model.rates[i % K], k, pat[i]);
                }
                (main[k] ? observed_main : observed_sub) += missing[s] ? -1 : 1;
                pattern_value_valid = false;
            }
        }
        if (!pattern_value_valid) {
            for (int b = 0; b < nbase; b++) {
                pattern_value[b] = pattern_log_ratio(pat, b);
            }
            if (unknown) {
                pattern_value_u = pattern_log_ratio_unknown(pat);
            }
            pattern_value_valid = true;
        }
    };
    auto entry_log_ratio = [&](const Genotype_likelihood &g, double pos) {
        Node *query = ro.lower_node;
        Peel &p = ratio_peel;
        vector<Msg> &msg = ratio_cell;
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
                    msg[kl[j]] = a.leaf_message(h);
                    a.repeel_path(flat, t, unit, kl[j], msg);
                }
            }
            double empty = (1 - pen)*(prior[0]*(count[1] == 0) + prior[1]*(count[0] == 0));
            Msg up_o = bo.upper_node == a.root.get() ? zero : above(msg, kd_o, unit);
            Msg up_n = bn.upper_node == a.root.get() ? zero : above(msg, kd_n, unit);
            w_o += g.L[c]*(join(bo, ro.upper_node->time, query, msg[kq], msg[kd_o], up_o, unit) - empty);
            w_n += g.L[c]*(join(bn, rn.upper_node->time, query, msg[kq], msg[kd_n], up_n, unit) - empty);
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
            for (int b = 0; b < nbase; b++) {
                for (int k = 0; k < K; k++) {
                    a.reference_peel(flat, t, unit*model.rates[k], ref[b*K + k], b);
                }
            }
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
            z_diff = product ? log_z(bn, rn.upper_node->time, rn.lower_node) - log_z(bo, ro.upper_node->time, ro.lower_node) : 0.0;
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
            for (int b = 0; b < nbase; b++) {
                ref_value[b] = pattern_log_ratio(ref, b);
            }
            if (unknown) {
                ref_value_u = pattern_log_ratio_unknown(ref);
            }
            swap(observed_main, all_main);
            swap(observed_sub, all_sub);
            ref_valid = true;
        }
        double rest = a.assayed[i];
        array<double, max_states> rest_b = {};
        const Bin_gaps *gp = (a.any_missing or S != 2) ? &a.gaps_of(i) : nullptr;
        if (S != 2) {
            rest_b = gp->all_bases;
        }
        auto score = [&](double count, const array<double, max_states> &bases) {
            apply_pattern();
            if (S == 2) {
                ll += count*pattern_value[0];
                rest -= count;
            } else if (unknown) {
                ll += count*pattern_value_u;
                rest -= count;
            } else {
                for (int b = 0; b < S; b++) {
                    ll += bases[b]*pattern_value[b];
                    rest_b[b] -= bases[b];
                }
            }
        };
        if (a.any_missing) {
            const Bin_gaps &g = *gp;
            if (g.uniform and g.full == missing_list) {
                touched.clear();
                score(g.count, g.count_bases);
            } else if (!g.sample.empty()) {
                touched = missing_list;
                for (int s : touched) {
                    target[s] = 0;
                    seen[s] = 1;
                }
                for (size_t j = 0; j < g.sample.size(); j++) {
                    if (g.before[j] > 0) {
                        score(g.before[j], g.bases[j]);
                    }
                    int s = g.sample[j];
                    target[s] = g.d[j] < 0;
                    if (!seen[s]) {
                        seen[s] = 1;
                        touched.push_back(s);
                    }
                }
                if (g.tail > 0) {
                    score(g.tail, g.tail_bases);
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
                double w_o = 0, w_n = 0;
                auto add = [&](const Site_messages &m) {
                    w_o += join(bo, ro.upper_node->time, query, m.nu, m.in[kd_o], m.above[kd_o], m.unit);
                    w_n += join(bn, rn.upper_node->time, query, m.nu, m.in[kd_n], m.above[kd_n], m.unit);
                };
                add(s);
                for (const Site_messages &part : s.parts) {
                    add(part);
                }
                ll += log((w_n - s.empty_term)/(w_o - s.empty_term));
            }
            rest -= 1;
            ++mut_it;
        }
        if (S == 2) {
            ll += rest*ref_value[0];
        } else if (unknown) {
            ll += rest*ref_value_u;
        } else {
            for (int b = 0; b < S; b++) {
                ll += rest_b[b]*ref_value[b];
            }
        }
    }
    return ll;
}

// the proposal reads the moved individual's genotype cells, which the ratio sums out, so they are redrawn on the proposed ARG before the decision and the reverse densities read the redrawn cells
double Threader_smc::redrawn_acceptance_ratio(ARG &a, tuple<double, Branch, double> cut_point) {
    double jac = cut_ratio(a);
    double q_new = bsp.branch_log_q(new_joining_branches, start_index, a.coordinates);
    double h_new = tsp.sel_log_q + tsp.time_log_q;
    check_drawn_density(q_new, h_new);
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

void Threader_smc::prepare_emission_tree(ARG &a, bool &stale, const Tree *tree) {
    if (!stale) {
        return;
    }
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

void Threader_smc::set_emission_groups(ARG &a, Emission &e, bool &stale, int bin, Node *query_node, const Tree *tree) {
    e.groups.clear();
    e.group_weights.clear();
    if (a.model.branch_product) {
        return;
    }
    prepare_emission_tree(a, stale, tree);
    double unit = a.penalty*a.thetas[bin]/(a.coordinates[bin + 1] - a.coordinates[bin]);
    int S = a.model.states;
    e.unit = unit;
    if (a.augment) {
        vector<Pattern_group> &groups = emit_groups;
        a.pattern_groups(bin, groups);
        e.unknown = a.unknown_base();
        if (base_tree != emit_trees - 1 or base_query != query_node or base_unit != unit) {
            a.history_pattern(emit_flat, emit_order, unit, groups[0].forced, query_node, base_pool[0]);
            base_pool[0].tree = emit_trees - 1;
            base_tree = emit_trees - 1;
            base_query = query_node;
            base_unit = unit;
        }
        e.main_pattern = &base_pool[0];
        for (int b = 0; b < S; b++) {
            e.base_counts[b] = e.unknown ? 0.0 : groups[0].count[b];
        }
        e.unknown_count = e.unknown ? groups[0].unknown : 0.0;
        e.patterns.clear();
        if (group_pool.size() < groups.size()) {
            group_pool.resize(groups.size());
        }
        for (size_t j = 1; j < groups.size(); j++) {
            a.history_pattern(emit_flat, emit_order, unit, groups[j].forced, query_node, group_pool[j - 1]);
            array<double, max_states + 1> n = {};
            copy(groups[j].count.begin(), groups[j].count.begin() + S, n.begin());
            n[S] = groups[j].unknown;
            e.patterns.push_back({&group_pool[j - 1], n});
        }
        return;
    }
    const Bin_gaps &g = a.gaps_of(bin);
    int n = 0;
    for (Node *s : a.sample_nodes) {
        n = max(n, s->index + 1);
    }
    group_missing.assign(n, 0);
    if (base_tree != emit_trees - 1 or base_query != query_node or base_unit != unit) {
        for (int b = 0; b < S; b++) {
            a.emission_site(emit_flat, emit_order, unit, 0.0, query_node, emit_peel, base_pool[b], b, &group_missing);
            base_pool[b].tree = emit_trees - 1;
        }
        base_tree = emit_trees - 1;
        base_query = query_node;
        base_unit = unit;
    }
    for (int b = 0; b < S; b++) {
        e.base_groups[b] = &base_pool[b];
    }
    bool unknown = a.unknown_base();
    array<double, max_states> rest = g.all_bases;
    double rest_unknown = unknown ? g.all_count : 0.0;
    size_t used = 0;
    vector<size_t> unknown_slots;
    vector<double> unknown_weights;
    auto add_group = [&](const array<double, max_states> &n, double count) {
        if (unknown) {
            unknown_slots.push_back(used);
            unknown_weights.push_back(count);
        }
        for (int b = 0; b < S; b++) {
            if (!unknown and n[b] <= 0) {
                continue;
            }
            if (used == group_pool.size()) {
                group_pool.emplace_back();
            }
            a.emission_site(emit_flat, emit_order, unit, 0.0, query_node, emit_peel, group_pool[used], b, &group_missing);
            group_pool[used].tree = emit_trees - 1;
            if (!unknown) {
                e.group_weights.push_back(n[b]);
            }
            used++;
        }
    };
    for (size_t j = 0; j < g.sample.size(); j++) {
        if (g.before[j] > 0) {
            add_group(g.bases[j], g.before[j]);
            for (int b = 0; b < S; b++) {
                rest[b] -= g.bases[j][b];
            }
            rest_unknown -= g.before[j];
        }
        group_missing[g.sample[j]] = g.d[j] < 0;
    }
    if (g.tail > 0) {
        add_group(g.tail_bases, g.tail);
        for (int b = 0; b < S; b++) {
            rest[b] -= g.tail_bases[b];
        }
        rest_unknown -= g.tail;
    }
    for (int b = 0; b < S; b++) {
        e.base_counts[b] = rest[b];
    }
    e.unknown_count = unknown ? rest_unknown : 0.0;
    e.unknown_groups.clear();
    if (unknown) {
        for (size_t i = 0; i < unknown_slots.size(); i++) {
            array<const Site_messages *, max_states> set = {};
            for (int b = 0; b < S; b++) {
                set[b] = &group_pool[unknown_slots[i] + b];
            }
            e.unknown_groups.push_back({set, unknown_weights[i]});
        }
    } else {
        for (size_t i = 0; i < used; i++) {
            e.groups.push_back(&group_pool[i]);
        }
    }
}

void Threader_smc::set_emission_sites(ARG &a, Emission &e, bool &stale, int bin, vector<pair<double, int>> &mut_set, Node *query_node, const Tree *tree) {
    double unit = a.penalty*a.thetas[bin]/(a.coordinates[bin + 1] - a.coordinates[bin]);
    for (auto &m : mut_set) {
        if (site_slot.count(m.first) == 0) {
            prepare_emission_tree(a, stale, tree);
            int slot = (int) site_slot.size();
            if (slot == (int) site_pool.size()) {
                site_pool.emplace_back();
            }
            if (a.augment) {
                a.history_site(emit_flat, emit_order, unit, m.first, query_node, site_pool[slot]);
            } else {
                a.emission_site(emit_flat, emit_order, unit, m.first, query_node, emit_peel, site_pool[slot]);
            }
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
    pe->model = &a.model;
    pe->augment = a.augment;
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
        if (a.augment) {
            for (auto &e : a.exception_sites[i]) {
                mut_set.push_back({e.first, -1});
            }
        }
        set_emission_groups(a, *pe, stale, i, query_node, nullptr);
        if (mut_set.size() > 0 and !no_data) {
            set_emission_sites(a, *pe, stale, i, mut_set, query_node, nullptr);
        } else {
            pe->sites.clear();
        }
        if ((mut_set.size() > 0 or !pe->groups.empty() or !pe->unknown_groups.empty() or !pe->patterns.empty()) and !no_data) {
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
    be->model = &a.model;
    be->augment = a.augment;
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
        if (a.augment) {
            for (auto &e : a.exception_sites[i]) {
                mut_set.push_back({e.first, -1});
            }
        }
        set_emission_groups(a, *be, stale, i, query_node, &tree);
        if (mut_set.size() > 0 and !no_data) {
            set_emission_sites(a, *be, stale, i, mut_set, query_node, &tree);
            if (stale) {
                Node *n = next_branch.lower_node;
                const vector<Node *> &nodes = emit_nodes[be->sites[0]->tree];
                n->flat_index = (int) (lower_bound(nodes.begin(), nodes.end(), n, compare_node()) - nodes.begin());
            }
        } else {
            be->sites.clear();
        }
        if ((mut_set.size() > 0 or !be->groups.empty() or !be->unknown_groups.empty() or !be->patterns.empty()) and !no_data) {
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
    check_drawn_density(q_new, h_new);
    return exp(log_pi_ratio(a) + (q_old - q_new) + (h_old - h_new) + log(jac));
}

double Threader_smc::random() {
    return uniform_random();
}

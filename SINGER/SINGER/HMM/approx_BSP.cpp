//
//  approx_BSP.cpp
//  SINGER
//
//  Created by Yun Deng on 8/15/23.
//

#include "approx_BSP.hpp"

template <class V, class K>
static typename V::value_type::second_type &keyed_slot(V &v, const K &key) {
    auto it = lower_bound(v.begin(), v.end(), key,
                          [](const typename V::value_type &e, const K &k) { return e.first < k; });
    if (it != v.end() and !(key < it->first)) {
        return it->second;
    }
    return v.insert(it, {key, {}})->second;
}

template <class V, class K>
static typename V::value_type::second_type &recycled_slot(V &v, size_t &n, const K &key) {
    auto it = lower_bound(v.begin(), v.begin() + n, key,
                          [](const typename V::value_type &e, const K &k) { return e.first < k; });
    if (it != v.begin() + n and !(key < it->first)) {
        return it->second;
    }
    size_t p = it - v.begin();
    if (n == v.size()) {
        v.emplace_back();
    }
    rotate(v.begin() + p, v.begin() + n, v.begin() + n + 1);
    v[p].first = key;
    v[p].second.clear();
    n++;
    return v[p].second;
}

template <class V>
static typename V::value_type::second_type &keyed_lookup(V &v, int key) {
    auto it = upper_bound(v.begin(), v.end(), key,
                          [](int k, const typename V::value_type &e) { return k < e.first; });
    --it;
    return it->second;
}

approx_BSP::approx_BSP() {}

approx_BSP::~approx_BSP() {
    vector<vector<double>>().swap(forward_probs);
    vector<pair<int, vector<Interval_ptr>>>().swap(state_spaces);
    vector<pair<int, vector<double>>>().swap(times);
    vector<pair<int, vector<double>>>().swap(weights);
}

void approx_BSP::push_row(const vector<double> &v) {
    if (n_rows < forward_probs.size()) {
        vector<double> &row = forward_probs[n_rows];
        if (row.capacity() < v.size()) {
            row.reserve(max(v.size(), row.capacity() + row.capacity()/2));
        }
        row.assign(v.begin(), v.end());
    } else {
        forward_probs.push_back(v);
    }
    n_rows++;
}

void approx_BSP::reset() {
    n_rows = 0;
    curr_index = 0;
    state_spaces.assign(1, {INT_MAX, {}});
    times.assign(1, {INT_MAX, {}});
    weights.assign(1, {INT_MAX, {}});
    curr_intervals.clear();
    temp_intervals.clear();
    rhos.clear();
    recomb_sums.clear();
    weight_sums.clear();
    n_transfer_intervals = 0;
    n_transfer_weights = 0;
    temp.clear();
    trace_back_probs.clear();
    valid_branches.clear();
    check_points.clear();
    prev_rho = -1;
    prev_theta = -1;
    prev_node = nullptr;
    dim = 0;
    recomb_sum = 0;
    weight_sum = 0;
    sample_index = -1;
    states_change = false;
    cc = nullptr;
    n_arena = 0;
}

void approx_BSP::reserve_memory(int length) {
    forward_probs.reserve(length);
}

void approx_BSP::start(set<Branch> &branches, double t) {
    cut_time = t;
    curr_index = 0;
    for (Branch b : branches) {
        if (b.upper_node->time > cut_time) {
            valid_branches.insert(b);
        }
    }
    double lb = 0;
    double ub = 0;
    double p = 0;
    Interval_ptr new_interval = nullptr;
    cc = make_shared<approx_coalescent_calculator>(cut_time);
    cc->start(valid_branches);
    for (const Branch &b : branches) {
        if (b.upper_node->time > cut_time) {
            lb = max(b.lower_node->time, cut_time);
            ub = b.upper_node->time;
            p = cc->prob(lb, ub);
            new_interval = make_interval(b, lb, ub, curr_index);
            new_interval->source_pos = curr_index;
            curr_intervals.push_back(new_interval);
            temp.push_back(p);
        }
    }
    cutoff = min(0.01, cutoff/curr_intervals.size()); // adjust cutoff based on number of states;
    push_row(temp);
    weight_sums.push_back(0.0);
    set_dimensions();
    compute_interval_info();
    keyed_slot(state_spaces, curr_index) = curr_intervals;
    temp.clear();
}

void approx_BSP::start(Tree &tree, double t) {
    cut_time = t;
    curr_index = 0;
    for (auto &x : tree.parents) {
        if (x.second->time > cut_time) {
            valid_branches.insert(Branch(x.first, x.second));
        }
    }
    double lb = 0;
    double ub = 0;
    double p = 0;
    Interval_ptr new_interval = nullptr;
    cc = make_shared<approx_coalescent_calculator>(cut_time);
    cc->start(valid_branches);
    for (auto &x : tree.parents) {
        if (x.second->time > cut_time) {
            lb = max(x.first->time, cut_time);
            ub = x.second->time;
            p = cc->prob(lb, ub);
            new_interval = make_interval(Branch(x.first, x.second), lb, ub, curr_index);
            new_interval->source_pos = curr_index;
            curr_intervals.push_back(new_interval);
            temp.push_back(p);
        }
    }
    cutoff = min(0.01, cutoff/curr_intervals.size()); // adjust cutoff based on number of states;
    push_row(temp);
    weight_sums.push_back(0.0);
    set_dimensions();
    compute_interval_info();
    keyed_slot(state_spaces, curr_index) = curr_intervals;
    temp.clear();
}

void approx_BSP::set_cutoff(double x) {
    cutoff = x;
}

void approx_BSP::set_emission(shared_ptr<Emission> e) {
    eh = e;
}

void approx_BSP::set_check_points(set<double> &p) {
    check_points = p;
}

void approx_BSP::forward(double rho) {
    rhos.push_back(rho);
    compute_recomb_probs(rho);
    compute_recomb_weights(rho);
    prev_rho = rho;
    curr_index += 1;
    {
        const double *rp = recomb_probs.data();
        const double *fp = forward_probs[curr_index - 1].data();
        int n = (int) recomb_probs.size();
        double a0 = 0, a1 = 0, a2 = 0, a3 = 0;
        int k = 0;
        for (; k + 3 < n; k += 4) {
            a0 += rp[k]*fp[k];     a1 += rp[k+1]*fp[k+1];
            a2 += rp[k+2]*fp[k+2]; a3 += rp[k+3]*fp[k+3];
        }
        for (; k < n; k++) {
            a0 += rp[k]*fp[k];
        }
        recomb_sum = (a0 + a1) + (a2 + a3);
    }
    push_row(recomb_probs);
    double *curr_probs = forward_probs[curr_index].data();
    const double *prev_probs = forward_probs[curr_index - 1].data();
    const double *rprobs = recomb_probs.data();
    const double *rweights = recomb_weights.data();
    for (int i = 0; i < dim; i++) {
        curr_probs[i] = prev_probs[i]*(1 - rprobs[i]) + recomb_sum*rweights[i];
    }
    recomb_sums.push_back(recomb_sum);
    weight_sums.push_back(weight_sum);
}

void approx_BSP::transfer(Recombination &r) {
    rhos.push_back(0);
    prev_rho = -1;
    prev_theta = -1;
    recomb_sums.push_back(0);
    weight_sums.push_back(0);
    sanity_check(r);
    curr_index += 1;
    n_transfer_weights = 0;
    n_transfer_intervals = 0;
    temp.clear();
    temp_intervals.clear();
    for (int i = 0; i < curr_intervals.size(); i++) {
        process_interval(r, i);
    }
    add_new_branches(r);
    generate_intervals(r);
    set_dimensions();
    cc->update(r);
    compute_interval_info();
    keyed_slot(state_spaces, curr_index) = curr_intervals;
}

double approx_BSP::get_recomb_prob(double rho, double t) {
    double p = rho*(t - cut_time)*exp(-rho*(t - cut_time));
    return p;
}

void approx_BSP::null_emit(double theta, Node *query_node) {
    compute_null_emit_prob(theta, query_node);
    prev_theta = theta;
    prev_node = query_node;
    double ws = 0;
    double *curr_probs = forward_probs[curr_index].data();
    const double *emit_probs = null_emit_probs.data();
    double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    int i = 0;
    for (; i + 3 < dim; i += 4) {
        double p0 = curr_probs[i], p1 = curr_probs[i+1];
        double p2 = curr_probs[i+2], p3 = curr_probs[i+3];
        p0 = (p0 > 0) ? max(epsilon, p0*emit_probs[i]) : 0.0;
        p1 = (p1 > 0) ? max(epsilon, p1*emit_probs[i+1]) : 0.0;
        p2 = (p2 > 0) ? max(epsilon, p2*emit_probs[i+2]) : 0.0;
        p3 = (p3 > 0) ? max(epsilon, p3*emit_probs[i+3]) : 0.0;
        curr_probs[i] = p0; curr_probs[i+1] = p1;
        curr_probs[i+2] = p2; curr_probs[i+3] = p3;
        s0 += p0; s1 += p1; s2 += p2; s3 += p3;
    }
    for (; i < dim; i++) {
        double p = curr_probs[i];
        p = (p > 0) ? max(epsilon, p*emit_probs[i]) : 0.0;
        curr_probs[i] = p;
        s0 += p;
    }
    ws = (s0 + s1) + (s2 + s3);
    assert(ws > 0);
    double inv_ws = 1.0/ws;
    for (int j = 0; j < dim; j++) {
        curr_probs[j] *= inv_ws;
    }
}

void approx_BSP::mut_emit(double theta, double bin_size, vector<double> &mut_set, Node *query_node) {
    compute_mut_emit_probs(theta, bin_size, mut_set, query_node);
    double ws = 0;
    double *curr_probs = forward_probs[curr_index].data();
    const double *emit_probs = mut_emit_probs.data();
    double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    int i = 0;
    for (; i + 3 < dim; i += 4) {
        double p0 = curr_probs[i], p1 = curr_probs[i+1];
        double p2 = curr_probs[i+2], p3 = curr_probs[i+3];
        p0 = (p0 > 0) ? max(epsilon, p0*emit_probs[i]) : 0.0;
        p1 = (p1 > 0) ? max(epsilon, p1*emit_probs[i+1]) : 0.0;
        p2 = (p2 > 0) ? max(epsilon, p2*emit_probs[i+2]) : 0.0;
        p3 = (p3 > 0) ? max(epsilon, p3*emit_probs[i+3]) : 0.0;
        curr_probs[i] = p0; curr_probs[i+1] = p1;
        curr_probs[i+2] = p2; curr_probs[i+3] = p3;
        s0 += p0; s1 += p1; s2 += p2; s3 += p3;
    }
    for (; i < dim; i++) {
        double p = curr_probs[i];
        p = (p > 0) ? max(epsilon, p*emit_probs[i]) : 0.0;
        curr_probs[i] = p;
        s0 += p;
    }
    ws = (s0 + s1) + (s2 + s3);
    assert(ws > 0);
    double inv_ws = 1.0/ws;
    for (int j = 0; j < dim; j++) {
        curr_probs[j] *= inv_ws;
    }
}

map<double, Branch> approx_BSP::sample_joining_branches(int start_index, vector<double> &coordinates) {
    prev_rho = -1;
    map<double, Branch> joining_branches = {};
    int x = curr_index;
    int y = 0;
    double pos = coordinates[x + start_index + 1];
    Interval_ptr interval = sample_curr_interval(x);
    Branch b = interval->branch;
    joining_branches[pos] = b;
    while (x >= 0) {
        vector<Interval_ptr> &intervals = get_state_space(x);
        assert(intervals[sample_index] == interval);
        x = trace_back_helper(interval, x);
        b = interval->branch;
        pos = coordinates[x + start_index];
        joining_branches[pos] = b;
        y = get_prev_breakpoint(x);
        if (x == 0) {
            break;
        } else if (x == y) {
            x -= 1;
            interval = sample_source_interval(interval, x);
            b = interval->branch;
        } else {
            x -= 1;
            interval = sample_prev_interval(x);
            b = interval->branch;
        }
    }
    simplify(joining_branches);
    return joining_branches;
}

Interval_ptr approx_BSP::make_interval(Branch b, double tl, double tu, int init_pos) {
    if (n_arena == arena.size()) {
        arena.emplace_back(b, tl, tu, init_pos);
        n_arena++;
        return &arena.back();
    }
    Interval &x = arena[n_arena++];
    assert(tl <= tu);
    assert(b.lower_node->time <= tl and b.upper_node->time >= tu);
    x.branch = b;
    x.lb = tl;
    x.ub = tu;
    x.start_pos = init_pos;
    x.weight = 0.0;
    x.time = 0.0;
    x.source_pos = 0;
    x.node = nullptr;
    x.reduction = 1.0;
    x.source_weights.clear();
    x.source_intervals.clear();
    return &x;
}

void approx_BSP::set_dimensions() {
    dim = (int) curr_intervals.size();
    time_points.resize(dim); time_points.assign(dim, 0);
    raw_weights.resize(dim); raw_weights.assign(dim, 0);
    recomb_probs.resize(dim); recomb_probs.assign(dim, 0);
    recomb_weights.resize(dim); recomb_weights.assign(dim, 0);
    null_emit_probs.resize(dim); null_emit_probs.assign(dim, 0);
    mut_emit_probs.resize(dim); mut_emit_probs.assign(dim, 0);
}

void approx_BSP::compute_recomb_probs(double rho) {
    if (prev_rho == rho) {
        return;
    }
    double rb = 0;
    for (int i = 0; i < dim; i++) {
        rb = get_recomb_prob(rho, time_points[i]);
        recomb_probs[i] = rb;
    }
}

void approx_BSP::compute_recomb_weights(double rho) {
    if (prev_rho == rho) {
        return;
    }
    for (int i = 0; i < dim; i++) {
        if (curr_intervals[i]->full(cut_time)) {
            recomb_weights[i] = recomb_probs[i]*raw_weights[i];
        }
    }
    weight_sum = accumulate(recomb_weights.begin(), recomb_weights.end(), 0.0);
    for (int i = 0; i < dim; i++) {
        recomb_weights[i] /= weight_sum;
    }
}

void approx_BSP::compute_null_emit_prob(double theta, Node *query_node) {
    if (theta == prev_theta and query_node == prev_node) {
        return;
    }
    for (int i = 0; i < dim; i++) {
        null_emit_probs[i] = eh->null_emit(curr_intervals[i]->branch, time_points[i], theta, query_node);
    }
}

void approx_BSP::compute_mut_emit_probs(double theta, double bin_size, vector<double> &mut_set, Node *query_node) {
    for (int i = 0; i < dim; i++) {
        mut_emit_probs[i] = eh->mut_emit(curr_intervals[i]->branch, time_points[i], theta, bin_size, mut_set, query_node);
    }
}

void approx_BSP::transfer_helper(Interval_info &next_interval, Interval_ptr &prev_interval, double w) {
    recycled_slot(transfer_weights, n_transfer_weights, next_interval).push_back(w);
    recycled_slot(transfer_intervals, n_transfer_intervals, next_interval).push_back(prev_interval);
}

void approx_BSP::transfer_helper(Interval_info &next_interval) {
    recycled_slot(transfer_weights, n_transfer_weights, next_interval);
    recycled_slot(transfer_intervals, n_transfer_intervals, next_interval);
}

void approx_BSP::add_new_branches(Recombination &r) { // add recombined branch and merging branch, if legal
    Interval_info next_interval;
    double lb = 0;
    double ub = 0;
    if (r.merging_branch != Branch() and r.merging_branch.upper_node->time > cut_time) {
        lb = max(cut_time, r.merging_branch.lower_node->time);
        ub = r.merging_branch.upper_node->time;
        next_interval = Interval_info(r.merging_branch, lb, ub);
        transfer_helper(next_interval);
    }
    if (r.recombined_branch != Branch() and r.recombined_branch.upper_node->time > cut_time) {
        lb = max(cut_time, r.recombined_branch.lower_node->time);
        ub = r.recombined_branch.upper_node->time;
        next_interval = Interval_info(r.recombined_branch, lb, ub);
        transfer_helper(next_interval);
    }
}

/*
void approx_BSP::compute_interval_info() {
    double t;
    double p;
    for (int i = 0; i < curr_intervals.size(); i++) {
        Interval_ptr interval = curr_intervals[i];
        p = cc->prob(interval->lb, interval->ub);
        t = cc->find_median(interval->lb, interval->ub);
        interval->weight = p;
        interval->time = t;
        raw_weights[i] = p;
        time_points[i] = t;
    }
    keyed_slot(times, curr_index) = time_points;
    keyed_slot(weights, curr_index) = raw_weights;
}
 */

void approx_BSP::compute_interval_info() {
    double t;
    double p;
    for (int i = 0; i < curr_intervals.size(); i++) {
        const Interval_ptr &interval = curr_intervals[i];
        if (interval->start_pos == curr_index) {
            p = cc->prob(interval->lb, interval->ub);
            t = cc->find_median(interval->lb, interval->ub);
            interval->weight = p;
            interval->time = t;
            raw_weights[i] = p;
            time_points[i] = t;
        } else {
            raw_weights[i] = interval->weight;
            time_points[i] = interval->time;
        }
    }
    keyed_slot(times, curr_index) = time_points;
    keyed_slot(weights, curr_index) = raw_weights;
}

void approx_BSP::sanity_check(Recombination &r) {
    for (int i = 0; i < curr_intervals.size(); i++) {
        const Interval_ptr& interval = curr_intervals[i];
        if (interval->lb == interval->ub and interval->lb == r.inserted_node->time and interval->branch != r.target_branch) {
            forward_probs[curr_index][i] = 0;
        }
        if (interval->lb == interval->ub and interval->lb == r.inserted_node->time and interval->branch == r.target_branch and interval->node != r.inserted_node) {
            forward_probs[curr_index][i] = 0;
        }
    }
}

void approx_BSP::generate_intervals(Recombination &r) {
    Branch b;
    double lb;
    double ub;
    double p;
    vector<Interval_ptr> intervals;
    vector<double> weights;
    Interval_info interval;
    Interval_ptr new_interval = nullptr;
    auto y = transfer_intervals.begin();
    auto x_end = transfer_weights.begin() + n_transfer_weights;
    for (auto x = transfer_weights.begin(); x != x_end; ++x, ++y) {
        interval = x->first;
        auto &weights = x->second;
        auto &intervals = y->second;
        b = interval.branch;
        lb = interval.lb;
        ub = interval.ub;
        p = accumulate(weights.begin(), weights.end(), 0.0);
        assert(!isnan(p));
        if (lb == max(cut_time, b.lower_node->time)) { // full intervals
            new_interval = make_interval(b, lb, ub, curr_index);
            temp_intervals.push_back(new_interval);
            temp.push_back(p);
            if (weights.size() > 0) {
                new_interval->source_weights.assign(weights.begin(), weights.end());
                new_interval->source_intervals.assign(intervals.begin(), intervals.end());
            }
        } else if (p > cutoff) { // partial intervals
            new_interval = make_interval(b, lb, ub, curr_index);
            temp_intervals.push_back(new_interval);
            temp.push_back(p);
            if (weights.size() > 0) {
                new_interval->source_weights.assign(weights.begin(), weights.end());
                new_interval->source_intervals.assign(intervals.begin(), intervals.end());
            }
            if (lb == ub) { // Need to find out where the point mass is from
                if (b == r.merging_branch and lb == r.deleted_node->time) {
                    new_interval->node = r.deleted_node; // creation of a new point mass
                } else {
                    assert(new_interval->source_intervals.size() == 1);
                    new_interval->node = new_interval->source_intervals.front()->node;
                }
            }
            if (new_interval->lb < new_interval->ub) {
                assert(new_interval->node == nullptr);
            } else {
                assert(new_interval->node != nullptr);
            }
        }
    }
    push_row(temp);
    curr_intervals = move(temp_intervals);
}

double approx_BSP::get_overwrite_prob(Recombination &r, double lb, double ub) {
    if (check_points.count(r.pos) > 0) {
        return 0.0;
    }
    double join_time = r.inserted_node->time;
    double p1 = cc->prob(lb, ub);
    double p2 = cc->prob(max(cut_time, r.start_time), join_time);
    if (p1 == 0 and p2 == 0) {
        return 1.0;
    }
    double overwrite_prob = p2/(p1 + p2);
    assert(!isnan(overwrite_prob));
    return overwrite_prob;
}

void approx_BSP::process_interval(Recombination &r, int i) {
    const Branch &prev_branch = curr_intervals[i]->branch;
    if (prev_branch == r.source_branch) {
        process_source_interval(r, i);
    } else if (prev_branch == r.target_branch) {
        process_target_interval(r, i);
    } else {
        process_other_interval(r, i);
    }
}

void approx_BSP::process_source_interval(Recombination &r, int i) {
    double w1, w2, lb, ub = 0;
    Interval_ptr prev_interval = curr_intervals[i];
    double p = forward_probs[curr_index - 1][i];
    double point_time = r.source_branch.upper_node->time;
    double break_time = r.start_time;
    Branch next_branch;
    Interval_info next_interval;
    if (prev_interval->ub <= break_time) {
        lb = prev_interval->lb;
        ub = prev_interval->ub;
        next_branch = r.recombined_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, p);
    } else if (prev_interval->lb >= break_time) {
        lb = point_time;
        ub = point_time;
        next_branch = r.merging_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, p);
    } else {
        w1 = cc->prob(prev_interval->lb, break_time);
        w2 = cc->prob(break_time, prev_interval->ub);
        if (w1 == 0 and w2 == 0) {
            w1 = 1;
            w2 = 0;
        } else {
            w1 = w1/(w1 + w2);
            w2 = 1 - w1;
        }
        lb = prev_interval->lb;
        ub = break_time;
        next_branch = r.recombined_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, w1*p);
        lb = point_time;
        ub = point_time;
        next_branch = r.merging_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, w2*p);
    }
}

void approx_BSP::process_target_interval(Recombination &r, int i) {
    double w0, w1, w2, lb, ub = 0;
    Interval_ptr prev_interval = curr_intervals[i];
    double p = forward_probs[curr_index - 1][i];
    double join_time = r.inserted_node->time;
    Branch next_branch;
    Interval_info next_interval;
    if (prev_interval->lb == prev_interval->ub and prev_interval->lb == join_time) {
        lb = max(cut_time, r.start_time);
        ub = r.recombined_branch.upper_node->time;
        next_branch = r.recombined_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, p);
    } else if (prev_interval->lb >= join_time) {
        lb = prev_interval->lb;
        ub = prev_interval->ub;
        next_branch = r.upper_transfer_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, p);
    } else if (prev_interval->ub <= join_time) {
        lb = prev_interval->lb;
        ub = prev_interval->ub;
        next_branch = r.lower_transfer_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, p);
    } else {
        w0 = get_overwrite_prob(r, prev_interval->lb, prev_interval->ub);
        w1 = cc->prob(prev_interval->lb, join_time);
        w2 = cc->prob(join_time, prev_interval->ub);
        if (w1 + w2 == 0) {
            w1 = 0;
            w2 = 0;
            w0 = 1.0;
        } else {
            w1 = w1/(w1 + w2);
            w2 = 1 - w1;
            w1 *= 1 - w0;
            w2 *= 1 - w0;
        }
        lb = prev_interval->lb;
        ub = join_time;
        next_branch = r.lower_transfer_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, w1*p);
        lb = join_time;
        ub = prev_interval->ub;
        next_branch = r.upper_transfer_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, w2*p);
        lb = max(r.start_time, cut_time);
        ub = r.recombined_branch.upper_node->time;;
        next_branch = r.recombined_branch;
        next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, w0*p);
    }
}

void approx_BSP::process_other_interval(Recombination &r, int i) {
    double lb, ub = 0;
    Interval_ptr prev_interval = curr_intervals[i];
    double p = forward_probs[curr_index - 1][i];
    if (prev_interval->branch != r.source_sister_branch and prev_interval->branch != r.source_parent_branch) {
        // in other words, not affected by recombination
        if (prev_interval->full(cut_time)) {
            temp_intervals.push_back(prev_interval);
            temp.push_back(p);
        } else if (p > cutoff) {
            temp_intervals.push_back(prev_interval);
            temp.push_back(p);
        }
    } else if (p > cutoff) { // will not create a full branch, so we need to prune
        lb = prev_interval->lb;
        ub = prev_interval->ub;
        Branch &next_branch = r.merging_branch;
        Interval_info next_interval = Interval_info(next_branch, lb, ub);
        transfer_helper(next_interval, prev_interval, p);
    }
}

double approx_BSP::random() {
    double p = uniform_random();
    return p;
}

int approx_BSP::get_prev_breakpoint(int x) {
    auto state_it = upper_bound(state_spaces.begin(), state_spaces.end(), x,
                                [](int k, const pair<int, vector<Interval_ptr>> &e) { return k < e.first; });
    state_it--;
    return state_it->first;
}

vector<Interval_ptr> &approx_BSP::get_state_space(int x) {
    return keyed_lookup(state_spaces, x);
}

vector<double> &approx_BSP::get_time_points(int x) {
    return keyed_lookup(times, x);
}

vector<double> &approx_BSP::get_raw_weights(int x) {
    return keyed_lookup(weights, x);
}

int approx_BSP::get_interval_index(Interval_ptr interval, vector<Interval_ptr > &intervals) {
    auto it = find(intervals.begin(), intervals.end(), interval);
    int index = (int) distance(intervals.begin(), it);
    return index;
}
 
void approx_BSP::simplify(map<double, Branch> &joining_branches) {
    map<double, Branch> simplified_joining_branches = {};
    Branch curr_branch = joining_branches.begin()->second;
    simplified_joining_branches[joining_branches.begin()->first] = curr_branch;
    for (auto x : joining_branches) {
        if (x.second != curr_branch) {
            simplified_joining_branches.insert(x);
            curr_branch = x.second;
        }
    }
    simplified_joining_branches[joining_branches.rbegin()->first] = joining_branches.rbegin()->second;
    joining_branches = simplified_joining_branches;
}

Interval_ptr approx_BSP::sample_curr_interval(int x) {
    vector<Interval_ptr > &intervals = get_state_space(x);
    double ws = accumulate(forward_probs[x].begin(), forward_probs[x].end(), 0.0);
    double q = random();
    double w = ws*q;
    for (int i = 0; i < intervals.size(); i++) {
        w -= forward_probs[x][i];
        if (w <= 0) {
            sample_index = i;
            return intervals[i];
        }
    }
    cerr << "approx_BSP sample_curr_interval failed" << endl;
    exit(1);
}

Interval_ptr approx_BSP::sample_prev_interval(int x) {
    vector<Interval_ptr > &intervals = get_state_space(x);
    vector<double> &prev_times = get_time_points(x);
    double rho = rhos[x];
    double ws = recomb_sums[x];
    double q = random();
    double w = ws*q;
    double rb = 0;
    for (int i = 0; i < intervals.size(); i++) {
        rb = get_recomb_prob(rho, prev_times[i]);
        w -= rb*forward_probs[x][i];
        if (w <= 0) {
            sample_index = i;
            return intervals[i];
        }
    }
    cerr << "approx_BSP sample_prev_interval failed" << endl;
    exit(1);
}

Interval_ptr approx_BSP::sample_source_interval(Interval_ptr interval, int x) {
    vector<Interval_ptr> &intervals = interval->source_intervals;
    vector<double> &weights = interval->source_weights;
    vector<Interval_ptr> &prev_intervals = get_state_space(x);
    if (x == interval->start_pos - 1) {
        double q = random();
        double ws = accumulate(weights.begin(), weights.end(), 0.0);
        double w = ws*q;
        for (int i = 0; i < weights.size(); i++) {
            w -= weights[i];
            if (w <= 0) {
                sample_index = get_interval_index(intervals[i], prev_intervals);
                return intervals[i];
            }
        }
        cerr << "approx bsp sample_source_interval failed" << endl;
        exit(1);
    } else {
        sample_index = get_interval_index(interval, prev_intervals);
        assert(prev_intervals[sample_index] == interval);
        return interval;
    }
}

int approx_BSP::trace_back_helper(Interval_ptr interval, int x) {
    int y = get_prev_breakpoint(x);
    if (!interval->full(cut_time)) {
        return y;
    }
    vector<double> &ts = get_time_points(x);
    vector<double> &ws = get_raw_weights(x);
    assert(sample_index < ts.size());
    double t = ts[sample_index];
    double w = ws[sample_index];
    double p = random();
    double q = 1;
    double shrinkage = 0;
    double recomb_prob = 0;
    double non_recomb_prob = 0;
    double all_prob = 0;
    while (x > y) {
        recomb_sum = recomb_sums[x - 1];
        weight_sum = weight_sums[x];
        if (recomb_sum == 0) {
            shrinkage = 1;
        } else {
            recomb_prob = get_recomb_prob(rhos[x - 1], t);
            non_recomb_prob = (1 - recomb_prob)*forward_probs[x - 1][sample_index];
            all_prob = non_recomb_prob + recomb_sum*w*recomb_prob/weight_sum;
            shrinkage = non_recomb_prob/all_prob;
            assert(shrinkage >= 0 and shrinkage <= 1);
        }
        q *= shrinkage;
        if (p >= q) {
            return x;
        }
        x -= 1;
    }
    return y;
}

double approx_BSP::avg_num_states() {
    int span = 0;
    double count = 0;
    auto x = state_spaces.begin();
    ++x;
    while (x->first != INT_MAX) {
        count += x->second.size()*(x->first - (x - 1)->first);
        span = x->first;
        ++x;
    }
    double avg = (double) count/span;
    return avg;
}

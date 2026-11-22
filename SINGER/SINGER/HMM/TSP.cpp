//
//  TSP.cpp
//  SINGER
//
//  Created by Yun Deng on 9/25/23.
//

#include "TSP.hpp"
#include <cstring>

int TSP::counter = 0;

TSP::TSP() {
}

TSP::~TSP() {
    vector<vector<double>>().swap(forward_probs);
    vector<pair<int, vector<Interval *>>>().swap(state_spaces);
}

void TSP::set_gap(double q) {
    gap = q;
}

void TSP::set_emission(shared_ptr<Emission> e) {
    eh = e;
}

void TSP::set_check_points(set<double> &p) {
    check_points = p;
}

Interval *TSP::make_interval(Branch b, double tl, double tu, int init_pos) {
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
    x.source_interval = nullptr;
    x.source_weights.clear();
    x.source_intervals.clear();
    return &x;
}

void TSP::push_row(const vector<double> &v) {
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

void TSP::reset() {
    state_spaces[0].first = INT_MAX;
    state_spaces[0].second.clear();
    n_spaces = 1;
    n_rows = 0;
    curr_index = 0;
    curr_branch = Branch();
    curr_intervals.clear();
    swaps.clear();
    applied = 0;
    n_kernels = 0;
    kernel_at.clear();
    lower_sums.clear();
    upper_sums.clear();
    diagonals.clear();
    factors.clear();
    lower_diagonals.clear();
    temp.clear();
    null_emit_probs.clear();
    mut_emit_probs.clear();
    trace_back_probs.clear();
    check_points.clear();
    prev_rho = -1;
    prev_theta = -1;
    prev_node = nullptr;
    null_logs_valid = false;
    tb_kernel = -1;
    dim = 0;
    sample_index = -1;
    n_arena = 0;
}

void TSP::save_state(TSP &s) {
    deque<Interval> a;
    vector<vector<double>> f;
    vector<pair<int, vector<Interval *>>> p;
    vector<Kernel> k;
    Store st;
    auto park = [&]() {
        arena.swap(a);
        forward_probs.swap(f);
        state_spaces.swap(p);
        kernels.swap(k);
        swap(store, st);
    };
    park();
    s = *this;
    park();
    s.cc = make_shared<coalescent_calculator>(*cc);
    s.cc1 = make_shared<coalescent_calculator>(*cc1);
    s.saved_row = forward_probs[curr_index];
    s.saved_intervals.clear();
    for (Interval *iv : curr_intervals) {
        s.saved_intervals.push_back(*iv);
    }
}

void TSP::restore_state(TSP &s) {
    deque<Interval> a;
    vector<vector<double>> f;
    vector<pair<int, vector<Interval *>>> p;
    vector<Kernel> k;
    Store st;
    auto park = [&]() {
        arena.swap(a);
        forward_probs.swap(f);
        state_spaces.swap(p);
        kernels.swap(k);
        swap(store, st);
    };
    park();
    *this = s;
    park();
    cc = make_shared<coalescent_calculator>(*s.cc);
    cc1 = make_shared<coalescent_calculator>(*s.cc1);
    forward_probs[curr_index] = saved_row;
    for (size_t i = 0; i < curr_intervals.size(); i++) {
        *curr_intervals[i] = saved_intervals[i];
    }
    state_spaces[n_spaces - 1].first = INT_MAX;
    state_spaces[n_spaces - 1].second.clear();
}

double TSP::per_length(double rho, double L) {
    return -expm1(-rho*L)/L;
}

void TSP::reserve_memory(int length) {
    forward_probs.reserve(length);
}

void TSP::start(Branch &branch, double t) {
    cut_time = t;
    cd_valid = false;
    null_logs_valid = false;
    cc_lowest_change = numeric_limits<double>::infinity();
    curr_index = 0;
    curr_branch = branch;
    generate_intervals(branch, branch.lower_node->time, branch.upper_node->time);
    set_dimensions();
    compute_factors();
    for (int i = 0; i < curr_intervals.size(); i++) {
        temp[i] = cc->prob(curr_intervals[i]->lb, curr_intervals[i]->ub);
    }
    set_state_space(0);
    push_row(temp);
    temp.clear();
}

void TSP::transfer(Recombination &r, Branch &prev_branch, Branch &next_branch) {
    prev_rho = -1;
    prev_theta = -1;
    prev_node = nullptr;
    null_logs_valid = false;
    sanity_check(r);
    size_t num_times = cc->num_times();
    cc->update(r);
    cc1->update(r);
    if (cc->num_times() != num_times) {
        cc_lowest_change = -numeric_limits<double>::infinity();
    } else {
        if (r.deleted_node->time > cut_time) {
            cc_lowest_change = min(cc_lowest_change, r.deleted_node->time);
        }
        if (r.inserted_node->time > cut_time) {
            cc_lowest_change = min(cc_lowest_change, r.inserted_node->time);
        }
    }
    curr_index += 1;
    swaps.push_back({curr_index, {r.deleted_node->time, r.inserted_node->time}});
    applied = swaps.size();
    curr_branch = next_branch;
    if (prev_branch == r.source_branch and next_branch == r.merging_branch) {
        set_interval_constraint(r);
    } else if (prev_branch == r.target_branch and next_branch == r.recombined_branch) {
        set_point_constraint(r);
    }
    curr_intervals.clear();
    if (prev_branch == r.source_branch and next_branch == r.merging_branch) {
        double t = r.deleted_node->time;
        generate_intervals(next_branch, next_branch.lower_node->time, t);
        generate_intervals(next_branch, t, t);
        temp.back() = 1.0f;
        curr_intervals.back()->node = r.deleted_node;
        generate_intervals(next_branch, t, next_branch.upper_node->time);
    } else if (prev_branch == r.target_branch and next_branch == r.recombined_branch) {
        generate_intervals(next_branch, next_branch.lower_node->time, r.start_time);
        generate_intervals(next_branch, r.start_time, next_branch.upper_node->time);
        for (int i = 0; i < curr_intervals.size(); i++) {
            if (curr_intervals[i]->lb >= r.start_time) {
                temp[i] = 1.0;
            }
        }
    } else {
        generate_intervals(next_branch, next_branch.lower_node->time, max(prev_branch.lower_node->time, next_branch.lower_node->time));
        transfer_intervals(r, prev_branch, next_branch);
        generate_intervals(next_branch, min(curr_intervals.back()->ub, next_branch.upper_node->time), next_branch.upper_node->time);
    }
    set_state_space(curr_index);
    push_row(temp);
    temp.clear();
    set_dimensions();
    compute_factors();
}

void TSP::recombine(Branch &prev_branch, Branch &next_branch) {
    assert(next_branch != Branch());
    vector<Interval *> prev_intervals = curr_intervals;
    curr_intervals.clear();
    prev_rho = -1;
    prev_theta = -1;
    prev_node = nullptr;
    null_logs_valid = false;
    curr_branch = next_branch;
    curr_index += 1;
    generate_intervals(next_branch, next_branch.lower_node->time, next_branch.upper_node->time);
    push_row(temp);
    set_state_space(curr_index);
    set_dimensions();
    compute_factors();
    double mass = accumulate(forward_probs[curr_index-1].begin(),
                             forward_probs[curr_index-1].begin() + prev_intervals.size(), 0.0);
    vector<double> w(curr_intervals.size(), 0.0);
    double tot = 0;
    for (int j = 0; j < curr_intervals.size(); j++) {
        double lb = max(curr_intervals[j]->lb, cut_time);
        double ub = curr_intervals[j]->ub;
        if (ub <= lb) continue;
        for (int i = 0; i < prev_intervals.size(); i++) {
            double p = forward_probs[curr_index-1][i];
            if (p > 0) w[j] += p*own_mass(prev_intervals[i]->time, lb, ub);
        }
        tot += w[j];
    }
    for (int j = 0; j < curr_intervals.size(); j++) {
        forward_probs[curr_index][j] += mass*w[j]/tot;
    }
    for (int i = 0; i < forward_probs[curr_index].size(); i++) {
        assert(forward_probs[curr_index][i] >= 0);
    }
    temp.clear();
}

vector<double> TSP::generate_grid(double lb, double ub) {
    assert(lb < ub);
    vector<double> points = {lb};
    double ls = cc->surv(lb);
    double us = cc->surv(ub);
    double q = ls - us;
    int n = max((int) ceil((1 - us/ls)/gap), min_num);
    for (int i = 1; i < n; i++) {
        points.emplace_back(cc->surv_inv(ls - i*q/n));
    }
    points.emplace_back(ub);
    return points;
}

void TSP::forward(double rho) {
    compute_diagonals(rho);
    if (rho != prev_rho) {
        if (n_kernels == (int) kernels.size()) {
            kernels.emplace_back();
        }
        Kernel &k = kernels[n_kernels++];
        k.lower = lower_diagonals;
        k.factors = factors;
        k.stays = stays;
        k.selfs = selfs;
        k.masses = masses;
    }
    compute_lower_sums();
    compute_upper_sums();
    curr_index += 1;
    prev_rho = rho;
    if ((int) kernel_at.size() <= curr_index) {
        kernel_at.resize(curr_index + 1, -1);
    }
    kernel_at[curr_index] = n_kernels - 1;
    push_row(lower_sums);
    for (int i = 0; i < dim; i++) {
        assert(forward_probs[curr_index][i] >= 0);
        forward_probs[curr_index][i] += diagonals[i]*forward_probs[curr_index-1][i] + lower_diagonals[i]*upper_sums[i];
        if (rho > 0 and (curr_intervals[i]->lb != curr_intervals[i]->ub or forward_probs[curr_index][i] > 0)) {
            forward_probs[curr_index][i] = max(epsilon, forward_probs[curr_index][i]);
        }
    }
}

void TSP::null_emit(double theta, Node *query_node) {
    compute_null_emit_probs(theta, query_node);
    prev_theta = theta;
    prev_node = query_node;
    double ws = 0;
    assert(dim == forward_probs[curr_index].size());
    double *curr_probs = forward_probs[curr_index].data();
    const double *emit_probs = null_emit_probs.data();
    double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    int i = 0;
    for (; i + 3 < dim; i += 4) {
        assert(curr_probs[i] >= 0 and curr_probs[i+1] >= 0);
        assert(curr_probs[i+2] >= 0 and curr_probs[i+3] >= 0);
        curr_probs[i] *= emit_probs[i];     curr_probs[i+1] *= emit_probs[i+1];
        curr_probs[i+2] *= emit_probs[i+2]; curr_probs[i+3] *= emit_probs[i+3];
        s0 += curr_probs[i];     s1 += curr_probs[i+1];
        s2 += curr_probs[i+2];   s3 += curr_probs[i+3];
    }
    for (; i < dim; i++) {
        assert(curr_probs[i] >= 0);
        curr_probs[i] *= emit_probs[i];
        s0 += curr_probs[i];
    }
    ws = (s0 + s1) + (s2 + s3);
    if (!(ws > 0)) {
        cerr << "Error: the time HMM's forward probabilities sum to " << ws << " at step " << curr_index << "." << endl;
        exit(1);
    }
    double inv_ws = 1.0/ws;
    for (int j = 0; j < dim; j++) {
        curr_probs[j] *= inv_ws;
    }
}

void TSP::mut_emit(Node *query_node) {
    compute_mut_emit_probs(query_node);
    double ws = 0;
    double *curr_probs = forward_probs[curr_index].data();
    const double *emit_probs = mut_emit_probs.data();
    double s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    int i = 0;
    for (; i + 3 < dim; i += 4) {
        curr_probs[i] *= emit_probs[i];     curr_probs[i+1] *= emit_probs[i+1];
        curr_probs[i+2] *= emit_probs[i+2]; curr_probs[i+3] *= emit_probs[i+3];
        s0 += curr_probs[i];     s1 += curr_probs[i+1];
        s2 += curr_probs[i+2];   s3 += curr_probs[i+3];
    }
    for (; i < dim; i++) {
        curr_probs[i] *= emit_probs[i];
        s0 += curr_probs[i];
    }
    ws = (s0 + s1) + (s2 + s3);
    if (!(ws > 0)) {
        cerr << "Error: the time HMM's forward probabilities sum to " << ws << " at step " << curr_index << "." << endl;
        exit(1);
    }
    double inv_ws = 1.0/ws;
    for (int j = 0; j < dim; j++) {
        curr_probs[j] *= inv_ws;
    }
}

void TSP::rewind_to(int x) {
    while (applied > 0 and swaps[applied - 1].first > x) {
        cc->update(swaps[applied - 1].second.second, swaps[applied - 1].second.first);
        cc1->update(swaps[applied - 1].second.second, swaps[applied - 1].second.first);
        applied -= 1;
    }
}

void TSP::forward_to(int x) {
    while (applied < swaps.size() and swaps[applied].first <= x) {
        cc->update(swaps[applied].second.first, swaps[applied].second.second);
        cc1->update(swaps[applied].second.first, swaps[applied].second.second);
        applied += 1;
    }
}

map<double, Node *> TSP::sample_joining_nodes(int start_index, vector<double> &coordinates) {
    forward_to(curr_index);
    sel_log_q = 0;
    time_log_q = 0;
    map<double, Node *> joining_nodes = {};
    int x = curr_index;
    pin_active = false;
    Interval *interval = sample_curr_interval(x);
    Node *n = sample_joining_node(interval);
    joining_nodes[coordinates[x + start_index + 1]] = nullptr;
    while (true) {
        x = trace_back_helper(interval, x);
        joining_nodes[coordinates[x + start_index]] = n;
        if (x == 0) {
            break;
        } if (x == interval->start_pos) {
            if (interval->source_interval != nullptr) {
                x -= 1;
                interval = sample_source_interval(interval, x);
                rewind_to(x);
            } else {
                x -= 1;
                interval = sample_recomb_interval(interval, x);
                rewind_to(x);
                pin_active = interval->node != nullptr;
                n = sample_joining_node(interval);
            }
        } else {
            x -= 1;
            interval = sample_prev_interval(interval, x);
            rewind_to(x);
            pin_active = false;
            n = sample_joining_node(interval);
        }
    }
    return joining_nodes;
}

double TSP::non_recomb_prob(double rho, double s) {
    return exp(-rho*(s - cut_time));
}

void TSP::fill_interval_time(Interval *iv) {
    iv->fill_time();
}

double TSP::own_mass(double t, double lb, double ub) {
    lb = max(lb, cut_time);
    if (ub <= lb) {
        return 0;
    }
    double m = 0;
    double tb = min(ub, t);
    if (tb > lb) {
        m += cc1->recomb_mass(t, tb) - cc1->recomb_mass(t, lb);
    }
    if (ub > t) {
        double den = cc1->prob(t, ub);
        if (den > 0) {
            m += (cc1->recomb_mass(t, ub) - cc1->recomb_mass(t, t))/den*cc1->surv(t)*cc->prob(max(lb, t), ub)/cc->surv(t);
        }
    }
    return m;
}

void TSP::generate_intervals(Branch &next_branch, double lb, double ub) {
    Interval *new_interval = nullptr;
    lb = max(cut_time, lb);
    ub = max(cut_time, ub);
    if (lb == ub) {
        if (lb == max(cut_time, next_branch.lower_node->time) or lb == next_branch.upper_node->time) {
            return;
        }
        else {
            new_interval = make_interval(next_branch, lb, ub, curr_index);
            fill_interval_time(new_interval);
            curr_intervals.emplace_back(new_interval);
            temp.emplace_back(0);
            return;
        }
    }
    vector<double> points = generate_grid(lb, ub);
    for (int i = 0; i < points.size() - 1; i++) {
        new_interval = make_interval(next_branch, points[i], points[i+1], curr_index);
        fill_interval_time(new_interval);
        curr_intervals.emplace_back(new_interval);
        temp.emplace_back(0);
    }
}

void TSP::transfer_intervals(Recombination &r, Branch &prev_branch, Branch &next_branch) {
    double lb;
    double ub;
    double p;
    vector<Interval *> &prev_intervals = get_state_space(curr_index - 1);
    Interval *interval = nullptr;
    Interval *new_interval = nullptr;
    for (int i = 0; i < prev_intervals.size(); i++) {
        interval = prev_intervals[i];
        lb = max(interval->lb, next_branch.lower_node->time);
        ub = min(interval->ub, next_branch.upper_node->time);
        if (prev_branch == r.source_branch) {
            ub = min(ub, r.start_time);
            if (lb == r.start_time) {
                continue;
            }
        }
        if (lb == ub and ub == next_branch.upper_node->time) {
            continue;
        }
        if (lb == ub and ub == next_branch.lower_node->time) {
            continue;
        }
        if (ub >= lb) {
            double w = get_prop(lb, ub, interval->lb, interval->ub);
            if (w > 0 and forward_probs[curr_index - 1][i] > 0) {
                p = max(epsilon, w*forward_probs[curr_index-1][i]);
            } else {
                p = 0;
            }
            new_interval = make_interval(next_branch, lb, ub, curr_index);
            fill_interval_time(new_interval);
            new_interval->node = interval->node;
            new_interval->source_interval = interval;
            curr_intervals.emplace_back(new_interval);
            temp.emplace_back(p);
        }
    }
}

void TSP::set_dimensions() {
    dim = (int) curr_intervals.size();
    diagonals.resize(dim); diagonals.assign(dim, 0);
    lower_diagonals.resize(dim); lower_diagonals.assign(dim, 0);
    lower_sums.resize(dim); lower_sums.assign(dim, 0);
    upper_sums.resize(dim); upper_sums.assign(dim, 0);
    null_emit_probs.resize(dim); null_emit_probs.assign(dim, 0);
    mut_emit_probs.resize(dim); mut_emit_probs.assign(dim, 0);
    factors.resize(dim); factors.assign(dim, 0);
    masses.resize(dim); masses.assign(dim, 0);
    stays.assign(dim, 0);
    selfs.assign(dim, 0);
}

double TSP::get_prop(double lb1, double ub1, double lb2, double ub2) {
    if (ub2 - lb2 < 1e-6) {
        return 1;
    }
    return cc->prob(lb1, ub1)/cc->prob(lb2, ub2);
}

void TSP::compute_null_emit_probs(double theta, Node *query_node) {
    if (!eh->model->branch_product) {
        update_norms(query_node);
        copy(norms.begin(), norms.end(), null_emit_probs.begin());
        return;
    }
    if (theta == prev_theta and query_node == prev_node and eh->bin_width == prev_width and eh->crho == prev_crho and eh->norm_scale == prev_scale) {
        return;
    }
    update_norms(query_node);
    prev_width = eh->bin_width;
    prev_crho = eh->crho;
    prev_scale = eh->norm_scale;
    copy(norms.begin(), norms.end(), null_emit_probs.begin());
}

void TSP::update_norms(Node *query_node) {
    if (!eh->model->branch_product) {
        int S = eh->model->states;
        if (!(null_logs_valid and query_node == log_node and eh->unit == log_unit)) {
            int K = (int) eh->model->rates.size();
            null_logs.resize((S + 1)*dim);
            edge_cache.resize(K*dim);
            for (int i = 0; i < dim; i++) {
                eh->edge_matrices(curr_branch, curr_intervals[i]->time, query_node, &edge_cache[K*i]);
                eh->null_logs_base(curr_branch, &edge_cache[K*i], &null_logs[(S + 1)*i]);
            }
            null_logs_valid = true;
            log_node = query_node;
            log_unit = eh->unit;
        }
        norms.resize(dim);
        const double *c = eh->base_counts;
        double u = eh->unknown_count;
        for (int i = 0; i < dim; i++) {
            const double *l = &null_logs[(S + 1)*i];
            double v = 0;
            for (int b = 0; b < S; b++) {
                v += c[b]*l[b];
            }
            norms[i] = exp(-(v + u*l[S]));
        }
        return;
    }
    if (!(null_logs_valid and query_node == log_node and eh->crho == log_crho and eh->norm_scale == log_scale)) {
        null_logs.resize(dim);
        for (int i = 0; i < dim; i++) {
            null_logs[i] = eh->null_log(curr_branch, curr_intervals[i]->time, query_node);
        }
        null_logs_valid = true;
        log_node = query_node;
        log_crho = eh->crho;
        log_scale = eh->norm_scale;
        logs_stamp++;
    }
    if (norm_stamp == logs_stamp and norm_width == eh->bin_width) {
        return;
    }
    norms.swap(alt_norms);
    swap(norm_stamp, alt_stamp);
    swap(norm_width, alt_width);
    if (norm_stamp == logs_stamp and norm_width == eh->bin_width) {
        return;
    }
    norms.resize(dim);
    for (int i = 0; i < dim; i++) {
        norms[i] = exp(-eh->bin_width*null_logs[i]);
    }
    norm_stamp = logs_stamp;
    norm_width = eh->bin_width;
}

void TSP::compute_mut_emit_probs(Node *query_node) {
    update_norms(query_node);
    if (!eh->model->branch_product) {
        int K = (int) eh->model->rates.size();
        for (int i = 0; i < dim; i++) {
            mut_emit_probs[i] = eh->mut_emit(curr_branch, norms[i], &edge_cache[K*i]);
        }
        return;
    }
    for (int i = 0; i < dim; i++) {
        mut_emit_probs[i] = eh->mut_emit(curr_branch, curr_intervals[i]->time, query_node, norms[i]);
    }
}

void TSP::compute_diagonals(double rho) {
    if (rho == prev_rho) {
        return;
    }
    if (cd_valid and rho == cd_rho and sister_mass == cd_sister_mass and cc_lowest_change > cd_top and (int) cd_lb.size() == dim) {
        bool same = true;
        for (int i = 0; i < dim and same; i++) {
            Interval *iv = curr_intervals[i];
            same = iv->lb == cd_lb[i] and iv->ub == cd_ub[i] and iv->time == cd_time[i];
        }
        if (same) {
            diagonals = cd_diagonals;
            lower_diagonals = cd_lower;
            factors = cd_factors;
            stays = cd_stays;
            selfs = cd_selfs;
            return;
        }
    }
    if (!find_stored(rho)) {
        assemble(*find_masses(), rho);
        keep_stored(rho);
    }
    cd_lb.resize(dim); cd_ub.resize(dim); cd_time.resize(dim);
    for (int i = 0; i < dim; i++) {
        cd_lb[i] = curr_intervals[i]->lb;
        cd_ub[i] = curr_intervals[i]->ub;
        cd_time[i] = curr_intervals[i]->time;
    }
    cd_diagonals = diagonals;
    cd_lower = lower_diagonals;
    cd_factors = factors;
    cd_stays = stays;
    cd_selfs = selfs;
    cd_rho = rho;
    cd_sister_mass = sister_mass;
    cd_top = curr_intervals[dim-1]->ub;
    cc_lowest_change = numeric_limits<double>::infinity();
    cd_valid = true;
}

static bool same_bits(double x, double y) {
    return memcmp(&x, &y, sizeof(double)) == 0;
}

static bool same_bits(const vector<double> &x, const vector<double> &y) {
    return x.size() == y.size() and memcmp(x.data(), y.data(), x.size()*sizeof(double)) == 0;
}

bool TSP::find_stored(double rho) {
    key.clear();
    for (Interval *iv : curr_intervals) {
        key.push_back(iv->lb);
        key.push_back(iv->ub);
        key.push_back(iv->time);
    }
    for (size_t k = 0; k < store.n; k++) {
        Stored &e = store.entries[k];
        if (e.branch == curr_branch and same_bits(e.rho, rho) and same_bits(e.sister_mass, sister_mass)
            and same_bits(e.cut_time, cut_time) and same_bits(e.key, key)
            and same_bits(e.times, cc->times()) and e.signs == cc->signs()) {
            diagonals = e.diagonals;
            lower_diagonals = e.lower;
            factors = e.factors;
            stays = e.stays;
            selfs = e.selfs;
            return true;
        }
    }
    return false;
}

void TSP::keep_stored(double rho) {
    if (store.n == store.entries.size()) {
        store.entries.emplace_back();
    }
    Stored &e = store.entries[store.n++];
    e.branch = curr_branch;
    e.rho = rho;
    e.sister_mass = sister_mass;
    e.cut_time = cut_time;
    e.key = key;
    e.times = cc->times();
    e.signs = cc->signs();
    e.diagonals = diagonals;
    e.lower = lower_diagonals;
    e.factors = factors;
    e.stays = stays;
    e.selfs = selfs;
}

TSP::Masses *TSP::find_masses() {
    for (size_t k = 0; k < mass_store.n; k++) {
        Masses &e = mass_store.entries[k];
        if (e.branch == curr_branch and same_bits(e.sister_mass, sister_mass) and same_bits(e.cut_time, cut_time)
            and same_bits(e.key, key) and same_bits(e.times, cc->times()) and e.signs == cc->signs()) {
            return &e;
        }
    }
    if (mass_store.n == mass_store.entries.size()) {
        mass_store.entries.emplace_back();
    }
    Masses &e = mass_store.entries[mass_store.n++];
    e.branch = curr_branch;
    e.sister_mass = sister_mass;
    e.cut_time = cut_time;
    e.key = key;
    e.times = cc->times();
    e.signs = cc->signs();
    build_masses(e);
    return &e;
}

void TSP::assemble(const Masses &e, double rho) {
    double lo = max(curr_branch.lower_node->time, cut_time);
    for (int i = 0; i < dim; i++) {
        double t = e.time[i];
        double stay;
        stay = non_recomb_prob(rho, t);
        bool small = t - cut_time < 0.005;
        double pref = small ? per_length(rho, t - cut_time)*(t - cut_time) : per_length(rho, t - cut_time);
        double pref2 = small ? per_length(rho, t - lo)*(t - lo) : per_length(rho, t - lo);
        auto jump = [&](const Jump &j) {
            if (rho == 0) return 0.0;
            if (j.empty) return epsilon;
            double p;
            if (small) {
                p = pref*j.d;
                p += lo == cut_time ? p : pref2*j.d;
            } else {
                p = pref*j.m;
                if (lo == cut_time) {
                    p += p;
                } else {
                    p += pref2*j.m2;
                }
            }
            if (sister_mass > 0) {
                p += rho*sister_mass*j.q;
            }
            return max(p, epsilon);
        };
        double base = stay + jump(e.base[i]);
        double self = jump(e.self[i]);
        diagonals[i] = (stay + self)/base;
        stays[i] = stay/base;
        selfs[i] = self/base;
        if (i > 0) {
            lower_diagonals[i-1] = jump(e.lower[i])/base;
        }
        if (i + 1 < dim) {
            double m = e.fm[i];
            factors[i] = (m > 0) ? jump(e.factor[i])/base/m : 0;
        }
        assert(!isnan(diagonals[i]));
    }
    lower_diagonals[dim-1] = 0;
    factors[dim-1] = 0;
}

void TSP::build_masses(Masses &e) {
    e.time.assign(dim, 0.0);
    e.fm.assign(dim, 0.0);
    e.base.assign(dim, Jump());
    e.self.assign(dim, Jump());
    e.lower.assign(dim, Jump());
    e.factor.assign(dim, Jump());
    for (vector<double> *v : {&A, &Bv, &BA, &BB, &GA, &GB, &HA, &HB, &SA, &SB}) {
        v->assign(dim, 0.0);
    }
    double R;
    int p0 = cc->locate(max(curr_intervals[0]->lb, cut_time));
    int p = p0;
    for (int j = 0; j < dim; j++) {
        A[j] = max(curr_intervals[j]->lb, cut_time);
        Bv[j] = curr_intervals[j]->ub;
        if (j > 0 and A[j] == Bv[j-1]) {
            BA[j] = BB[j-1];
            GA[j] = GB[j-1];
            HA[j] = HB[j-1];
            SA[j] = SB[j-1];
        } else {
            p = cc->advance(A[j], p);
            BA[j] = cc1->recomb_mass(A[j], p, 0, R);
            GA[j] = cc->tail(A[j], p, SA[j]);
            HA[j] = cc1->tail(A[j], p, R);
        }
        if (Bv[j] == A[j]) {
            BB[j] = BA[j];
            GB[j] = GA[j];
            HB[j] = HA[j];
            SB[j] = SA[j];
        } else {
            p = cc->advance(Bv[j], p);
            BB[j] = cc1->recomb_mass(Bv[j], p, 0, R);
            GB[j] = cc->tail(Bv[j], p, SB[j]);
            HB[j] = cc1->tail(Bv[j], p, R);
        }
    }
    double top = Bv[dim-1];
    double lo = max(curr_branch.lower_node->time, cut_time);
    double c = 0;
    if (lo > cut_time) {
        double b1, r1;
        b1 = cc1->recomb_mass(lo, cc->locate(lo), 1.0, r1);
        c = r1 - b1;
    }
    p = p0;
    for (int i = 0; i < dim; i++) {
        double t = curr_intervals[i]->time;
        e.time[i] = t;
        p = cc->advance(t, p);
        double S, S1, Rt;
        double Gt = cc->tail(t, p, S);
        double Ht = cc1->tail(t, p, S1);
        double d = top > t ? max(Ht - HB[dim-1], 0.0) : 0;
        double Bt = cc1->recomb_mass(t, p, d, Rt);
        double W = d > 0 ? (Rt - Bt)/d*S1 : 0;
        bool small = t - cut_time < 0.005;
        double W2 = W - S1*c;
        double sis = exp(-(t - cut_time));
        auto mass = [&](double a, double b, double Ba, double Bb, double Ga, double Gb, double Ha, double Hb, double Sa, double Sb) {
            Jump j;
            j.empty = b <= a;
            if (j.empty) return j;
            double tb = min(b, t);
            if (small) {
                j.d = Sa - Sb;
            } else {
                double m = 0;
                if (tb > a) m += (tb == b ? Bb : Bt) - Ba;
                if (b > t) m += W*((a >= t ? Ga : Gt) - Gb)/S;
                j.m = m;
                if (lo != cut_time) {
                    double m2 = 0;
                    if (tb > a) m2 += (tb == b ? Bb : Bt) - Ba - c*(Ha - (tb == b ? Hb : Ht));
                    if (b > t) m2 += W2*((a >= t ? Ga : Gt) - Gb)/S;
                    j.m2 = m2;
                }
            }
            if (sister_mass > 0) {
                double q = 0;
                if (tb > a) q += Ha - (tb == b ? Hb : Ht);
                if (b > t) q += sis*((a >= t ? Ga : Gt) - Gb);
                j.q = q;
            }
            return j;
        };
        e.base[i] = mass(A[0], top, BA[0], BB[dim-1], GA[0], GB[dim-1], HA[0], HB[dim-1], SA[0], SB[dim-1]);
        e.self[i] = mass(A[i], Bv[i], BA[i], BB[i], GA[i], GB[i], HA[i], HB[i], SA[i], SB[i]);
        if (i > 0) {
            e.lower[i] = mass(A[i-1], Bv[i-1], BA[i-1], BB[i-1], GA[i-1], GB[i-1], HA[i-1], HB[i-1], SA[i-1], SB[i-1]);
        }
        if (i + 1 < dim) {
            e.fm[i] = GB[i] - GB[dim-1];
            e.factor[i] = mass(max(Bv[i], cut_time), top, BB[i], BB[dim-1], GB[i], GB[dim-1], HB[i], HB[dim-1], SB[i], SB[dim-1]);
        }
    }
}

void TSP::compute_lower_sums() {
    lower_sums[0] = 0;
    double u = 0;
    for (int i = 1; i < dim; i++) {
        u += factors[i-1]*forward_probs[curr_index][i-1];
        lower_sums[i] = u*masses[i];
        assert(!isnan(lower_sums[i]));
    }
}

void TSP::compute_upper_sums() {
    partial_sum(forward_probs[curr_index].rbegin(), forward_probs[curr_index].rend()-1, upper_sums.rbegin()+1);
}

void TSP::compute_factors() {
    for (int i = 0; i < dim; i++) {
        masses[i] = cc->prob(max(curr_intervals[i]->lb, cut_time), curr_intervals[i]->ub);
    }
    prev_rho = -1;
}

void TSP::compute_trace_back_probs(int x, Interval *interval, vector<Interval *> &intervals) {
    Node *keep = pinned(interval) ? interval->node : nullptr;
    int kx = kernel_at[x];
    if (tb_kernel == kx and tb_interval == interval and tb_keep == keep) {
        return;
    }
    Kernel &K = kernels[kx];
    int k = get_interval_index(interval, intervals);
    trace_back_probs.resize(intervals.size());
    for (int i = 0; i < trace_back_probs.size(); i++) {
        trace_back_probs[i] = i == k ? K.selfs[k] : (i > k ? K.lower[k] : K.factors[i]*K.masses[k]);
        if (keep != nullptr and intervals[i]->node == keep and intervals[i]->lb == intervals[i]->ub) {
            trace_back_probs[i] = 0;
        }
    }
    tb_stay = K.stays[k];
    tb_kernel = kx;
    tb_interval = interval;
    tb_keep = keep;
}

void TSP::sanity_check(Recombination &r) {
    for (int i = 0; i < curr_intervals.size(); i++) {
        if (curr_intervals[i]->lb == curr_intervals[i]->ub and curr_intervals[i]->lb == r.inserted_node->time and curr_intervals[i]->branch != r.target_branch) {
            forward_probs[curr_index][i] = 0;
        }
    }
}

void TSP::set_state_space(int x) {
    if (n_spaces == (int) state_spaces.size()) {
        state_spaces.emplace_back();
    }
    swap(state_spaces[n_spaces - 1], state_spaces[n_spaces]);
    state_spaces[n_spaces - 1].first = x;
    state_spaces[n_spaces - 1].second = curr_intervals;
    n_spaces++;
}

vector<Interval *> &TSP::get_state_space(int x) {
    auto it = upper_bound(state_spaces.begin(), state_spaces.begin() + n_spaces, x,
                          [](int k, const pair<int, vector<Interval *>> &e) { return k < e.first; });
    --it;
    return it->second;
}

int TSP::get_interval_index(Interval *interval, vector<Interval *> &intervals) {
    return (int) distance(intervals.begin(), find(intervals.begin(), intervals.end(), interval));
}

int TSP::get_prev_breakpoint(int x) {
    auto it = upper_bound(state_spaces.begin(), state_spaces.begin() + n_spaces, x,
                          [](int k, const pair<int, vector<Interval *>> &e) { return k < e.first; });
    --it;
    return it->first;
}

Interval *TSP::sample_curr_interval(int x) {
    vector<Interval *> &intervals = get_state_space(x);
    double ws = accumulate(forward_probs[x].begin(), forward_probs[x].end(), 0.0);
    double w = ws*uniform_random();
    int i = 0;
    for (; i + 1 < (int) intervals.size(); i++) {
        w -= forward_probs[x][i];
        if (w <= 0) break;
    }
    sample_index = i;
    sel_log_q += log(forward_probs[x][i]/ws);
    return intervals[i];
}

Interval *TSP::sample_prev_interval(Interval *interval, int x) {
    vector<Interval *> &intervals = get_state_space(x);
    double ws = 0;
    compute_trace_back_probs(x + 1, interval, intervals);
    vector<double> &probs = forward_probs[x];
    ws = jump_mass(interval, intervals, probs);
    double w = ws*uniform_random();
    assert(ws > 0);
    int i = 0;
    for (; i + 1 < (int) intervals.size(); i++) {
        w -= trace_back_probs[i]*probs[i];
        if (w <= 0) break;
    }
    sample_index = i;
    sel_log_q += log(trace_back_probs[i]*probs[i]/ws);
    return intervals[i];
}

Interval *TSP::sample_source_interval(Interval *interval, int x) {
    Interval *sample_interval = interval->source_interval;
    sample_index = get_interval_index(sample_interval, get_state_space(x));
    return sample_interval;
}

Interval *TSP::sample_recomb_interval(Interval *interval, int x) {
    if (interval->lb == interval->ub) {
        return sample_curr_interval(x);
    }
    vector<Interval *> &intervals = get_state_space(x);
    double ws = 0;
    for (int i = 0; i < (int) intervals.size(); i++) {
        ws += max(own_mass(intervals[i]->time, interval->lb, interval->ub), epsilon)*forward_probs[x][i];
    }
    assert(ws > 0);
    double w = ws*uniform_random();
    int i = 0;
    for (; i + 1 < (int) intervals.size(); i++) {
        w -= max(own_mass(intervals[i]->time, interval->lb, interval->ub), epsilon)*forward_probs[x][i];
        if (w <= 0) break;
    }
    sample_index = i;
    sel_log_q += log(max(own_mass(intervals[i]->time, interval->lb, interval->ub), epsilon)*forward_probs[x][i]/ws);
    return intervals[i];
}

int TSP::trace_back_helper(Interval *interval, int x) {
    int y = get_prev_breakpoint(x);
    double non_recomb_prob = 0;
    double all_prob = 0;
    double q = uniform_random();
    double p = 1;
    double shrinkage;
    vector<Interval *> &intervals = get_state_space(x);
    while (p > q and x > y) {
        compute_trace_back_probs(x, interval, intervals);
        vector<double> &prev_probs = forward_probs[x - 1];
        non_recomb_prob = stay_mass(interval, intervals, prev_probs);
        all_prob = non_recomb_prob + jump_mass(interval, intervals, prev_probs);
        assert(all_prob > 0);
        shrinkage = non_recomb_prob/all_prob;
        assert(!isnan(shrinkage));
        p *= shrinkage;
        if (p <= q) {
            sel_log_q += ffbs_jump_logq(x, interval);
            return x;
        }
        sel_log_q += ffbs_stay_logq(x, interval);
        x -= 1;
    }
    assert(forward_probs[y][sample_index] > 0);
    return y;
}

void TSP::set_interval_constraint(Recombination &r) {
    vector<Interval *> &intervals = get_state_space(curr_index - 1);
    Interval *interval;
    for (int i = 0; i < intervals.size(); i++) {
        interval = intervals[i];
        if (interval->ub <= r.start_time) {
            forward_probs[curr_index-1][i] = 0;        } else {
            interval->lb = max(r.start_time, interval->lb);
            fill_interval_time(interval);
        }
    }
}

void TSP::set_point_constraint(Recombination &r) {
    Interval *interval;
    Interval *point_interval = search_point_interval(r);
    vector<Interval *> &intervals = get_state_space(curr_index - 1);
    for (int i = 0; i < intervals.size(); i++) {
        interval = intervals[i];
        if (interval == point_interval) {
            forward_probs[curr_index-1][i] = 1;            interval->node = r.inserted_node;
        } else {
            forward_probs[curr_index-1][i] = 0;
        }
    }
}

Interval *TSP::search_point_interval(Recombination &r) {
    double t = r.inserted_node->time;
    Interval *point_interval = nullptr;
    for (Interval *i : curr_intervals) {
        if (i->ub > t and i->lb < t) {
            point_interval = i;
        }
    }
    for (Interval *i : curr_intervals) {
        if (i->ub == i->lb and i->lb == t and i->node == r.inserted_node) {
            point_interval = i;
        }
    }
    if (point_interval != nullptr) {
        return point_interval;
    }
    vector<Interval *> candidate_point_intervals = {};
    for (Interval *i : curr_intervals) {
        if (i->lb <= t and i->ub >= t and i->lb != i->ub) {
            candidate_point_intervals.emplace_back(i);
        }
    }
    assert(candidate_point_intervals.size() == 2);
    Interval *test_interval = candidate_point_intervals[0];
    while (test_interval->source_interval != nullptr) {
        test_interval = test_interval->source_interval;
        if (test_interval->branch.upper_node == r.inserted_node or test_interval->branch.lower_node == r.inserted_node) {
            return candidate_point_intervals[1];
        }
    }
    return candidate_point_intervals[0];
}

double TSP::sample_time(double lb, double ub) {
    return cc->prob_inv(lb, ub, uniform_random());
}

bool TSP::pinned(Interval *iv) {
    return iv->node != nullptr and (iv->lb == iv->ub or pin_active);
}

double TSP::stay_mass(Interval *iv, vector<Interval *> &intervals, vector<double> &probs) {
    return tb_stay*probs[get_interval_index(iv, intervals)];
}

double TSP::jump_mass(Interval *iv, vector<Interval *> &intervals, vector<double> &probs) {
    return inner_product(trace_back_probs.begin(), trace_back_probs.end(), probs.begin(), 0.0);
}

double TSP::ffbs_stay_logq(int x, Interval *iv) {
    vector<Interval *> &pv = get_state_space(x);
    compute_trace_back_probs(x, iv, pv);
    double nr = stay_mass(iv, pv, forward_probs[x - 1]);
    double ap = nr + jump_mass(iv, pv, forward_probs[x - 1]);
    if (ap <= 0 or nr <= 0) return -numeric_limits<double>::infinity();
    return log(nr) - log(ap);
}

double TSP::ffbs_jump_logq(int x, Interval *iv) {
    vector<Interval *> &pv = get_state_space(x);
    compute_trace_back_probs(x, iv, pv);
    double wsj = jump_mass(iv, pv, forward_probs[x - 1]);
    double ap = stay_mass(iv, pv, forward_probs[x - 1]) + wsj;
    if (ap <= 0 or wsj <= 0) return -numeric_limits<double>::infinity();
    return log(wsj) - log(ap);
}

int TSP::find_old_interval(int x, const Branch &b, double t) {
    vector<Interval *> &intervals = get_state_space(x);
    for (int i = 0; i < (int) intervals.size(); i++) {
        Interval *iv = intervals[i];
        if (iv->branch != b) continue;
        if (pinned(iv)) {
            if (iv->node == find_node) return i;
        } else if (iv->lb < t and t < iv->ub) {
            return i;
        }
    }
    return -1;
}

void TSP::eval_time_at(Interval *interval, double t) {
    if (pinned(interval)) return;
    time_log_q += cc->log_density(interval->lb, interval->ub, t);
}

double TSP::eval_joining_nodes(map<double, Branch> &old_jb, map<double, Branch> &old_ab,
                              int start_index, vector<double> &coordinates) {
    forward_to(curr_index);
    sel_log_q = 0;
    time_log_q = 0;
    auto old_branch_at = [&](int xx) {
        auto it = old_jb.upper_bound(coordinates[xx + start_index]);
        if (it != old_jb.begin()) --it;
        return it->second;
    };
    auto old_time_at = [&](int xx) {
        auto it = old_ab.upper_bound(coordinates[xx + start_index]);
        if (it != old_ab.begin()) --it;
        return it->second.upper_node->time;
    };
    auto old_node_at = [&](int xx) {
        auto it = old_ab.upper_bound(coordinates[xx + start_index]);
        if (it != old_ab.begin()) --it;
        return it->second.upper_node;
    };
    double neg_inf = -numeric_limits<double>::infinity();
    auto eval_curr_sel = [&](int xx, int idx) {
        double ws = accumulate(forward_probs[xx].begin(), forward_probs[xx].end(), 0.0);
        if (ws > 0 and forward_probs[xx][idx] > 0) {
            sel_log_q += log(forward_probs[xx][idx] / ws);
        } else {
            sel_log_q = neg_inf;
        }
    };
    int x = curr_index;
    double t = old_time_at(x);
    find_node = old_node_at(x);
    pin_active = false;
    int oi = find_old_interval(x, old_branch_at(x), t);
    if (oi < 0) return neg_inf;
    Interval *interval = get_state_space(x)[oi];
    eval_curr_sel(x, oi);
    sample_index = oi;
    eval_time_at(interval, t);
    while (true) {
        int y = get_prev_breakpoint(x);
        Node *cur_node = old_node_at(x);
        int xstop = y;
        for (int xx = x; xx > y; xx--) {
            if (old_node_at(xx - 1) != cur_node) { xstop = xx; break; }
        }
        for (int xx = x; xx > xstop; xx--) {
            sel_log_q += ffbs_stay_logq(xx, interval);
        }
        x = xstop;
        if (x == 0) {
            break;
        }
        if (x == interval->start_pos and interval->source_interval != nullptr) {
            if (old_node_at(x - 1) != cur_node) return neg_inf;
            x -= 1;
            interval = interval->source_interval;
            vector<Interval *> &pv = get_state_space(x);
            sample_index = get_interval_index(interval, pv);
            rewind_to(x);
        } else if (x == interval->start_pos) {
            x -= 1;
            double nt = old_time_at(x);
            find_node = old_node_at(x);
            pin_active = true;
            int ni = find_old_interval(x, old_branch_at(x), nt);
            if (ni < 0) return neg_inf;
            vector<Interval *> &pv = get_state_space(x);
            pin_active = pv[ni]->node != nullptr;
            if (!pinned(pv[ni]) and find_node == cur_node) return neg_inf;
            if (interval->lb == interval->ub) {
                eval_curr_sel(x, ni);
            } else {
                double ws = 0;
                for (int i = 0; i < (int) pv.size(); i++)
                    ws += max(own_mass(pv[i]->time, interval->lb, interval->ub), epsilon) * forward_probs[x][i];
                double num = max(own_mass(pv[ni]->time, interval->lb, interval->ub), epsilon) * forward_probs[x][ni];
                if (ws > 0 and num > 0) {
                    sel_log_q += log(num / ws);
                } else {
                    sel_log_q = neg_inf;
                }
            }
            rewind_to(x);
            sample_index = ni;
            interval = pv[ni];
            eval_time_at(interval, nt);
        } else {
            sel_log_q += ffbs_jump_logq(x, interval);
            x -= 1;
            vector<Interval *> &pv = get_state_space(x);
            compute_trace_back_probs(x + 1, interval, pv);
            double nt = old_time_at(x);
            find_node = old_node_at(x);
            pin_active = false;
            int ni = find_old_interval(x, old_branch_at(x), nt);
            if (ni < 0) return neg_inf;
            double ws = jump_mass(interval, pv, forward_probs[x]);
            double num = trace_back_probs[ni] * forward_probs[x][ni];
            if (ws > 0 and num > 0) {
                sel_log_q += log(num / ws);
            } else {
                sel_log_q = neg_inf;
            }
            rewind_to(x);
            sample_index = ni;
            interval = pv[ni];
            eval_time_at(interval, nt);
        }
    }
    return sel_log_q + time_log_q;
}

Node *TSP::sample_joining_node(Interval *interval) {
    Node *n = nullptr;
    double t;
    if (pinned(interval)) {
        n = interval->node;
    } else {
        t = sample_time(interval->lb, interval->ub);
        time_log_q += cc->log_density(interval->lb, interval->ub, t);
        node_owner.push_back(new_node(t));
        n = node_owner.back().get();
        n->set_index(counter);
        counter += 1;
    }
    return n;
}

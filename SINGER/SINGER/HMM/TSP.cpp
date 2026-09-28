//
//  TSP.cpp
//  SINGER
//
//  Created by Yun Deng on 9/25/23.
//

#include "TSP.hpp"

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
    x.reduction = 1.0;
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
    state_spaces.clear();
    state_spaces.push_back({INT_MAX, {}});
    n_rows = 0;
    curr_index = 0;
    curr_branch = Branch();
    curr_intervals.clear();
    rhos.clear();
    sister_masses.clear();
    thetas.clear();
    lower_sums.clear();
    upper_sums.clear();
    diagonals.clear();
    factors.clear();
    lower_diagonals.clear();
    upper_diagonals.clear();
    temp.clear();
    null_emit_probs.clear();
    mut_emit_probs.clear();
    trace_back_probs.clear();
    check_points.clear();
    prev_rho = -1;
    prev_theta = -1;
    prev_node = nullptr;
    tb_generation = -1;
    cc_generation = 0;
    dim = 0;
    sample_index = -1;
    lower_bound = 0;
    n_arena = 0;
}

double TSP::per_length(double rho, double L) {
    return -expm1(-rho*L)/L;
}

void TSP::reserve_memory(int length) {
    forward_probs.reserve(length);
    rhos.reserve(length);
}

void TSP::start(Branch &branch, double t) {
    cut_time = t;
    cd_valid = false;
    cc_lowest_change = numeric_limits<double>::infinity();
    curr_index = 0;
    curr_branch = branch;
    lower_bound = max(cut_time, branch.lower_node->time);
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
    rhos.emplace_back(0);
    sister_masses.emplace_back(sister_mass);
    prev_rho = -1;
    prev_theta = -1;
    prev_node = nullptr;
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
    cc_generation += 1;
    curr_index += 1;
    curr_branch = next_branch;
    lower_bound = max(cut_time, next_branch.lower_node->time);
    if (prev_branch == r.source_branch and next_branch == r.merging_branch) {
        set_interval_constraint(r);
    } else if (prev_branch == r.target_branch and next_branch == r.recombined_branch) {
        set_point_constraint(r);
    }
    curr_intervals.clear();
    if (prev_branch == r.source_branch and next_branch == r.merging_branch) { // switch to a point mass
        double t = r.deleted_node->time;
        generate_intervals(next_branch, next_branch.lower_node->time, t);
        generate_intervals(next_branch, t, t);
        temp.back() = 1.0f; // set point mass prob marker
        curr_intervals.back()->node = r.deleted_node; // set point mass node marker
        generate_intervals(next_branch, t, next_branch.upper_node->time);
    } else if (prev_branch == r.target_branch and next_branch == r.recombined_branch) { // switch from a point mass
        generate_intervals(next_branch, next_branch.lower_node->time, r.start_time);
        generate_intervals(next_branch, r.start_time, next_branch.upper_node->time);
        for (int i = 0; i < curr_intervals.size(); i++) {
            if (curr_intervals[i]->lb >= r.start_time) {
                temp[i] = 1.0;
            }
        }
    } else {
        double lb;
        double ub;
        lb = next_branch.lower_node->time;
        ub = max(prev_branch.lower_node->time, next_branch.lower_node->time);
        generate_intervals(next_branch, lb, ub);
        transfer_intervals(r, prev_branch, next_branch);
        lb = min(curr_intervals.back()->ub, next_branch.upper_node->time);
        ub = next_branch.upper_node->time;
        generate_intervals(next_branch, lb, ub);
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
    rhos.emplace_back(0);
    sister_masses.emplace_back(sister_mass);
    prev_rho = -1;
    prev_theta = -1;
    prev_node = nullptr;
    curr_branch = next_branch;
    curr_index += 1;
    lower_bound = max(cut_time, next_branch.lower_node->time);
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

double TSP::get_exp_quantile(double p) {
    assert(p >= 0 and p <= 1);
    if (p < 1e-6) {
        return 0;
    }
    if (1 - p < 1e-6) {
        return numeric_limits<double>::infinity();
    }
    return -log(1 - p);
}

vector<double> TSP::generate_grid(double lb, double ub) {
    assert(lb < ub);
    vector<double> points = {lb};
    double ls = cc->surv(lb);
    double us = cc->surv(ub);
    double q = ls - us;
    int n = max((int) ceil((1 - us/ls)/gap), min_num);
    double l;
    for (int i = 1; i < n; i++) {
        l = cc->surv_inv(ls - i*q/n);
        points.emplace_back(l);
    }
    points.emplace_back(ub);
    return points;
}

double TSP::recomb_cdf(double s, double t) {
    if (isinf(t)) {
        return 1;
    }
    if (t == 0) {
        return 0;
    }
    double cdf = cc->recomb_mass(s, t)/(s - cut_time);
    assert(!isnan(cdf));
    return cdf;
}

void TSP::forward(double rho) {
    rhos.emplace_back(rho);
    sister_masses.emplace_back(sister_mass);
    compute_diagonals(rho);
    compute_lower_diagonals(rho);
    compute_upper_diagonals(rho);
    compute_lower_sums();
    compute_upper_sums();
    curr_index += 1;
    prev_rho = rho;
    push_row(lower_sums);
    for (int i = 0; i < dim; i++) {
        assert(forward_probs[curr_index][i] >= 0);
        forward_probs[curr_index][i] += diagonals[i]*forward_probs[curr_index-1][i] + lower_diagonals[i]*upper_sums[i];
        if (curr_intervals[i]->lb != curr_intervals[i]->ub or forward_probs[curr_index][i] > 0) {
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
    if (ws > 0) {
        double inv_ws = 1.0/ws;
        for (int j = 0; j < dim; j++) {
            curr_probs[j] *= inv_ws;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            curr_probs[i] = 1.0/forward_probs[curr_index].size();
        }
    }
}

void TSP::mut_emit(double theta, double bin_size, vector<double> &mut_set, Node *query_node) {
    compute_mut_emit_probs(theta, bin_size, mut_set, query_node);
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
    assert(ws > 0);
    double inv_ws = 1.0/ws;
    for (int j = 0; j < dim; j++) {
        curr_probs[j] *= inv_ws;
    }
}

map<double, Node *> TSP::sample_joining_nodes(int start_index, vector<double> &coordinates) {
    prev_rho = -1;
    log_h = 0;
    sel_log_q = 0;
    time_log_q = 0;
    map<double, Node *> joining_nodes = {};
    int x = curr_index;
    double pos = coordinates[x + start_index + 1];
    pin_active = false;
    Interval *interval = sample_curr_interval(x);
    Node *n = sample_joining_node(interval);
    joining_nodes[pos] = nullptr;
    while (x >= 0) {
        x = trace_back_helper(interval, x);
        pos = coordinates[x + start_index];
        joining_nodes[pos] = n;
        assert(x >= interval->start_pos);
        if (x == 0) {
            break;
        } if (x == interval->start_pos) {
            if (interval->source_interval != nullptr) {
                x -= 1;
                interval = sample_source_interval(interval, x);
            } else {
                x -= 1;
                interval = sample_recomb_interval(interval, x);
                pin_active = interval->node != nullptr;
                n = sample_joining_node(interval);
            }
        } else {
            x -= 1;
            interval = sample_prev_interval(interval, x);
            pin_active = false;
            n = sample_joining_node(interval);
        }
        prev_rho = -1;
    }
    return joining_nodes;
}

double TSP::non_recomb_prob(double rho, double s) {
    return exp(-rho*(s - cut_time));
}

double TSP::sister_factor(double s) {
    return 2.0;
}

void TSP::fill_interval_time(Interval *iv) {
    iv->s_lb = cc->surv(iv->lb);
    iv->s_ub = cc->surv(iv->ub);
    if (iv->ub - iv->lb < 1e-3) {
        iv->fill_time();
        return;
    }
    iv->time = cc->surv_inv(0.5*(iv->s_lb + iv->s_ub));
}

double TSP::recomb_prob(double s, double t1, double t2) {
    assert(t1 <= t2);
    assert(t1 >= cut_time and s >= cut_time);
    if (t1 == t2) {
        return 0;
    }
    if (s - cut_time < 0.005) {
        return max(epsilon, cc->surv(t1) - cc->surv(t2));
    }
    double pl = recomb_cdf(s, t1);
    double pu = recomb_cdf(s, t2);
    double p = pu - pl;
    assert(pu >= pl);
    p = max(p, epsilon);
    return p;
}

double TSP::psmc_cdf(double rho, double s, double t) {
    if (s <= cut_time) {
        return 0;
    }
    return -expm1(-rho*(s - cut_time))*recomb_cdf(s, t);
}

double TSP::standard_recomb_cdf(double rho, double s, double t) {
    double l = 2*s;
    double integral;
    double cdf;
    if (t <= s) {
        integral = t + 0.5*exp(-2*t) - 0.5;
    } else {
        integral = s + 0.5*exp(-2*s) - 0.5 + (1 - exp(s-t))*(1 - exp(-2*s));
    }
    cdf = integral/l;
    return cdf;
}

double TSP::psmc_prob(double rho, double s, double t1, double t2) {
    assert(s != numeric_limits<double>::infinity());
    assert(t1 <= t2);
    assert(t1 >= lower_bound and s >= lower_bound);
    if (t1 == t2) {
        return 0;
    }
    double prob;
    double uq = 0;
    double lq = 0;
    uq = psmc_cdf(rho, s, t2);
    lq = psmc_cdf(rho, s, t1);
    prob = uq - lq;
    assert(!isnan(prob));
    prob = max(prob, epsilon);
    assert(prob <= 1);
    return prob;
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
            double w = (cc1->recomb_mass(t, ub) - cc1->recomb_mass(t, t))/den*cc1->surv(t);
            m += w*cc->prob(max(lb, t), ub)/cc->surv(t);
        }
    }
    return m;
}

double TSP::jump_prob(double rho, double s, double t1, double t2) {
    double lb = max(t1, cut_time);
    if (t2 <= lb) {
        return epsilon;
    }
    double p;
    double fs = sister_factor(s);
    if (s - cut_time < 0.005) {
        p = fs*per_length(rho, s - cut_time)*(s - cut_time)*(cc->surv(lb) - cc->surv(t2));
    } else {
        p = fs*per_length(rho, s - cut_time)*own_mass(s, lb, t2);
    }
    if (sister_mass > 0) {
        double q = 0;
        double tb = min(t2, s);
        if (tb > lb) {
            q += cc1->prob(lb, tb);
        }
        if (t2 > s) {
            q += exp(-(s - cut_time))*cc->prob(max(lb, s), t2);
        }
        p += rho*sister_mass*q;
    }
    return max(p, epsilon);
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
    double l;
    double u;
    for (int i = 0; i < points.size() - 1; i++) {
        l = points[i];
        u = points[i+1];
        new_interval = make_interval(next_branch, l, u, curr_index);
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
            assert(!isnan(p));
            new_interval = make_interval(next_branch, lb, ub, curr_index);
            fill_interval_time(new_interval);
            new_interval->node = interval->node;
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
    upper_diagonals.resize(dim); upper_diagonals.assign(dim, 0);
    lower_sums.resize(dim); lower_sums.assign(dim, 0);
    upper_sums.resize(dim); upper_sums.assign(dim, 0);
    null_emit_probs.resize(dim); null_emit_probs.assign(dim, 0);
    mut_emit_probs.resize(dim); mut_emit_probs.assign(dim, 0);
    factors.resize(dim); factors.assign(dim, 0);
    masses.resize(dim); masses.assign(dim, 0);
}

double TSP::random() {
    double p = uniform_random();
    return p;
}

double TSP::get_prop(double lb1, double ub1, double lb2, double ub2) {
    double p;
    if (ub2 - lb2 < 1e-6) {
        p = 1;
    } else {
        double p1 = cc->surv(lb1) - cc->surv(ub1);
        double p2 = cc->surv(lb2) - cc->surv(ub2);
        p = p1/p2;
    }
    return p;
}

void TSP::compute_null_emit_probs(double theta, Node *query_node) {
    if (theta == prev_theta and query_node == prev_node) {
        return;
    }
    for (int i = 0; i < dim; i++) {
        null_emit_probs[i] = eh->null_emit(curr_branch, curr_intervals[i]->time, theta, query_node);
    }
}

void TSP::compute_mut_emit_probs(double theta, double bin_size, vector<double> &mut_set, Node *query_node) {
    for (int i = 0; i < dim; i++) {
        mut_emit_probs[i] = eh->mut_emit(curr_branch, curr_intervals[i]->time, theta, bin_size, mut_set, query_node);
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
            upper_diagonals = cd_upper;
            factors = cd_factors;
            return;
        }
    }
    vector<double> A(dim), Bv(dim), BA(dim), BB(dim), GA(dim), GB(dim), HA(dim), HB(dim), SA(dim), SB(dim);
    for (int j = 0; j < dim; j++) {
        A[j] = max(curr_intervals[j]->lb, cut_time);
        Bv[j] = curr_intervals[j]->ub;
        BA[j] = cc1->recomb_mass(A[j] + 1, A[j]);
        BB[j] = cc1->recomb_mass(Bv[j] + 1, Bv[j]);
        GA[j] = cc->prob(cut_time, A[j]);
        GB[j] = cc->prob(cut_time, Bv[j]);
        HA[j] = cc1->prob(cut_time, A[j]);
        HB[j] = cc1->prob(cut_time, Bv[j]);
        SA[j] = cc->surv(A[j]);
        SB[j] = cc->surv(Bv[j]);
    }
    double top = Bv[dim-1];
    double a0 = A[0];
    for (int i = 0; i < dim; i++) {
        double t = curr_intervals[i]->time;
        double stay;
        stay = non_recomb_prob(rho, t);
        double fs = sister_factor(t);
        double Bt = cc1->recomb_mass(t + 1, t);
        double Gt = cc->prob(cut_time, t);
        double Ht = cc1->prob(cut_time, t);
        double S = cc->surv(t);
        double W = 0;
        if (top > t) {
            double den = cc1->prob(t, top);
            if (den > 0) W = (cc1->recomb_mass(t, top) - cc1->recomb_mass(t, t))/den*cc1->surv(t);
        }
        bool small = t - cut_time < 0.005;
        double pref = small ? 0 : per_length(rho, t - cut_time);
        double pref0 = per_length(rho, t - cut_time)*(t - cut_time);
        double sis = exp(-(t - cut_time));
        auto jump = [&](double a, double b, double Ba, double Bb, double Ga, double Gb, double Ha, double Hb, double Sa, double Sb) {
            if (b <= a) return epsilon;
            double tb = min(b, t);
            double p;
            if (small) {
                p = fs*pref0*(Sa - Sb);
            } else {
                double m = 0;
                if (tb > a) m += (tb == b ? Bb : Bt) - Ba;
                if (b > t) m += W*(Gb - (a >= t ? Ga : Gt))/S;
                p = fs*pref*m;
            }
            if (sister_mass > 0) {
                double q = 0;
                if (tb > a) q += (tb == b ? Hb : Ht) - Ha;
                if (b > t) q += sis*(Gb - (a >= t ? Ga : Gt));
                p += rho*sister_mass*q;
            }
            return max(p, epsilon);
        };
        double full = jump(a0, top, BA[0], BB[dim-1], GA[0], GB[dim-1], HA[0], HB[dim-1], SA[0], SB[dim-1]);
        double base = stay + full;
        double self = jump(A[i], Bv[i], BA[i], BB[i], GA[i], GB[i], HA[i], HB[i], SA[i], SB[i]);
        diagonals[i] = (stay + self)/base;
        if (i > 0) {
            lower_diagonals[i-1] = jump(A[i-1], Bv[i-1], BA[i-1], BB[i-1], GA[i-1], GB[i-1], HA[i-1], HB[i-1], SA[i-1], SB[i-1])/base;
        }
        if (i + 1 < dim) {
            upper_diagonals[i+1] = jump(A[i+1], Bv[i+1], BA[i+1], BB[i+1], GA[i+1], GB[i+1], HA[i+1], HB[i+1], SA[i+1], SB[i+1])/base;
            double a = max(Bv[i], cut_time);
            double m = GB[dim-1] - GB[i];
            factors[i] = (m > 0) ? jump(a, top, BB[i], BB[dim-1], GB[i], GB[dim-1], HB[i], HB[dim-1], SB[i], SB[dim-1])/base/m : 0;
        }
        assert(!isnan(diagonals[i]));
    }
    lower_diagonals[dim-1] = 0;
    upper_diagonals[0] = 0;
    factors[dim-1] = 0;
    cd_lb.resize(dim); cd_ub.resize(dim); cd_time.resize(dim);
    for (int i = 0; i < dim; i++) {
        cd_lb[i] = curr_intervals[i]->lb;
        cd_ub[i] = curr_intervals[i]->ub;
        cd_time[i] = curr_intervals[i]->time;
    }
    cd_diagonals = diagonals;
    cd_lower = lower_diagonals;
    cd_upper = upper_diagonals;
    cd_factors = factors;
    cd_rho = rho;
    cd_sister_mass = sister_mass;
    cd_top = top;
    cc_lowest_change = numeric_limits<double>::infinity();
    cd_valid = true;
}

void TSP::compute_lower_diagonals(double rho) {
}

void TSP::compute_upper_diagonals(double rho) {
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

void TSP::compute_emissions(vector<double> &mut_set, const Branch &branch, Node *node) {
    fill(emissions.begin(), emissions.end(), 0);
    double sl, su, s0, sm = 0;
    for (double x : mut_set) {
        sl = branch.lower_node->get_state(x);
        su = branch.upper_node->get_state(x);
        s0 = node->get_state(x);
        if (sl + su + s0 > 1.5) {
            sm = 1;
        } else {
            sm = 0;
        }
        emissions[0] += abs(sm - sl);
        emissions[1] += abs(sm - su);
        emissions[2] += abs(sm - s0);
        emissions[3] += abs(sl - su);
    }
}

void TSP::compute_trace_back_probs(double rho, Interval *interval, vector<Interval *> &intervals) {
    if (rho == prev_rho) {
        return;
    }
    if (tb_generation == cc_generation and tb_interval == interval and tb_lb == interval->lb
        and tb_ub == interval->ub and tb_sister_mass == sister_mass
        and tb_states == (const Interval *const *) intervals.data() and tb_nstates == intervals.size()
        and trace_back_probs.size() == intervals.size() and rho == tb_rho) {
        prev_rho = rho;
        return;
    }
    for (int i = 0; i < trace_back_probs.size(); i++) {
        trace_back_probs[i] = jump_prob(rho, intervals[i]->time, interval->lb, interval->ub);
    }
    tb_generation = cc_generation;
    tb_interval = interval;
    tb_lb = interval->lb;
    tb_ub = interval->ub;
    tb_sister_mass = sister_mass;
    tb_states = (const Interval *const *) intervals.data();
    tb_nstates = intervals.size();
    tb_rho = rho;
}

void TSP::sanity_check(Recombination &r) {
    for (int i = 0; i < curr_intervals.size(); i++) {
        if (curr_intervals[i]->lb == curr_intervals[i]->ub and curr_intervals[i]->lb == r.inserted_node->time and curr_intervals[i]->branch != r.target_branch) {
            forward_probs[curr_index][i] = 0;
        }
        /*
        if (curr_intervals[i]->lb == curr_intervals[i]->ub and curr_intervals[i]->lb == r.inserted_node->time and curr_intervals[i]->node != r.deleted_node) {
            forward_probs[curr_index][i] = 0;
        }
         */
    }
}

void TSP::set_state_space(int x) {
    state_spaces.insert(prev(state_spaces.end()), {x, curr_intervals});
}

vector<Interval *> &TSP::get_state_space(int x) {
    auto it = upper_bound(state_spaces.begin(), state_spaces.end(), x,
                          [](int k, const pair<int, vector<Interval *>> &e) { return k < e.first; });
    --it;
    return it->second;
}

int TSP::get_interval_index(Interval *interval, vector<Interval *> &intervals) {
    auto it = find(intervals.begin(), intervals.end(), interval);
    int index = (int) distance(intervals.begin(), it);
    return index;
}

int TSP::get_prev_breakpoint(int x) {
    auto it = upper_bound(state_spaces.begin(), state_spaces.end(), x,
                          [](int k, const pair<int, vector<Interval *>> &e) { return k < e.first; });
    --it;
    return it->first;
}

Interval *TSP::sample_curr_interval(int x) {
    vector<Interval *> &intervals = get_state_space(x);
    double ws = accumulate(forward_probs[x].begin(), forward_probs[x].end(), 0.0);
    double q = random();
    double w = ws*q;
    for (int i = 0; i < intervals.size(); i++) {
        w -= forward_probs[x][i];
        if (w <= 0) {
            sample_index = i;
            sel_log_q += log(forward_probs[x][i]/ws);
            return intervals[i];
        }
    }
    cerr << "tsp sample curr interval failed" << endl;
    exit(1);
}

Interval *TSP::sample_prev_interval(Interval *interval, int x) {
    vector<Interval *> &intervals = get_state_space(x);
    lower_bound = intervals.front()->lb;
    double ws = 0;
    double rho = rhos[x];
    sister_mass = sister_masses[x];
    prev_rho = -1;
    compute_trace_back_probs(rho, interval, intervals);
    vector<double> &probs = forward_probs[x];
    ws = jump_mass(interval, intervals, probs);
    double q = random();
    double w = ws*q;
    assert(ws > 0);
    for (int i = 0; i < intervals.size(); i++) {
        w -= trace_back_probs[i]*probs[i];
        if (w <= 0) {
            sample_index = i;
            sel_log_q += log(trace_back_probs[i]*probs[i]/ws);
            return intervals[i];
        }
    }
    cerr << "tsp prev sampling failed" << endl;
    exit(1);
}

Interval *TSP::sample_source_interval(Interval *interval, int x) {
    Interval *sample_interval = interval->source_interval;
    vector<Interval *> &intervals = get_state_space(x);
    sample_index = get_interval_index(sample_interval, intervals);
    return sample_interval;
}

Interval *TSP::sample_recomb_interval(Interval *interval, int x) {
    if (interval->lb == interval->ub) { // handles point mass case nicely
        return sample_curr_interval(x);
    }
    vector<Interval *> &intervals = get_state_space(x);
    double ws = 0;
    Interval *prev_interval = nullptr;
    for (int i = 0; i < intervals.size(); i++) {
        prev_interval = intervals[i];
        ws += max(own_mass(prev_interval->time, interval->lb, interval->ub), epsilon)*forward_probs[x][i];
    }
    assert(ws > 0);
    double q = random();
    double w = ws*q;
    for (int i = 0; i< intervals.size(); i++) {
        prev_interval = intervals[i];
        w -= max(own_mass(prev_interval->time, interval->lb, interval->ub), epsilon)*forward_probs[x][i];
        if (w <= 0) {
            sample_index = i;
            sel_log_q += log(max(own_mass(intervals[i]->time, interval->lb, interval->ub), epsilon)*forward_probs[x][i]/ws);
            return intervals[i];
        }
    }
    cerr << "tsp recomb sampling failed" << endl;
    exit(1);
}

int TSP::trace_back_helper(Interval *interval, int x) {
    int y = get_prev_breakpoint(x);
    double non_recomb_prob = 0;
    double all_prob = 0;
    double q = random();
    double p = 1;
    double shrinkage;
    double rho;
    vector<Interval *> &intervals = get_state_space(x);
    lower_bound = intervals.front()->lb;
    if (trace_back_probs.size() != intervals.size()) trace_back_probs.assign(intervals.size(), 0.0);
    if (x > y) {
        sister_mass = sister_masses[x-1];
        prev_rho = -1;
    }
    while (p > q and x > y) {
        rho = rhos[x-1];
        compute_trace_back_probs(rho, interval, intervals);
        prev_rho = rho;
        vector<double> &prev_probs = forward_probs[x - 1];
        non_recomb_prob = stay_mass(interval, intervals, prev_probs, rho);
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
            forward_probs[curr_index-1][i] = 0; // curr_index - 1 because the index has already moved forward
        } else {
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
            forward_probs[curr_index-1][i] = 1; // curr_index - 1 because the index has already moved forward
            interval->node = r.inserted_node;
        } else {
            forward_probs[curr_index-1][i] = 0;
        }
    }
}

Interval *TSP::search_point_interval(Recombination &r) {
    double t = r.inserted_node->time;
    Interval *point_interval = nullptr;
    for (Interval *i : curr_intervals) { // cross interval is always valid
        if (i->ub > t and i->lb < t) {
            point_interval = i;
        }
    }
    for (Interval *i : curr_intervals) {
        if (i->ub == i->lb and i->lb == t and i->node == r.inserted_node) { // point interval valid only the node is correct
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
    double ls = cc->surv(lb);
    double us = cc->surv(ub);
    double m = cc->surv_inv(ls - uniform_random()*(ls - us));
    if (m < lb) m = lb;
    if (m > ub) m = ub;
    return m;
}

double TSP::exp_median(double lb, double ub) {
    assert(lb <= ub);
    if (ub - lb <= 0.005) {
        return (0.45 + 0.1*random())*(ub - lb) + lb;
    }
    double ls = cc->surv(lb);
    double us = cc->surv(ub);
    double m = cc->surv_inv(ls - (0.45 + 0.1*random())*(ls - us));
    assert(m >= lb and m <= ub);
    return m;
}


bool TSP::pinned(Interval *iv) {
    return iv->node != nullptr and (iv->lb == iv->ub or pin_active);
}

double TSP::stay_mass(Interval *iv, vector<Interval *> &intervals, vector<double> &probs, double rho) {
    int fi = get_interval_index(iv, intervals);
    return non_recomb_prob(rho, iv->time)*probs[fi];
}

double TSP::jump_mass(Interval *iv, vector<Interval *> &intervals, vector<double> &probs) {
    return inner_product(trace_back_probs.begin(), trace_back_probs.end(), probs.begin(), 0.0);
}

double TSP::ffbs_stay_logq(int x, Interval *iv) {
    vector<Interval *> &pv = get_state_space(x);
    if (trace_back_probs.size() != pv.size()) trace_back_probs.assign(pv.size(), 0.0);
    sister_mass = sister_masses[x - 1];
    prev_rho = -1;
    compute_trace_back_probs(rhos[x - 1], iv, pv);
    double nr = stay_mass(iv, pv, forward_probs[x - 1], rhos[x - 1]);
    double ap = nr + jump_mass(iv, pv, forward_probs[x - 1]);
    if (ap <= 0 or nr <= 0) return -numeric_limits<double>::infinity();
    return log(nr) - log(ap);
}

double TSP::ffbs_jump_logq(int x, Interval *iv) {
    vector<Interval *> &pv = get_state_space(x);
    if (trace_back_probs.size() != pv.size()) trace_back_probs.assign(pv.size(), 0.0);
    sister_mass = sister_masses[x - 1];
    prev_rho = -1;
    compute_trace_back_probs(rhos[x - 1], iv, pv);
    double nr = stay_mass(iv, pv, forward_probs[x - 1], rhos[x - 1]);
    double wsj = jump_mass(iv, pv, forward_probs[x - 1]);
    double ap = nr + wsj;
    if (ap <= 0 or wsj <= 0) return -numeric_limits<double>::infinity();
    return log(wsj) - log(ap);
}

Node *TSP::sample_joining_node(Interval *interval) {
    Node *n = nullptr;
    double t;
    if (pinned(interval)) {
        n = interval->node;
    } else {
        if (interval->ub - interval->lb > 0.005) {
            double sd = cc->surv(interval->lb) - cc->surv(interval->ub);
            if (sd > 0) log_h += log(sd);
        }
        t = sample_time(interval->lb, interval->ub);
        {
            double sv = cc->surv(t);
            double sd = cc->surv(interval->lb) - cc->surv(interval->ub);
            if (sd > 0 and sv > 0) time_log_q += log(sv) - log(sd) + log(cc->rate(t));
        }
        node_owner.push_back(new_node(t));
        n = node_owner.back().get();
        n->set_index(counter);
        counter += 1;
    }
    assert(n != nullptr);
    return n;
}

//
//  RSP_smc.cpp
//  SINGER
//
//  Created by Yun Deng on 3/13/23.
//

#include "RSP_smc.hpp"

RSP_smc::RSP_smc() {
}

void RSP_smc::set_tree(Tree &tree) {
    vector<double> node_times;
    for (auto &x : tree.parents) {
        if (x.first->time > 0) {
            node_times.push_back(x.first->time);
        }
    }
    sort(node_times.begin(), node_times.end());
    level_times.assign(1, 0.0);
    for (double t : node_times) {
        if (t != level_times.back()) {
            level_times.push_back(t);
        }
    }
    int m = (int) level_times.size();
    level_rates.resize(m);
    level_lambda.assign(m, 0.0);
    for (int j = 0; j < m; j++) {
        level_rates[j] = 1.0 + (node_times.end() - upper_bound(node_times.begin(), node_times.end(), level_times[j]));
        if (j > 0) {
            level_lambda[j] = level_lambda[j-1] + level_rates[j-1]*(level_times[j] - level_times[j-1]);
        }
    }
}

int RSP_smc::level_of(double s) {
    const double *base = level_times.data();
    int len = (int) level_times.size();
    while (len > 1) {
        int half = len/2;
        base = (base[half] <= s) ? base + half : base;
        len -= half;
    }
    return (int) (base - level_times.data());
}

double RSP_smc::lambda(double s) {
    int j = level_of(s);
    return level_lambda[j] + level_rates[j]*(s - level_times[j]);
}

double RSP_smc::exp_lambda_integral(double lo, double hi) {
    double lam_hi = lambda(hi);
    double total = 0;
    int j = level_of(lo);
    double a = lo;
    while (a < hi) {
        double b = (j + 1 < (int) level_times.size()) ? min(level_times[j+1], hi) : hi;
        double k = level_rates[j];
        total += (exp(lambda(b) - lam_hi) - exp(lambda(a) - lam_hi))/k;
        a = b;
        j += 1;
    }
    return total;
}

double RSP_smc::sister_mass(double lo, double hi) {
    double lam_hi = lambda(hi) + hi;
    double total = 0;
    int j = level_of(lo);
    double a = lo;
    while (a < hi) {
        double b = (j + 1 < (int) level_times.size()) ? min(level_times[j+1], hi) : hi;
        double k = level_rates[j] + 1;
        total += (exp(lambda(b) + b - lam_hi) - exp(lambda(a) + a - lam_hi))/k;
        a = b;
        j += 1;
    }
    return total;
}

double RSP_smc::sample_start_time(Branch b, int density, double join_time, double cut_time) {
    double lb = b.lower_node->time;
    double ub = b.upper_node->time;
    lb = max(cut_time, lb);
    ub = min(join_time, ub);
    assert(lb < ub);
    double p;
    double t;
    double w;
    double start_time = 0;
    vector<double> start_times = {};
    vector<double> weights = {};
    double weight_sum = 0;
    for (int i = 0; i < density; i++) {
        t = random_time(lb, ub);
        w = recomb_pdf(t, ub);
        start_times.push_back(t);
        weights.push_back(w);
        weight_sum += w;
    }
    p = uniform_random();
    weight_sum = weight_sum*p;
    for (int i = 0; i < density; i++) {
        weight_sum -= weights[i];
        if (weight_sum <= 0) {
            start_time = start_times[i];
            break;
        }
    }
    assert(start_time >= cut_time);
    assert(start_time <= join_time);
    assert(start_time <= ub);
    return start_time;
}

pair<Branch, double> RSP_smc::sample_start_time(Branch b1, Branch b2, int density, double join_time, double cut_time) {
    double lb1 = max(b1.lower_node->time, cut_time);
    double ub1 = min(b1.upper_node->time, join_time);
    double lb2 = max(b2.lower_node->time, cut_time);
    double ub2 = min(b2.upper_node->time, join_time);
    assert(lb1 < ub1);
    assert(lb2 < ub2);
    double q = (ub1 - lb1)/(ub1 + ub2 - lb1 - lb2);
    int n1 = round(density*q);
    int n2 = density - n1;
    double p;
    double t;
    double w;
    Branch source_branch;
    double start_time = 0;
    vector<double> start_times(0);
    vector<double> weights(0);
    vector<int> branch_indices(0);
    double weight_sum = 0;
    for (int i = 0; i < n1; i++) {
        t = random_time(lb1, ub1);
        w = recomb_pdf(t, ub1);
        start_times.push_back(t);
        weights.push_back(w);
        weight_sum += w;
        branch_indices.push_back(1);
    }
    for (int i = 0; i < n2; i++) {
        t = random_time(lb2, ub2);
        w = recomb_pdf(t, ub2);
        start_times.push_back(t);
        weights.push_back(w);
        weight_sum += w;
        branch_indices.push_back(2);
    }
    p = uniform_random();
    weight_sum = weight_sum*p;
    for (int i = 0; i < density; i++) {
        weight_sum -= weights[i];
        if (weight_sum <= 0) {
            start_time = start_times[i];
            if (branch_indices[i] == 1) {
                source_branch = b1;
            } else {
                source_branch = b2;
            }
            break;
        }
    }
    assert(weight_sum <= 0);
    assert(start_time >= cut_time);
    assert(start_time <= join_time);
    assert(start_time <= min(ub1, ub2));
    return {source_branch, start_time};
}

void RSP_smc::sample_recombination(Recombination &r, double cut_time, Tree &tree) {
    if (r.pos == 0) {
        return;
    }
    if (r.start_time > 0) {
        return;
    }
    if (r.deleted_branches.size() == 0) {
        return;
    }
    get_coalescence_rate(tree, r, cut_time);
    vector<Branch> source_candidates;
    for (Branch b : r.deleted_branches) {
        if (b.upper_node == r.deleted_node and b.lower_node->time < r.inserted_node->time) {
            Branch candidate_recombined_branch = Branch(b.lower_node, r.inserted_node);
            if (r.create(candidate_recombined_branch)) {
                source_candidates.push_back(b);
            }
        }
    }
    if (source_candidates.size() == 1) {
        r.source_branch = source_candidates[0];
        r.start_time = sample_start_time(r.source_branch, 20, r.inserted_node->time, cut_time);
    } else if (source_candidates.size() == 2) {
        pair<Branch, double> breakpoint = sample_start_time(source_candidates[0], source_candidates[1], 40, r.inserted_node->time, cut_time);
        r.source_branch = breakpoint.first;
        r.start_time = breakpoint.second;
    } else {
        cout << r.pos << " " << source_candidates.size() << endl;
        cerr << "no candidates in smc sampling" << endl;
        exit(1);
    }
    r.find_target_branch();
    r.find_recomb_info();
    assert(r.target_branch != Branch());
    assert(r.merging_branch != Branch());
    assert(r.start_time <= r.inserted_node->time);
    assert(r.start_time >= cut_time);
}

double RSP_smc::draw_start_time(double lo, double hi) {
    double lam_hi = lambda(hi);
    vector<double> piece_lo, piece_hi, piece_mass;
    int j = level_of(lo);
    double a = lo;
    while (a < hi) {
        double b = (j + 1 < (int) level_times.size()) ? min(level_times[j+1], hi) : hi;
        double k = level_rates[j];
        piece_lo.push_back(a);
        piece_hi.push_back(b);
        piece_mass.push_back((exp(lambda(b) - lam_hi) - exp(lambda(a) - lam_hi))/k);
        a = b;
        j += 1;
    }
    double w = uniform_random()*accumulate(piece_mass.begin(), piece_mass.end(), 0.0);
    int p = 0;
    while (p + 1 < (int) piece_mass.size() and w > piece_mass[p]) {
        w -= piece_mass[p];
        p += 1;
    }
    double k = level_rates[level_of(piece_lo[p])];
    double u = piece_lo[p] + log1p(uniform_random()*expm1(k*(piece_hi[p] - piece_lo[p])))/k;
    return min(max(u, piece_lo[p]), piece_hi[p]);
}

vector<Branch> RSP_smc::source_candidates(Recombination &r) {
    vector<Branch> candidates;
    for (Branch b : r.deleted_branches) {
        if (b.upper_node == r.deleted_node and b.lower_node->time < r.inserted_node->time) {
            if (r.create(Branch(b.lower_node, r.inserted_node))) {
                candidates.push_back(b);
            }
        }
    }
    return candidates;
}

double RSP_smc::start_lower_bound(const Branch &candidate, double cut_time, const Branch &own) {
    double lo = candidate.lower_node->time;
    return (candidate == own) ? max(cut_time, lo) : lo;
}

void RSP_smc::sample_recombination(Recombination &r, double cut_time, Tree &tree, const Branch &own) {
    if (r.pos == 0) {
        return;
    }
    if (r.start_time > 0) {
        return;
    }
    if (r.deleted_branches.size() == 0) {
        return;
    }
    set_tree(tree);
    double join_time = r.inserted_node->time;
    double lam_join = lambda(join_time);
    vector<Branch> candidates = source_candidates(r);
    vector<double> weights(candidates.size(), 0.0);
    for (int i = 0; i < (int) candidates.size(); i++) {
        double lo = start_lower_bound(candidates[i], cut_time, own);
        double hi = min(join_time, candidates[i].upper_node->time);
        if (lo < hi) {
            weights[i] = exp_lambda_integral(lo, hi)*exp(lambda(hi) - lam_join);
        }
    }
    double total = accumulate(weights.begin(), weights.end(), 0.0);
    if (total <= 0) {
        cout << r.pos << " " << candidates.size() << endl;
        cerr << "no candidates in smc sampling" << endl;
        exit(1);
    }
    double w = uniform_random()*total;
    int c = 0;
    while (c + 1 < (int) candidates.size() and w > weights[c]) {
        w -= weights[c];
        c += 1;
    }
    r.source_branch = candidates[c];
    r.start_time = draw_start_time(start_lower_bound(r.source_branch, cut_time, own), min(join_time, r.source_branch.upper_node->time));
    r.find_target_branch();
    r.find_recomb_info();
}

double RSP_smc::log_start_density(Recombination &r) {
    double join_time = r.inserted_node->time;
    double lo = r.source_branch.lower_node->time;
    double hi = min(join_time, r.source_branch.upper_node->time);
    if (r.start_time < lo or r.start_time > hi) {
        return -numeric_limits<double>::infinity();
    }
    return lambda(r.start_time) - lambda(join_time);
}

double RSP_smc::log_start_marginal(Recombination &r, double cut_time, const Branch &own) {
    double join_time = r.inserted_node->time;
    double lam_join = lambda(join_time);
    double total = 0;
    for (Branch &b : source_candidates(r)) {
        double lo = start_lower_bound(b, cut_time, own);
        double hi = min(join_time, b.upper_node->time);
        if (lo < hi) {
            total += exp_lambda_integral(lo, hi)*exp(lambda(hi) - lam_join);
        }
    }
    return log(total);
}

void RSP_smc::set_tree(const Flat_tree &tree) {
    node_times.clear();
    for (auto &x : tree.parents) {
        if (x.first->time > 0) {
            node_times.push_back(x.first->time);
        }
    }
    int n = (int) node_times.size();
    level_times.assign(1, 0.0);
    level_rates.assign(1, 1.0 + n);
    level_lambda.assign(1, 0.0);
    int i = 0;
    while (i < n) {
        double t = node_times[i];
        int e = i + 1;
        while (e < n and node_times[e] == t) {
            e++;
        }
        int j = (int) level_times.size();
        level_times.push_back(t);
        level_rates.push_back(1.0 + (n - e));
        level_lambda.push_back(level_lambda[j-1] + level_rates[j-1]*(t - level_times[j-1]));
        i = e;
    }
}

double RSP_smc::unchanged_recomb_length(const Flat_tree &tree) {
    int m = (int) level_times.size();
    int f = 0;
    int m_prev = (int) prev_times.size();
    while (f < m and f < m_prev and level_times[f] == prev_times[f] and level_rates[f] == prev_rates[f] and level_lambda[f] == prev_lambda[f]) {
        f++;
    }
    int j0 = max(f - 1, 0);
    a_sum.resize(m);
    e_sum.resize(m);
    eg_sum.resize(m);
    g_tail.resize(m + 1);
    e_lev.resize(m);
    g_lev.resize(m);
    a_sum[0] = 0;
    e_sum[0] = 0;
    eg_sum[0] = 0;
    e_lev[m-1] = 0;
    g_lev[m-1] = 0;
    g_tail[m-1] = 0;
    g_tail[m] = 0;
    for (int j = j0; j + 1 < m; j++) {
        double k = level_rates[j];
        double d = level_times[j+1] - level_times[j];
        double decay = -expm1(-k*d);
        e_lev[j] = exp(level_lambda[j+1])*decay/k;
        g_lev[j] = exp(-level_lambda[j])*decay/k;
        a_sum[j+1] = a_sum[j] + d/k + (-decay)/k/k;
        e_sum[j+1] = e_sum[j] + e_lev[j];
    }
    for (int j = m - 2; j >= 0; j--) {
        g_tail[j] = g_tail[j+1] + g_lev[j];
    }
    for (int j = 0; j + 1 < m; j++) {
        eg_sum[j+1] = eg_sum[j] + e_lev[j]*g_tail[j+1];
    }
    double total = 0;
    int jc = 0;
    for (auto &x : tree.parents) {
        while (jc + 1 < m and level_times[jc+1] <= x.first->time) {
            jc++;
        }
        if (x.second->index == -1) {
            continue;
        }
        int jp = level_of(x.second->time);
        total += (a_sum[jp] - a_sum[jc]) + (eg_sum[jp] - eg_sum[jc]) - g_tail[jp]*(e_sum[jp] - e_sum[jc]);
    }
    prev_times = level_times;
    prev_rates = level_rates;
    prev_lambda = level_lambda;
    return total;
}

void RSP_smc::approx_sample_recombination(Recombination &r, double cut_time) {
    if (r.pos == 0) {
        return;
    }
    if (r.start_time > 0) {
        return;
    }
    if (r.deleted_branches.size() == 0) {
        return;
    }
    vector<Branch> source_candidates;
    for (Branch b : r.deleted_branches) {
        if (b.upper_node == r.deleted_node and b.lower_node->time < r.inserted_node->time) {
            Branch candidate_recombined_branch = Branch(b.lower_node, r.inserted_node);
            if (r.create(candidate_recombined_branch)) {
                source_candidates.push_back(b);
            }
        }
    }
    double lb1, lb2, ub1, ub2;
    if (source_candidates.size() == 1) {
        r.source_branch = source_candidates[0];
        lb1 = max(cut_time, r.source_branch.lower_node->time);
        ub1 = min(r.source_branch.upper_node->time, r.inserted_node->time);
        // r.start_time = random_time(lb1, ub1, 0.5);
        r.start_time = choose_time(lb1, ub1);
    } else if (source_candidates.size() == 2) {
        lb1 = max(cut_time, source_candidates[0].lower_node->time);
        lb2 = max(cut_time, source_candidates[1].lower_node->time);
        ub1 = min(source_candidates[0].upper_node->time, r.inserted_node->time);
        ub2 = min(source_candidates[1].upper_node->time, r.inserted_node->time);
        double q = (ub1 - lb1)/(ub1 + ub2 - lb1 - lb2);
        double p = 0.5;
        if (p <= q) {
            r.source_branch = source_candidates[0];
            // r.start_time = random_time(lb1, ub1, 0.5);
            r.start_time = choose_time(lb1, ub1);
        } else {
            r.source_branch = source_candidates[1];
            // r.start_time = random_time(lb2, ub2, 0.5);
            r.start_time = choose_time(lb2, ub2);
        }
    } else {
        cout << r.pos << " " << source_candidates.size() << endl;
        cerr << "no candidates in smc sampling" << endl;
        exit(1);
    }
    if (r.deleted_node->time == r.inserted_node->time) {
        r.inserted_node->time = nextafter(r.inserted_node->time, numeric_limits<double>::infinity());
    }
    r.find_target_branch();
    r.find_recomb_info();
    assert(r.target_branch != Branch());
    assert(r.merging_branch != Branch());
    assert(r.start_time >= cut_time);
    double ub = min(r.deleted_node->time, r.inserted_node->time);
    if (r.start_time >= ub) {
        r.start_time = nextafter(ub, -numeric_limits<double>::infinity());
    }
    assert(r.start_time <= ub);
}

void RSP_smc::adjust(Recombination &r, double cut_time) {
    if (r.pos == 0) {
        return;
    }
    if (r.start_time > 0) {
        return;
    }
    if (r.deleted_branches.size() == 0) {
        return;
    }
    double lb, ub;
    lb = max(cut_time, r.source_branch.lower_node->time);
    ub = min(r.deleted_node->time, r.inserted_node->time);
    r.start_time = choose_time(lb, ub);
    if (r.start_time >= ub or r.start_time <= lb) {
        r.start_time = 0.5*(lb + ub);
    }
    assert(r.start_time <= ub);
}

void RSP_smc::approx_sample_recombination(Recombination &r, double cut_time, double n) {
    if (r.pos == 0) {
        return;
    }
    if (r.start_time > 0) {
        return;
    }
    if (r.deleted_branches.size() == 0) {
        return;
    }
    vector<Branch> source_candidates;
    for (Branch b : r.deleted_branches) {
        if (b.upper_node == r.deleted_node and b.lower_node->time < r.inserted_node->time) {
            Branch candidate_recombined_branch = Branch(b.lower_node, r.inserted_node);
            if (r.create(candidate_recombined_branch)) {
                source_candidates.push_back(b);
            }
        }
    }
    double lb1, lb2, ub1, ub2;
    if (source_candidates.size() == 1) {
        r.source_branch = source_candidates[0];
        lb1 = max(cut_time, r.source_branch.lower_node->time);
        ub1 = min(r.source_branch.upper_node->time, r.inserted_node->time);
        r.start_time = choose_time(lb1, ub1, n);
    } else if (source_candidates.size() == 2) {
        lb1 = max(cut_time, source_candidates[0].lower_node->time);
        lb2 = max(cut_time, source_candidates[1].lower_node->time);
        ub1 = min(source_candidates[0].upper_node->time, r.inserted_node->time);
        ub2 = min(source_candidates[1].upper_node->time, r.inserted_node->time);
        double q = (ub1 - lb1)/(ub1 + ub2 - lb1 - lb2);
        double p = 0.5;
        if (p <= q) {
            r.source_branch = source_candidates[0];
            r.start_time = choose_time(lb1, ub1, n);
        } else {
            r.source_branch = source_candidates[1];
            r.start_time = choose_time(lb2, ub2, n);
        }
    } else {
        cout << r.pos << " " << source_candidates.size() << endl;
        cerr << "no candidates in smc sampling" << endl;
        exit(1);
    }
    if (r.deleted_node->time == r.inserted_node->time) {
        r.inserted_node->time = nextafter(r.inserted_node->time, numeric_limits<double>::infinity());
    }
    r.find_target_branch();
    r.find_recomb_info();
    assert(r.target_branch != Branch());
    assert(r.merging_branch != Branch());
    assert(r.start_time >= cut_time);
    double ub = min(r.deleted_node->time, r.inserted_node->time);
    if (r.start_time >= ub) {
        r.start_time = nextafter(ub, -numeric_limits<double>::infinity());
    }
    assert(r.start_time <= ub);
}

void RSP_smc::adjust(Recombination &r, double cut_time, double n) {
    if (r.pos == 0) {
        return;
    }
    if (r.start_time > 0) {
        return;
    }
    if (r.deleted_branches.size() == 0) {
        return;
    }
    double lb, ub;
    lb = max(cut_time, r.source_branch.lower_node->time);
    ub = min(r.deleted_node->time, r.inserted_node->time);
    r.start_time = choose_time(lb, ub, n);
    if (r.start_time >= ub or r.start_time <= lb) {
        r.start_time = 0.5*(lb + ub);
    }
    assert(r.start_time <= ub);
}

// private methods:

double RSP_smc::recomb_pdf(double s, double t) {
    double pdf = 1.0;
    double curr_time = s;
    double next_coalescence_time = s;
    map<double, int>::iterator rate_it = coalescence_rates.upper_bound(s);
    rate_it--;
    int rate;
    while (next_coalescence_time < t) {
        rate = rate_it->second;
        rate_it++;
        next_coalescence_time = rate_it->first;
        pdf *= exp(-rate*(-curr_time + min(t, next_coalescence_time)));
        curr_time = next_coalescence_time;
    }
    return pdf;
}

void RSP_smc::get_coalescence_rate(Tree &tree, Recombination &r, double cut_time) {
    coalescence_rates.clear();
    vector<double> coalescence_times = {cut_time};
    for (auto &x : tree.parents) {
        if (x.first->time > cut_time and x.second != r.deleted_node) {
            coalescence_times.push_back(x.first->time);
        }
    }
    coalescence_times.push_back(numeric_limits<double>::infinity());
    sort(coalescence_times.begin(), coalescence_times.end());
    int n = (int) coalescence_times.size();
    for (int i = 0; i < n; i++) {
        coalescence_rates[coalescence_times[i]] = n-i-1;
    }
}

double RSP_smc::random_time(double lb, double ub) {
    double t = (uniform_random()*0.01 + 0.99)*(ub - lb) + lb;
    return t;
}

double RSP_smc::random_time(double lb, double ub, double q) {
    double t = (q*uniform_random() + (1 - q))*(ub - lb) + lb;
    return t;
}

double RSP_smc::choose_time(double lb, double ub) {
    assert(!isinf(ub));
    double l = exp(lb);
    double u = exp(ub);
    double m = 0.5*(l + u);
    double mt = log(m);
    if (ub - lb < 0.01) {
        mt = 0.5*(lb + ub);
    }
    assert(mt >= lb and mt <= ub);
    return mt;
}

double RSP_smc::choose_time(double lb, double ub, double n) {
    double lambda = 0;
    double lambda_l = 0;
    double lambda_u = 0;
    double mt = 0;
    if (n > 1) {
        lambda_l = n/(n + (1 - n)*exp(-0.5*lb));
        lambda_u = n/(n + (1 - n)*exp(-0.5*ub));
        lambda = lambda_l;
    } else {
        lambda = 1;
    }
    double l = exp(lambda*lb - lambda*ub);
    double u = 1;
    double m = 0.5*(l + u);
    mt = ub + log(m)/lambda;
    if (ub - lb < 0.01) {
        mt = 0.5*(lb + ub);
    }
    assert(mt >= lb and mt <= ub);
    return mt;
}

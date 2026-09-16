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
    return (int) (upper_bound(level_times.begin(), level_times.end(), s) - level_times.begin()) - 1;
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

double RSP_smc::unchanged_recomb_length(Tree &tree) {
    int m = (int) level_times.size();
    vector<double> a_sum(m, 0.0), e_sum(m, 0.0), eg_sum(m, 0.0), g_tail(m + 1, 0.0);
    vector<double> e_lev(m, 0.0), g_lev(m, 0.0);
    for (int j = 0; j + 1 < m; j++) {
        double k = level_rates[j];
        double d = level_times[j+1] - level_times[j];
        double decay = -expm1(-k*d);
        e_lev[j] = exp(level_lambda[j+1])*decay/k;
        g_lev[j] = exp(-level_lambda[j])*decay/k;
    }
    for (int j = m - 2; j >= 0; j--) {
        g_tail[j] = g_tail[j+1] + g_lev[j];
    }
    for (int j = 0; j + 1 < m; j++) {
        double k = level_rates[j];
        double d = level_times[j+1] - level_times[j];
        a_sum[j+1] = a_sum[j] + d/k + expm1(-k*d)/k/k;
        e_sum[j+1] = e_sum[j] + e_lev[j];
        eg_sum[j+1] = eg_sum[j] + e_lev[j]*g_tail[j+1];
    }
    double total = 0;
    for (auto &x : tree.parents) {
        if (x.second->index == -1) {
            continue;
        }
        int jc = level_of(x.first->time);
        int jp = level_of(x.second->time);
        total += (a_sum[jp] - a_sum[jc]) + (eg_sum[jp] - eg_sum[jc]) - g_tail[jp]*(e_sum[jp] - e_sum[jc]);
    }
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

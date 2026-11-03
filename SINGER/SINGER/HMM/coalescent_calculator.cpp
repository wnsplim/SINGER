//
//  coalescent_calculator.cpp
//  SINGER
//
//  Created by Yun Deng on 6/25/23.
//

#include "coalescent_calculator.hpp"

coalescent_calculator::coalescent_calculator(double x) {
    cut_time = x;
    t.push_back(cut_time);
    sgn.push_back(0);
}

coalescent_calculator::~coalescent_calculator() {}

void coalescent_calculator::set_events(vector<pair<double, int>> &events) {
    sort(events.begin(), events.end());
    t.assign(1, cut_time);
    sgn.assign(1, 0);
    for (auto &e : events) {
        t.push_back(e.first);
        sgn.push_back(e.second);
    }
    rebuild_from = 0;
    memo_n = memo_next = 0;
}

void coalescent_calculator::start(set<Branch> &branches) {
    vector<pair<double, int>> events;
    for (const Branch &b : branches) {
        if (b.lower_node->time > cut_time) {
            events.push_back({b.lower_node->time, b.lower_node->is_sample ? 1 : -1});
        }
    }
    set_events(events);
}

void coalescent_calculator::start(Tree &tree) {
    vector<pair<double, int>> events;
    for (auto &x : tree.parents) {
        if (x.second->time > cut_time and x.first->time > cut_time) {
            events.push_back({x.first->time, x.first->is_sample ? 1 : -1});
        }
    }
    set_events(events);
}

void coalescent_calculator::update(Recombination &r) {
    update(r.deleted_node->time, r.inserted_node->time);
}

void coalescent_calculator::update(double t_old, double t_new) {
    memo_n = memo_next = 0;
    if (t_old > cut_time) {
        auto it = lower_bound(t.begin(), t.end(), t_old);
        if (it != t.end() and *it == t_old) {
            rebuild_from = min(rebuild_from, (int) (it - t.begin()) - 1);
            sgn.erase(sgn.begin() + (it - t.begin()));
            t.erase(it);
        }
    }
    if (t_new > cut_time) {
        auto it = lower_bound(t.begin(), t.end(), t_new);
        rebuild_from = min(rebuild_from, (int) (it - t.begin()) - 1);
        sgn.insert(sgn.begin() + (it - t.begin()), -1);
        t.insert(it, t_new);
    }
}

void coalescent_calculator::refresh() {
    int m = (int) t.size();
    if (rebuild_from >= m - 1 and (int) Lam.size() == m) {
        return;
    }
    int j0 = max(0, rebuild_from);
    if ((int) Lam.size() != m) {
        Lam.resize(m);
        E.resize(m);
        T.resize(m);
        W.resize(m);
        B.resize(m);
        kk.resize(m);
        j0 = 0;
    }
    if (j0 == 0) {
        E[0] = exp(-Lam[0]);
        kk[0] = 1 - accumulate(sgn.begin(), sgn.end(), 0);
    }
    for (int j = j0; j + 1 < m; j++) {
        kk[j+1] = kk[j] + sgn[j+1];
        double k = kk[j] + extra;
        double dt = t[j+1] - t[j];
        double ea = E[j];
        if (k == 0) {
            E[j+1] = ea;
            Lam[j+1] = Lam[j];
            B[j+1] = B[j] + W[j]*dt + 0.5*dt*dt;
            W[j+1] = W[j] + dt;
            continue;
        }
        E[j+1] = exp(-(Lam[j] + k*dt));
        Lam[j+1] = Lam[j] + k*dt;
        double e = -expm1(-k*dt);
        B[j+1] = B[j] + W[j]*e/k + dt/k - e/k/k;
        W[j+1] = W[j]*(1 - e) + e/k;
    }
    T[m-1] = E[m-1]/(kk[m-1] + extra);
    for (int j = m - 2; j >= 0; j--) {
        double k = kk[j] + extra;
        T[j] = T[j+1] + (k == 0 ? 0 : E[j]*(-expm1(-k*(t[j+1] - t[j])))/k);
    }
    rebuild_from = m;
}

double coalescent_calculator::tail_at(double x) {
    if (isinf(x)) {
        return 0;
    }
    int m = (int) t.size();
    int j = (int) (upper_bound(t.begin(), t.end(), x) - t.begin()) - 1;
    double k = kk[j] + extra;
    double above = j + 1 < m ? T[j+1] : 0;
    if (k == 0) {
        return above;
    }
    double ex = exp(-(Lam[j] + k*(x - t[j])));
    if (j + 1 == m) {
        return ex/k;
    }
    return above + ex*(-expm1(-k*(t[j+1] - x)))/k;
}

double coalescent_calculator::tail(double x) {
    refresh();
    for (int k = 0; k < memo_n; k++) {
        if (memo_x[k] == x) {
            return memo_t[k];
        }
    }
    double v = tail_at(x);
    memo_x[memo_next] = x;
    memo_t[memo_next] = v;
    memo_next = (memo_next + 1) % 4;
    memo_n = min(memo_n + 1, 4);
    return v;
}

double coalescent_calculator::prob(double x, double y) {
    return max(tail(x) - tail(y), 0.0);
}

double coalescent_calculator::rel_prob(double x, double y) {
    refresh();
    int m = (int) t.size();
    int j = (int) (upper_bound(t.begin(), t.end(), x) - t.begin()) - 1;
    double p = 0, dl = 0, a = x;
    for (; a < y and j < m; j++) {
        double k = kk[j] + extra;
        double b = (j + 1 < m) ? min(t[j+1], y) : y;
        if (k > 0) {
            p += exp(-dl)*(-expm1(-k*(b - a)))/k;
            dl += k*(b - a);
        }
        a = b;
    }
    return p;
}

double coalescent_calculator::prob_inv(double x, double y, double u) {
    refresh();
    int m = (int) t.size();
    int j = (int) (upper_bound(t.begin(), t.end(), x) - t.begin()) - 1;
    double k0 = kk[j] + extra;
    if (k0 > 0 and (j + 1 >= m or y <= t[j+1])) {
        double target = u*((-expm1(-k0*(y - x)))/k0);
        return min(y, x - log1p(-min(1.0, target*k0))/k0);
    }
    double target = u*rel_prob(x, y);
    double dl = 0, a = x;
    for (; j < m; j++) {
        double k = kk[j] + extra;
        double b = (j + 1 < m) ? min(t[j+1], y) : y;
        if (k > 0) {
            double ea = exp(-dl);
            double seg = ea*(-expm1(-k*(b - a)))/k;
            if (target <= seg or b >= y) {
                return min(y, a - log1p(-min(1.0, target*k/ea))/k);
            }
            target -= seg;
            dl += k*(b - a);
        }
        a = b;
    }
    return y;
}

double coalescent_calculator::log_density(double x, double y, double s) {
    double p = rel_prob(x, y);
    int jx = (int) (upper_bound(t.begin(), t.end(), x) - t.begin()) - 1;
    int js = (int) (upper_bound(t.begin(), t.end(), s) - t.begin()) - 1;
    double lx = Lam[jx] + (kk[jx] + extra)*(x - t[jx]);
    double ls = Lam[js] + (kk[js] + extra)*(s - t[js]);
    return p > 0 ? lx - ls - log(p) : -numeric_limits<double>::infinity();
}

double coalescent_calculator::surv(double x) {
    if (isinf(x)) {
        return 0;
    }
    refresh();
    int j = (int) (upper_bound(t.begin(), t.end(), x) - t.begin()) - 1;
    return exp(-(Lam[j] + (kk[j] + extra)*(x - t[j])));
}

double coalescent_calculator::surv_inv(double p) {
    if (p <= 0) {
        return numeric_limits<double>::infinity();
    }
    refresh();
    int m = (int) t.size();
    double l = -log(p);
    int j = (int) (upper_bound(Lam.begin(), Lam.begin() + m, l) - Lam.begin()) - 1;
    j = max(j, 0);
    return t[j] + (l - Lam[j])/(kk[j] + extra);
}

double coalescent_calculator::recomb_mass(double s, double v) {
    refresh();
    double x = min(s, v);
    int jx = (int) (upper_bound(t.begin(), t.end(), x) - t.begin()) - 1;
    double kx = kk[jx] + extra;
    double dx = x - t[jx];
    double ex = -expm1(-kx*dx);
    double b2 = B[jx] + W[jx]*ex/kx + dx/kx - ex/kx/kx;
    if (s > v) {
        return b2;
    }
    double wx = W[jx]*(1 - ex) + ex/kx;
    double lx = Lam[jx] + kx*dx;
    return b2 + wx*(tail(x) - tail(v))*exp(lx);
}

double coalescent_calculator::find_median(double x, double y) {
    if (x == y) {
        return x;
    }
    if (y - x < 0.001) {
        return 0.5*(x + y);
    }
    refresh();
    int m = (int) t.size();
    int j = (int) (upper_bound(t.begin(), t.end(), x) - t.begin()) - 1;
    double dl = 0; // referenced to x so the exponential factors out and cannot underflow
    double P = 0, Q1 = 0;
    double a = x;
    while (a < y and j < m) {
        double k = kk[j] + extra;
        double b = (j + 1 < m) ? min(t[j+1], y) : y;
        double ea = exp(-dl);
        if (isinf(b)) {
            P += ea/k;
            Q1 += ((a - x)*ea)/k + ea/k/k;
            break;
        }
        double eb = exp(-(dl + k*(b - a)));
        P += (ea - eb)/k;
        Q1 += ((a - x)*ea - (b - x)*eb)/k + (ea - eb)/k/k;
        dl += k*(b - a);
        a = b;
        j += 1;
    }
    return x + Q1/P;
}

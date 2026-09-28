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
        G.resize(m);
        Q.resize(m);
        W.resize(m);
        B.resize(m);
        kk.resize(m);
        j0 = 0; // k = m - j shifts for the whole array when the size changes
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
            G[j+1] = G[j];
            Q[j+1] = Q[j];
            B[j+1] = B[j] + W[j]*dt + 0.5*dt*dt;
            W[j+1] = W[j] + dt;
            continue;
        }
        double eb = exp(-(Lam[j] + k*dt));
        E[j+1] = eb;
        Lam[j+1] = Lam[j] + k*dt;
        G[j+1] = G[j] + (ea - eb)/k;
        Q[j+1] = Q[j] + ((t[j] - cut_time)*ea - (t[j+1] - cut_time)*eb)/k + (ea - eb)/k/k;
        double e = -expm1(-k*dt);
        B[j+1] = B[j] + W[j]*e/k + dt/k - e/k/k;
        W[j+1] = W[j]*(1 - e) + e/k;
    }
    double ez = exp(-Lam[m-1]);
    double kz = 1 + extra;
    tail_G = ez/kz;
    tail_Q = ((t[m-1] - cut_time)/kz + 1/(kz*kz))*ez;
    first_moment = G[m-1] + tail_G;
    rebuild_from = m;
}

void coalescent_calculator::at(double x, double &g, double &q) {
    int m = (int) t.size();
    if (isinf(x)) {
        g = G[m-1] + tail_G;
        q = Q[m-1] + tail_Q;
        return;
    }
    int j = (int) (upper_bound(t.begin(), t.end(), x) - t.begin()) - 1;
    double k = kk[j] + extra;
    if (k == 0) {
        g = G[j];
        q = Q[j];
        return;
    }
    double ea = E[j];
    double ex = exp(-(Lam[j] + k*(x - t[j])));
    g = G[j] + (ea - ex)/k;
    q = Q[j] + ((t[j] - cut_time)*ea - (x - cut_time)*ex)/k + (ea - ex)/k/k;
}

double coalescent_calculator::prob(double x, double y) {
    refresh();
    double gx, qx, gy, qy;
    at(x, gx, qx);
    at(y, gy, qy);
    return max(gy - gx, 0.0);
}

double coalescent_calculator::find_median(double x, double y) {
    return compute_time_weights(x, y).first;
}

double coalescent_calculator::surv(double x) {
    if (isinf(x)) {
        return 0;
    }
    refresh();
    int j = (int) (upper_bound(t.begin(), t.end(), x) - t.begin()) - 1;
    return exp(-(Lam[j] + (kk[j] + extra)*(x - t[j])));
}

double coalescent_calculator::rate(double x) {
    refresh();
    int j = (int) (upper_bound(t.begin(), t.end(), x) - t.begin()) - 1;
    return (double)(kk[j] + extra);
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
    int m = (int) t.size();
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
    double gx = G[jx] + (E[jx] - exp(-lx))/kx;
    double gv;
    if (isinf(v)) {
        gv = G[m-1] + tail_G;
    } else {
        int jv = (int) (upper_bound(t.begin(), t.end(), v) - t.begin()) - 1;
        double kv = kk[jv] + extra;
        double lv = Lam[jv] + kv*(v - t[jv]);
        gv = G[jv] + (E[jv] - exp(-lv))/kv;
    }
    return b2 + wx*(gv - gx)*exp(lx);
}

pair<double, double> coalescent_calculator::compute_time_weights(double x, double y) {
    if (x == y) {
        return {x, 0};
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
    double time = x + Q1/P;
    double w = E[j < m ? j : m-1]*((x - cut_time)*P + Q1)/first_moment;
    if (y - x < 0.001) {
        time = 0.5*(x + y);
        w = (time - cut_time)*prob(x, y)/first_moment;
    }
    return {time, w};
}

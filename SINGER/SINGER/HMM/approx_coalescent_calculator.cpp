//
//  approx_coalescent_calculator.cpp
//  SINGER
//
//  Created by Yun Deng on 6/25/23.
//

#include "approx_coalescent_calculator.hpp"

approx_coalescent_calculator::approx_coalescent_calculator(double x) {
    cut_time = x;
    t.push_back(cut_time);
}

approx_coalescent_calculator::~approx_coalescent_calculator() {}

void approx_coalescent_calculator::start(set<Branch> &branches) {
    for (const Branch &b : branches) {
        if (b.lower_node->time > cut_time) {
            t.push_back(b.lower_node->time);
        }
    }
    sort(t.begin(), t.end());
    rebuild_from = 0;
}

void approx_coalescent_calculator::start(Tree &tree) {
    for (auto &x : tree.parents) {
        if (x.second->time > cut_time and x.first->time > cut_time) {
            t.push_back(x.first->time);
        }
    }
    sort(t.begin(), t.end());
    rebuild_from = 0;
}

void approx_coalescent_calculator::update(Recombination &r) {
    double t_old = r.deleted_node->time;
    double t_new = r.inserted_node->time;
    if (t_old > cut_time) {
        auto it = lower_bound(t.begin(), t.end(), t_old);
        if (it != t.end() and *it == t_old) {
            rebuild_from = min(rebuild_from, (int) (it - t.begin()) - 1);
            t.erase(it);
        }
    }
    if (t_new > cut_time) {
        auto it = lower_bound(t.begin(), t.end(), t_new);
        rebuild_from = min(rebuild_from, (int) (it - t.begin()) - 1);
        t.insert(it, t_new);
    }
}

void approx_coalescent_calculator::refresh() {
    int m = (int) t.size();
    if (rebuild_from >= m - 1 and (int) Lam.size() == m) {
        return;
    }
    int j0 = max(0, rebuild_from);
    if ((int) Lam.size() != m) {
        Lam.resize(m);
        G.resize(m);
        Q.resize(m);
        j0 = 0; // k = m - j shifts for the whole array when the size changes
    }
    for (int j = j0; j + 1 < m; j++) {
        double k = m - j;
        double ea = exp(-Lam[j]);
        double eb = exp(-(Lam[j] + k*(t[j+1] - t[j])));
        Lam[j+1] = Lam[j] + k*(t[j+1] - t[j]);
        G[j+1] = G[j] + (ea - eb)/k;
        Q[j+1] = Q[j] + ((t[j] - cut_time)*ea - (t[j+1] - cut_time)*eb)/k + (ea - eb)/k/k;
    }
    double ez = exp(-Lam[m-1]);
    tail_G = ez;
    tail_Q = (t[m-1] - cut_time)*ez + ez;
    first_moment = G[m-1] + tail_G;
    rebuild_from = m;
}

void approx_coalescent_calculator::compute_first_moment() {}

void approx_coalescent_calculator::at(double x, double &g, double &q) {
    int m = (int) t.size();
    if (isinf(x)) {
        g = G[m-1] + tail_G;
        q = Q[m-1] + tail_Q;
        return;
    }
    int j = (int) (upper_bound(t.begin(), t.end(), x) - t.begin()) - 1;
    double k = m - j;
    double ea = exp(-Lam[j]);
    double ex = exp(-(Lam[j] + k*(x - t[j])));
    g = G[j] + (ea - ex)/k;
    q = Q[j] + ((t[j] - cut_time)*ea - (x - cut_time)*ex)/k + (ea - ex)/k/k;
}

double approx_coalescent_calculator::prob(double x, double y) {
    refresh();
    double gx, qx, gy, qy;
    at(x, gx, qx);
    at(y, gy, qy);
    return max(gy - gx, 0.0);
}

double approx_coalescent_calculator::find_median(double x, double y) {
    return compute_time_weights(x, y).first;
}

pair<double, double> approx_coalescent_calculator::compute_time_weights(double x, double y) {
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
        double k = m - j;
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
    double w = exp(-Lam[j < m ? j : m-1])*((x - cut_time)*P + Q1)/first_moment;
    if (y - x < 0.001) {
        time = 0.5*(x + y);
        w = (time - cut_time)*prob(x, y)/first_moment;
    }
    return {time, w};
}

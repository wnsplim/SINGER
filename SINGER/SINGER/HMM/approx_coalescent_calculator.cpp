#include "approx_coalescent_calculator.hpp"

approx_coalescent_calculator::approx_coalescent_calculator(double t) {
    cut_time = t;
}

void approx_coalescent_calculator::start(set<Branch> &branches) {
    for (auto &x : branches) {
        if (x.upper_node->time > cut_time and x.lower_node->time <= cut_time) {
            n0 += 1;
        }
    }
}

void approx_coalescent_calculator::update(Recombination &r) {
    if (r.deleted_node->time > cut_time) {
        n0 -= 1;
    }
    if (r.inserted_node->time > cut_time) {
        n0 += 1;
    }
}

double approx_coalescent_calculator::prob(double x, double y) {
    double p1 = prob_integral(x);
    double p2 = prob_integral(y);
    return max(p2 - p1, 0.0);
}

double approx_coalescent_calculator::prob_integral(double x) {
    x -= cut_time;
    if (n0 == 1) {
        return 1 - exp(-x);
    }
    double v = n0 + (1 - n0)*exp(-0.5*x);
    double p = -2*log(v) - 2*n0/v;
    return p/(1 - n0)/(1 - n0);
}

double approx_coalescent_calculator::find_median(double x, double y) {
    if (x == y) {
        return x;
    }
    if (y - x <= 0.01) {
        return 0.5*(x + y);
    }
    if (n0 == 1) {
        double q = 0.5*(exp(-x) + exp(-y));
        return -log(q);
    }
    x -= cut_time;
    y -= cut_time;
    double nx = n0 + (1 - n0)*exp(-0.5*x);
    double ny = n0 + (1 - n0)*exp(-0.5*y);
    double nm = sqrt(nx*ny);
    double z = 2*log(n0 - 1) - 2*log(n0 - nm);
    if (isinf(y)) {
        z = x + 1;
    }
    if (z < x or z > y) {
        z = 0.5*(x + y);
    }
    return z + cut_time;
}

//
//  Interval.cpp
//  SINGER
//
//  Created by Yun Deng on 4/9/22.
//

#include "Interval.hpp"

Interval::Interval(Branch b, double tl, double tu, int init_pos) {
    branch = b;
    lb = tl;
    ub = tu;
    assert(lb <= ub);
    assert(b.lower_node->time <= tl and b.upper_node->time >= tu);
    start_pos = init_pos;
}

void Interval::fill_time() {
    if (ub == numeric_limits<double>::infinity()) {
        time = lb + log(2);
    } else if (abs(lb - ub) < 1e-3) {
        time = 0.5*(lb + ub);
    } else {
        double lq = 1 - exp(-lb);
        double uq = 1 - exp(-ub);
        if (uq - lq < 1e-3) {
            time = 0.5*(lb + ub);
        } else {
            double q = 0.5*(lq + uq);
            time = -log(1 - q);
        }
    }
    assert(!isinf(time));
    assert(time >= lb and time <= ub);
}

bool Interval::full(double t) {
    assert(lb >= t);
    return lb == max(t, branch.lower_node->time) and ub == branch.upper_node->time;
}

Interval_info::Interval_info() {
}

Interval_info::Interval_info(Branch b, double tl, double tu) {
    assert(tl <= tu);
    assert(tl >= b.lower_node->time and tu <= b.upper_node->time);
    branch = b;
    lb = tl;
    ub = tu;
}

bool Interval_info::operator<(const Interval_info& other) const {
    if (seed_pos != other.seed_pos) {
        return seed_pos < other.seed_pos;
    }
    if (branch != other.branch) {
        return branch < other.branch;
    }
    if (ub != other.ub) {
        return ub < other.ub;
    }
    return lb < other.lb;
}

//
//  Emission.cpp
//  SINGER
//

#include "Emission.hpp"

static double join_value(const Site_messages &s, Branch &branch, double time, Node *node) {
    int k = branch.lower_node->flat_index;
    double l2 = s.unit*(time - branch.lower_node->time);
    double l3 = s.unit*(time - node->time);
    const Msg &md = s.in[k];
    double d0 = md[0] + l2*md[1], d1 = md[1] + l2*md[0];
    double m0 = d0*(s.nu[0] + l3*s.nu[1]), m1 = d1*(s.nu[1] + l3*s.nu[0]);
    const Msg &a = s.above[k];
    double w;
    if (branch.upper_node->index == -1) {
        w = a[0]*m0 + a[1]*m1;
    } else {
        double l1 = s.unit*(branch.upper_node->time - time);
        w = a[0]*(m0 + l1*m1) + a[1]*(m1 + l1*m0);
    }
    return w - s.empty_term;
}

static double site_ratio_two(const Site_messages &s, Branch &branch, double time, Node *node) {
    double v = join_value(s, branch, time, node);
    for (const Site_messages &x : s.alt) {
        v += join_value(x, branch, time, node);
    }
    return v/s.scale;
}

double Emission::mut_emit(Branch &branch, double time, Node *node, double emit) {
    for (const Site_messages *site : sites) {
        emit *= site_ratio_two(*site, branch, time, node);
    }
    return emit;
}

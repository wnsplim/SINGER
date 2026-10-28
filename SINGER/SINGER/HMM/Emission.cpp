//
//  Emission.cpp
//  SINGER
//

#include "Emission.hpp"

double Emission::mut_emit(Branch &branch, double time, Node *node, double norm) {
    int k = branch.lower_node->flat_index;
    bool at_root = branch.upper_node->index == -1;
    double emit = norm;
    for (const Site_messages *site : sites) {
        const Site_messages &s = *site;
        double l2 = s.unit*(time - branch.lower_node->time);
        double l3 = s.unit*(time - node->time);
        const array<double, 2> &md = s.in[k];
        double d0 = md[0] + l2*md[1], d1 = md[1] + l2*md[0];
        double m0 = d0*(s.nu[0] + l3*s.nu[1]), m1 = d1*(s.nu[1] + l3*s.nu[0]);
        const array<double, 2> &a = s.above[k];
        double w;
        if (at_root) {
            w = a[0]*m0 + a[1]*m1;
        } else {
            double l1 = s.unit*(branch.upper_node->time - time);
            w = a[0]*(m0 + l1*m1) + a[1]*(m1 + l1*m0);
        }
        emit *= (w - s.empty_term)/s.scale;
    }
    return emit;
}

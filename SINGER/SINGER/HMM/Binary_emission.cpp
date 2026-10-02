//
//  Binary_emission.cpp
//  SINGER
//
//  Created by Yun Deng on 4/6/23.
//

#include "Binary_emission.hpp"

Binary_emission::Binary_emission() {}

Binary_emission::~Binary_emission() {}

double Binary_emission::null_emit(Branch &branch, double time, double theta, Node *node) {
    double ll = time - branch.lower_node->time;
    double lu = branch.upper_node->time - time;
    double l0 = time - node->time;
    return norm_ratio(ll, lu, l0);
}

double Binary_emission::mut_emit(Branch &branch, double time, double theta, double bin_size, vector<double> &mut_set, Node *node) {
    double emit_prob = 1;
    double old_prob = 1;
    double ll = time - branch.lower_node->time;
    double lu = branch.upper_node->time - time;
    double l0 = time - node->time;
    double u0 = ll*theta/bin_size;
    double u1 = isinf(lu) ? 1.0 : lu*theta/bin_size;
    double u2 = l0*theta/bin_size;
    double uo = isinf(lu) ? 1.0 : (ll + lu)*theta/bin_size;
    double pen[4] = {1.0, penalty, penalty*penalty, penalty*penalty*penalty};
    bool at_root = branch.upper_node->index == -1;
    double w0 = at_root ? ancestral_prob : 1.0;
    double w1 = at_root ? 1 - ancestral_prob : 1.0;
    bool lower_any = any_missing and branch.lower_node->missing_sites.size() > 0;
    bool query_any = any_missing and node->missing_sites.size() > 0;
    for (double m : mut_set) {
        int sl = (int) branch.lower_node->get_state(m);
        int su = (int) branch.upper_node->get_state(m);
        int s0 = (int) node->get_state(m);
        bool ml = lower_any and branch.lower_node->is_missing(m);
        bool m0 = query_any and node->is_missing(m);
        if (ml or m0) {
            double p0 = u0*penalty, p1 = at_root ? 1.0 : u1*penalty, p2 = u2*penalty, po = at_root ? 1.0 : uo*penalty;
            double num = 0, den = 0;
            for (int a = ml ? 0 : sl; a <= (ml ? 1 : sl); a++) {
                den += (a ? w1 : w0)*((a == su) ? 1.0 : po);
                for (int b = m0 ? 0 : s0; b <= (m0 ? 1 : s0); b++) {
                    num += w0*((a == 0) ? 1.0 : p0)*((su == 0) ? 1.0 : p1)*((b == 0) ? 1.0 : p2)
                         + w1*((a == 1) ? 1.0 : p0)*((su == 1) ? 1.0 : p1)*((b == 1) ? 1.0 : p2);
                }
            }
            emit_prob *= num;
            old_prob *= den;
            continue;
        }
        int base = at_root ? 0 : abs(sl - su);
        int k0 = sl + s0 + (at_root ? 0 : su);
        double t0 = (sl ? u0 : 1.0)*(su ? u1 : 1.0)*(s0 ? u2 : 1.0)*pen[k0 - base];
        double t1 = (sl ? 1.0 : u0)*(su ? 1.0 : u1)*(s0 ? 1.0 : u2)*pen[(at_root ? 2 : 3) - k0 - base];
        if (at_root) {
            emit_prob *= w0*t0 + w1*t1;
            old_prob *= sl ? w1 : w0;
        } else {
            emit_prob *= t0 + t1;
        }
        if (base) {
            old_prob *= uo;
        }
    }
    emit_prob *= norm_ratio(ll, lu, l0);
    emit_prob /= old_prob;
    assert(emit_prob != 0);
    return emit_prob;
}

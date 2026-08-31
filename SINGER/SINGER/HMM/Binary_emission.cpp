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
    for (double m : mut_set) {
        int sl = (int) branch.lower_node->get_state(m);
        int su = (int) branch.upper_node->get_state(m);
        int s0 = (int) node->get_state(m);
        int base = abs(sl - su);
        int k0 = sl + su + s0;
        double t0 = (sl ? u0 : 1.0)*(su ? u1 : 1.0)*(s0 ? u2 : 1.0)*pen[k0 - base];
        double t1 = (sl ? 1.0 : u0)*(su ? 1.0 : u1)*(s0 ? 1.0 : u2)*pen[3 - k0 - base];
        emit_prob *= t0 + t1;
        if (base) {
            old_prob *= uo;
        }
    }
    emit_prob *= norm_ratio(ll, lu, l0);
    emit_prob /= old_prob;
    assert(emit_prob != 0);
    return emit_prob;
}

double Binary_emission::emit(Branch &branch, double time, double theta, double bin_size, vector<double> &emissions, Node *node) {
    double emit_prob = 1;
    double old_prob = 1;
    double ll = time - branch.lower_node->time;
    double lu = branch.upper_node->time - time;
    double l0 = time - node->time;
    emit_prob = calculate_prob(theta, bin_size, ll, lu, l0, emissions[0], emissions[1], emissions[2]);
    old_prob = calculate_prob(theta*(ll + lu), bin_size, emissions[3]);
    emit_prob /= old_prob;
    // assert(emit_prob != 0);
    return emit_prob;
}

double Binary_emission::calculate_prob(double theta, double bin_size, double ll, double lu, double l0, int sl, int su, int s0) {
    double prob = 1;
    prob *= calculate_prob(ll*theta, bin_size, sl);
    if (!isinf(lu)) {
        prob *= calculate_prob(lu*theta, bin_size, su);
    }
    prob *= calculate_prob(l0*theta, bin_size, s0);
    return prob;
}

double Binary_emission::calculate_prob(double theta, double bin_size, int s) {
    if (isinf(theta)) {
        return 1.0;
    }
    double unit_theta = theta/bin_size;
    if (s == 0) {
        return exp(-theta);
    }
    if (s == 1) {
        return exp(-theta)*unit_theta;
    }
    return exp(-theta)*pow(unit_theta, s);
}

void Binary_emission::get_diff(vector<double> &mut_set, const Branch &branch, Node *node) {
    double sl = 0;
    double su = 0;
    double s0 = 0;
    double sm = 0;
    fill(diff.begin(), diff.end(), 0);
    for (double x : mut_set) {
        sl = branch.lower_node->get_state(x);
        su = branch.upper_node->get_state(x);
        s0 = node->get_state(x);
        if (sl + su + s0 > 1.5) {
            sm = 1;
        } else {
            sm = 0;
        }
        diff[0] += abs(sm - sl);
        diff[1] += abs(sm - su);
        diff[2] += abs(sm - s0);
        diff[3] += abs(sl - su);
    }
}

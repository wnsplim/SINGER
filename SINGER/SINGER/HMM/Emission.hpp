//
//  Emission.hpp
//  SINGER
//
//  Created by Yun Deng on 4/3/23.
//

#ifndef Emission_hpp
#define Emission_hpp

#include <stdio.h>
#include "Branch.hpp"
#include "ARG.hpp"

class Emission {

public:

    double crho = 0;
    double norm_scale = 1.0;
    double bin_width = 1.0;
    bool any_missing = false;

    virtual double null_emit(Branch &branch, double time, double theta, Node *node) = 0;
    virtual double mut_emit(Branch &branch, double time, double theta, double bin_size, vector<double> &mut_set, Node *node) = 0;
    virtual double emit(Branch &branch, double time, double theta, double bin_size, vector<double> &emissions, Node *node) = 0;

    void set_tree_product(double c, double rho, double p, double w) {
        crho = c*rho;
        norm_scale = p/(c + p - 1);
        bin_width = w;
    }

    double norm_ratio(double ll, double lu, double l0) {
        double a = crho*ll;
        double d = crho*l0;
        double km1;
        if (isinf(lu)) {
            km1 = a + d + a*d;
        } else {
            double b = crho*lu;
            double s = a + b;
            km1 = d + a*b*(1 + d)*(1 - s*(1 - s*(1 - s)));
        }
        double x = norm_scale*km1;
        return exp(-bin_width*x*(1 - x*(0.5 - x*((1.0/3.0) - x*0.25))));
    }
};

inline double branch_product(Tree &tree, double c, double rho) {
    double p = 1.0;
    for (auto &x : tree.parents) {
        double l = x.second->time - x.first->time;
        if (!isinf(l)) {
            p *= 1 + c*rho*l;
        }
    }
    return p;
}

inline double called_width(Node *query_node, ARG &a, int i) {
    double w = a.coordinates[i + 1] - a.coordinates[i];
    if (!a.any_missing or query_node == nullptr or query_node->missing_sites.size() == 0) {
        return w;
    }
    vector<double> &ms = query_node->missing_sites;
    auto lo = lower_bound(ms.begin(), ms.end(), a.coordinates[i]);
    auto hi = lower_bound(lo, ms.end(), a.coordinates[i + 1]);
    return w - (double) (hi - lo);
}

inline bool varying_rate(ARG &a, int lo, int hi) {
    double rho = a.thetas[lo]/(a.coordinates[lo + 1] - a.coordinates[lo]);
    for (int i = lo + 1; i < hi; i++) {
        if (a.thetas[i] != rho*(a.coordinates[i + 1] - a.coordinates[i])) {
            return true;
        }
    }
    return false;
}

inline double update_branch_product(double p, Recombination &r, double crho) {
    for (const Branch &b : r.deleted_branches) {
        double l = b.upper_node->time - b.lower_node->time;
        if (!isinf(l)) {
            p /= 1 + crho*l;
        }
    }
    for (const Branch &b : r.inserted_branches) {
        double l = b.upper_node->time - b.lower_node->time;
        if (!isinf(l)) {
            p *= 1 + crho*l;
        }
    }
    return p;
}

#endif /* Emission_hpp */

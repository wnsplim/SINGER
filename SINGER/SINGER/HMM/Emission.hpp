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

    virtual double mut_emit(Branch &branch, double time, double theta, double bin_size, vector<pair<double, int>> &mut_set, Node *node, double norm) = 0;

    void set_tree_product(double c, double rho, double p, double w) {
        crho = c*rho;
        norm_scale = p/(c + p - 1);
        bin_width = w;
    }

    double null_log(Branch &branch, double time, Node *node) {
        double ll = time - branch.lower_node->time;
        double lu = branch.upper_node->time - time;
        double l0 = time - node->time;
        double a = crho*ll;
        double d = crho*l0;
        double km1;
        if (isinf(lu)) {
            km1 = a + d + a*d;
        } else {
            double b = crho*lu;
            km1 = d + a*b*(1 + d)/(1 + a + b);
        }
        return log1p(norm_scale*km1);
    }
};

inline double branch_product(Tree &tree, double c, double rho, Node *cut) {
    double p = 1.0;
    for (auto &x : tree.parents) {
        double l = x.second->time - x.first->time;
        if (!isinf(l) and x.second != cut) {
            p *= 1 + c*rho*l;
        }
    }
    return p;
}

inline double called_width(Node *query_node, ARG &a, int i) {
    double w = a.assayed[i];
    if (!a.any_missing or query_node == nullptr or !query_node->has_missing()) {
        return w;
    }
    return query_node->called_width(i, w);
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

inline double update_branch_product(double p, Recombination &r, double crho, Node *cut) {
    for (const Branch &b : r.deleted_branches) {
        double l = b.upper_node->time - b.lower_node->time;
        if (!isinf(l) and b.upper_node != cut) {
            p /= 1 + crho*l;
        }
    }
    for (const Branch &b : r.inserted_branches) {
        double l = b.upper_node->time - b.lower_node->time;
        if (!isinf(l) and b.upper_node != cut) {
            p *= 1 + crho*l;
        }
    }
    return p;
}

#endif /* Emission_hpp */

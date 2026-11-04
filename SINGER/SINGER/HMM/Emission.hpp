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

struct Edge_matrices {
    bool at_root = false;
    double time = 0;
    double query_length = 0;
    double lower_length = 0;
    double upper_length = 0;
    double root_length = 0;
    double lower[16];
    double query[16];
    double upper[16];
};

class Emission {

public:

    double crho = 0;
    double norm_scale = 1.0;
    double bin_width = 1.0;
    const Site_model *model = nullptr;
    double unit = 0;
    double base_counts[4] = {0, 0, 0, 0};
    const Site_messages *base_groups[4] = {nullptr, nullptr, nullptr, nullptr};
    double unknown_count = 0;
    vector<const Site_messages *> sites = {};
    vector<const Site_messages *> groups = {};
    vector<double> group_weights = {};
    vector<pair<array<const Site_messages *, 4>, double>> unknown_groups = {};
    bool augment = false;
    bool unknown = false;
    const Site_messages *main_pattern = nullptr;
    vector<pair<const Site_messages *, array<double, 5>>> patterns = {};

    double mut_emit(Branch &branch, double time, Node *node, double norm);

    double mut_emit(Branch &branch, double norm, const Edge_matrices *cached);

    void null_logs_base(Branch &branch, const Edge_matrices *cached, double *out);

    void edge_matrices(Branch &branch, double time, Node *node, Edge_matrices *out);

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

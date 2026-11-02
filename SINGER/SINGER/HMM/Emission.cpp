//
//  Emission.cpp
//  SINGER
//

#include "Emission.hpp"

static double site_ratio_two(const Site_messages &s, Branch &branch, double time, Node *node) {
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
    return (w - s.empty_term)/s.scale;
}

static double ratio_with(int S, const Site_messages &s, int k, const Edge_matrices &e) {
    Msg d, nu, m, up;
    transmit(e.lower, S, s.in[k], d);
    transmit(e.query, S, s.nu, nu);
    for (int q = 0; q < S; q++) {
        m[q] = d[q]*nu[q];
    }
    const Msg &a = s.above[k];
    double w = 0;
    if (e.at_root) {
        for (int q = 0; q < S; q++) {
            w += a[q]*m[q];
        }
    } else {
        transmit(e.upper, S, m, up);
        for (int q = 0; q < S; q++) {
            w += a[q]*up[q];
        }
    }
    return (w - s.empty_term)/s.scale;
}

static double state_ratio(const Site_messages &s, int k, const Edge_matrices &e, const double *pi) {
    unsigned mask = s.lower_mask[k];
    int x = s.upper_state[k];
    double w = 0;
    for (int q = 0; q < 4; q++) {
        double top = x < 0 ? pi[q] : e.upper[4*x + q];
        double down = masked_sum(e.lower + 4*q, mask);
        double query = 0;
        for (int a = 0; a < 4; a++) {
            query += e.query[4*q + a]*s.nu[a];
        }
        w += top*down*query;
    }
    return w/s.old_factor[k];
}

static double mixture_with(int S, const Site_messages &s, int k, const Edge_matrices *edge) {
    if (s.parts.empty()) {
        return ratio_with(S, s, k, edge[0]);
    }
    double total = s.z, value = s.z*ratio_with(S, s, k, edge[0]);
    for (size_t c = 0; c < s.parts.size(); c++) {
        total += s.parts[c].z;
        value += s.parts[c].z*ratio_with(S, s.parts[c], k, edge[c + 1]);
    }
    return value/total;
}

void Emission::edge_matrices(Branch &branch, double time, Node *node, Edge_matrices *out) {
    int K = (int) model->rates.size();
    for (int c = 0; c < K; c++) {
        double u = unit*model->rates[c];
        Edge_matrices &e = out[c];
        e.at_root = branch.upper_node->index == -1;
        model->transition(u*(time - branch.lower_node->time), e.lower);
        model->transition(u*(time - node->time), e.query);
        if (!e.at_root) {
            model->transition(u*(branch.upper_node->time - time), e.upper);
        }
    }
}

double Emission::mut_emit(Branch &branch, double norm, const Edge_matrices *cached) {
    double emit = norm;
    int S = model->states, k = branch.lower_node->flat_index;
    for (const Site_messages *site : sites) {
        emit *= site->augmented ? state_ratio(*site, k, cached[site->category], model->pi.data()) : mixture_with(S, *site, k, cached);
    }
    for (size_t i = 0; i < groups.size(); i++) {
        emit *= pow(mixture_with(S, *groups[i], k, cached), group_weights[i]);
    }
    return emit;
}

double Emission::mut_emit(Branch &branch, double time, Node *node, double norm) {
    if (model->states == 2) {
        double emit = norm;
        for (const Site_messages *site : sites) {
            emit *= site_ratio_two(*site, branch, time, node);
        }
        return emit;
    }
    edge.resize(model->rates.size());
    edge_matrices(branch, time, node, edge.data());
    return mut_emit(branch, norm, edge.data());
}

void Emission::null_logs_base(Branch &branch, const Edge_matrices *cached, double *out) {
    int S = model->states, k = branch.lower_node->flat_index;
    for (int b = 0; b < 4; b++) {
        out[b] = -log(mixture_with(S, *base_groups[b], k, cached));
    }
}

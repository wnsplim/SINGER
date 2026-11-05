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
    transmit(e.lower, S, false, s.in[k], d);
    transmit(e.query, S, false, s.nu, nu);
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
        transmit(e.upper, S, false, m, up);
        for (int q = 0; q < S; q++) {
            w += a[q]*up[q];
        }
    }
    return (w - s.empty_term)/s.scale;
}

template <int N>
static double split_value(const Site_model &m, const Site_messages &s, const Edge_matrices &e, int b) {
    if (m.continuous_gamma()) {
        double A = m.alpha - m.rate(b, b)*(s.T - e.upper_length - e.lower_length);
        return (double) m.node_expectation(b, A, e.upper_length, e.lower_length, e.query_length, 0.0);
    }
    const int S = N > 0 ? N : m.states;
    double v = 0;
    for (int y = 0; y < S; y++) {
        v += e.upper[S*b + y]*e.lower[S*y + b]*e.query[S*y + b];
    }
    return v;
}

template <int N>
static double pattern_value(const Site_model &m, const Site_messages &s, int k, const Edge_matrices &e, int b) {
    const int S = N > 0 ? N : m.states;
    if (s.kind[k] == E2) {
        return s.fixed[S*k + b];
    }
    if (s.kind[k] == E1) {
        return split_value<N>(m, s, e, b)/s.fixed[S*k + b];
    }
    if (m.continuous_gamma()) {
        return m.rate_integrated_entry(b, m.alpha - m.rate(b, b)*s.T, e.root_length);
    }
    return e.upper[(S + 1)*b];
}

template <int N>
static double pattern_mixture(const Site_model &m, const Site_messages &s, int k, const Edge_matrices &e) {
    const int S = N > 0 ? N : m.states;
    double v = 0, w = 0;
    if (s.kind[k] == E1) {
        for (int b = 0; b < S; b++) {
            v += s.fixed[S*k + b]*split_value<N>(m, s, e, b);
        }
        return v;
    }
    for (int b = 0; b < S; b++) {
        v += s.weight[b]*pattern_value<N>(m, s, k, e, b);
        w += s.weight[b];
    }
    return v/w;
}

template <int N>
static double history_value(const Site_model &m, const Site_messages &s, int k, const Edge_matrices &e, bool unknown) {
    if (s.base >= 0) {
        return unknown ? pattern_mixture<N>(m, s, k, e) : pattern_value<N>(m, s, k, e, s.base);
    }
    const int S = N > 0 ? N : m.states;
    if (s.kind[k] == E1) {
        int x = s.upper_state[k], xl = s.lower_state[k];
        if (m.continuous_gamma()) {
            return m.node_sum(s.rate*e.upper_length, s.rate*e.lower_length, s.rate*e.query_length, x, xl, s.nu.data())/s.fixed[k];
        }
        double v = 0;
        for (int y = 0; y < S; y++) {
            const double *q = e.query + S*y;
            double d = 0;
            for (int a = 0; a < S; a++) {
                d += q[a]*s.nu[a];
            }
            v += e.upper[S*x + y]*e.lower[S*y + xl]*d;
        }
        return v/s.fixed[k];
    }
    if (s.kind[k] == ROOT) {
        int x = s.upper_state[k];
        if (m.continuous_gamma()) {
            return m.transition_dot(s.rate*e.root_length, x, s.nu.data());
        }
        const double *row = e.upper + S*x;
        double d = 0;
        for (int a = 0; a < S; a++) {
            d += row[a]*s.nu[a];
        }
        return d;
    }
    return s.fixed[k];
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

static double total_z(const Site_messages &s) {
    double z = s.z;
    for (const Site_messages &part : s.parts) {
        z += part.z;
    }
    return z;
}

static double mixture_unknown(int S, const array<const Site_messages *, max_states> &set, int k, const Edge_matrices *cached) {
    double with = 0, without = 0;
    for (int b = 0; b < S; b++) {
        double z = total_z(*set[b]);
        with += z*mixture_with(S, *set[b], k, cached);
        without += z;
    }
    return with/without;
}

void Emission::edge_matrices(Branch &branch, double time, Node *node, Edge_matrices *out) {
    if (augment) {
        Edge_matrices &e = out[0];
        e.at_root = branch.upper_node->index == -1;
        e.time = time;
        e.query_length = unit*(time - node->time);
        e.lower_length = unit*(time - branch.lower_node->time);
        e.upper_length = e.at_root ? 0.0 : unit*(branch.upper_node->time - time);
        e.root_length = unit*(2*time - branch.lower_node->time - node->time);
        if (model->continuous_gamma()) {
            return;
        }
        model->transition(e.query_length, e.query);
        if (e.at_root) {
            model->transition(e.root_length, e.upper);
        } else {
            model->transition(e.lower_length, e.lower);
            model->transition(e.upper_length, e.upper);
        }
        return;
    }
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

template <int N>
static double history_emit(const Emission &h, int k, double emit, const Edge_matrices &e) {
    const Site_model &m = *h.model;
    const int S = N > 0 ? N : m.states;
    for (const Site_messages *site : h.sites) {
        if (!site->silent) {
            emit *= history_value<N>(m, *site, k, e, h.unknown);
        }
    }
    for (auto &[p, n] : h.patterns) {
        if (p->silent) {
            continue;
        }
        if (h.unknown) {
            emit *= pow(pattern_mixture<N>(m, *p, k, e), n[S]);
            continue;
        }
        for (int b = 0; b < S; b++) {
            if (n[b] > 0) {
                emit *= pow(pattern_value<N>(m, *p, k, e, b), n[b]);
            }
        }
    }
    return emit;
}

double Emission::mut_emit(Branch &branch, double norm, const Edge_matrices *cached) {
    double emit = norm;
    int S = model->states, k = branch.lower_node->flat_index;
    if (augment) {
        return S == 4 ? history_emit<4>(*this, k, emit, cached[0]) : S == 2 ? history_emit<2>(*this, k, emit, cached[0]) : history_emit<0>(*this, k, emit, cached[0]);
    }
    for (const Site_messages *site : sites) {
        emit *= mixture_with(S, *site, k, cached);
    }
    for (size_t i = 0; i < groups.size(); i++) {
        emit *= pow(mixture_with(S, *groups[i], k, cached), group_weights[i]);
    }
    for (auto &u : unknown_groups) {
        emit *= pow(mixture_unknown(S, u.first, k, cached), u.second);
    }
    return emit;
}

double Emission::mut_emit(Branch &branch, double time, Node *node, double norm) {
    double emit = norm;
    for (const Site_messages *site : sites) {
        emit *= site_ratio_two(*site, branch, time, node);
    }
    return emit;
}

void Emission::null_logs_base(Branch &branch, const Edge_matrices *cached, double *out) {
    int S = model->states, k = branch.lower_node->flat_index;
    if (augment) {
        const Site_messages &m = *main_pattern;
        for (int b = 0; b < S; b++) {
            out[b] = m.silent or unknown ? 0.0 : -log(pattern_value<0>(*model, m, k, cached[0], b));
        }
        out[S] = unknown and !m.silent ? -log(pattern_mixture<0>(*model, m, k, cached[0])) : 0.0;
        return;
    }
    for (int b = 0; b < S; b++) {
        out[b] = -log(mixture_with(S, *base_groups[b], k, cached));
    }
    out[S] = 0;
    if (unknown_count > 0) {
        array<const Site_messages *, max_states> set;
        copy(base_groups, base_groups + max_states, set.begin());
        out[S] = -log(mixture_unknown(S, set, k, cached));
    }
}

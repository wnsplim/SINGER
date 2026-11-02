//
//  Scaler.cpp
//  SINGER
//
//  Created by Yun Deng on 1/3/24.
//

#include "Scaler.hpp"

Scaler::Scaler() {}

void Scaler::reset() {
    rates.clear();
    accumulated_arg_length.clear();
    old_grid.assign(1, 0);
    new_grid.assign(1, 0);
    expected_arg_length.clear();
    observed_arg_length.clear();
    scaling_factors.clear();
}

void Scaler::compute_deltas(ARG &a, Rate_map &mm, double start, double m) {
    auto span_length = [&](double x, double y) {
        return mm.coordinates.empty() ? y - x : mm.segment_distance(start + x, start + y)/m;
    };
    map<Node *, double, compare_node> delta;
    a.create_node_set();
    for (Node *n : a.node_set) {
        delta[n] = 0;
    }
    a.informed_edges([&](Node *child, Node *parent, double x0, double x1) {
        double w = span_length(x0, x1);
        delta[child] += w;
        delta[parent] -= w;
    });
    heterochronous = (*a.sample_nodes.rbegin())->time > 0;
    sorted_nodes.clear();
    node_deltas.clear();
    for (Node *n : a.sample_nodes) {
        sorted_nodes.push_back(n);
        node_deltas.push_back(delta[n]);
    }
    for (auto &x : delta) {
        if (!x.first->is_sample) {
            sorted_nodes.push_back(x.first);
            node_deltas.push_back(x.second);
        }
    }
    children.clear();
    if (heterochronous) {
        vector<int> order(sorted_nodes.size());
        iota(order.begin(), order.end(), 0);
        compare_node cmp;
        sort(order.begin(), order.end(), [&](int i, int j) { return cmp(sorted_nodes[i], sorted_nodes[j]); });
        vector<Node *> nodes_sorted;
        vector<double> deltas_sorted;
        for (int i : order) {
            nodes_sorted.push_back(sorted_nodes[i]);
            deltas_sorted.push_back(node_deltas[i]);
        }
        sorted_nodes = nodes_sorted;
        node_deltas = deltas_sorted;
        for (auto &x : a.recombinations) {
            for (const Branch &b : x.second.inserted_branches) {
                if (b.upper_node != a.root.get()) {
                    children[b.upper_node].push_back(b.lower_node);
                }
            }
        }
    }
    times.clear();
    position.clear();
    for (int i = 0; i < (int) sorted_nodes.size(); i++) {
        times.push_back(sorted_nodes[i]->time);
        position[sorted_nodes[i]] = i;
    }
}

void Scaler::compute_old_grid() {
    expected_arg_length.resize(num_windows);
    rates.resize(sorted_nodes.size());
    accumulated_arg_length.resize(sorted_nodes.size());
    partial_sum(node_deltas.begin(), node_deltas.end(), rates.begin());
    assert(rates.back() > -1);
    for (int i = 1; i < sorted_nodes.size(); i++) {
        accumulated_arg_length[i] = accumulated_arg_length[i-1] + rates[i-1]*(times[i] - times[i-1]);
    }
    double unit_arg_length = accumulated_arg_length.back()/num_windows;
    double partial_arg_length = 0;
    int new_index = 0;
    double rate = 0;
    double residue = 0;
    for (int i = 1; i <= num_windows; i++) {
        partial_arg_length = accumulated_arg_length.back()*i/num_windows;
        auto it = upper_bound(accumulated_arg_length.begin(), accumulated_arg_length.end(), partial_arg_length);
        new_index = (int) distance(accumulated_arg_length.begin(), it);
        new_index = min((int) sorted_nodes.size() - 1, new_index);
        rate = rates[new_index - 1];
        residue = accumulated_arg_length[new_index] - partial_arg_length;
        residue = max(0.0, residue);
        expected_arg_length[i-1] = unit_arg_length;
        old_grid.push_back(times[new_index] - residue/rate);
    }
    old_grid.back() = nextafter(times.back(), INT_MAX);
    assert(old_grid.size() == num_windows + 1);
}

void Scaler::compute_new_grid(double theta) {
    for (auto &x : observed_arg_length) {
        x /= theta;
    }
    double base_time = 0;
    double old_window_width = 0, scaling_factor = 0;
    for (int i = 1; i < old_grid.size(); i++) {
        old_window_width = old_grid[i] - old_grid[i-1];
        scaling_factor = observed_arg_length[i-1]/expected_arg_length[i-1];
        scaling_factors.push_back(scaling_factor);
        base_time += old_window_width*scaling_factor;
        new_grid.push_back(base_time);
    }
}

void Scaler::map_mutations(ARG &a) {
    observed_arg_length.assign(num_windows, 0.0);
    function<double(Node *)> time_of = [&](Node *n) { return n == a.root.get() ? n->time : times[position.at(n)]; };
    a.expected_counts(&time_of, [&](Node *lower, Node *upper, double e) { add_mutation(e, time_of(lower), time_of(upper)); }, nullptr, nullptr);
}

void Scaler::add_mutation(double w, double lb, double ub) {
    if (isinf(ub)) {
        return;
    }
    double x, y, l;
    int index;
    auto it = upper_bound(old_grid.begin(), old_grid.end(), lb);
    it--;
    index = (int) distance(old_grid.begin(), it);
    while (old_grid[index] < ub) {
        x = old_grid[index];
        y = old_grid[index + 1];
        l = min(ub, y) - max(lb, x);
        observed_arg_length[index] += w*l/(ub - lb);
        index++;
    }
}

void Scaler::rescale(ARG &a, double theta) {
    compute_old_grid();
    map_mutations(a);
    compute_new_grid(theta);
    int k = 0;
    double old_t = 0;
    double new_t = 0;
    for (int i = 0; i < sorted_nodes.size(); i++) {
        if (sorted_nodes[i]->is_sample) {
            continue;
        }
        while (times[i] > old_grid[k+1]) {
            k++;
            assert(k <= num_windows);
        }
        new_t = scaling_factors[k]*(times[i] - old_grid[k]) + new_grid[k];
        assert(new_t >= old_t);
        times[i] = new_t;
        old_t = new_t;
    }
    if (heterochronous) {
        for (int i = 0; i < sorted_nodes.size(); i++) {
            if (sorted_nodes[i]->is_sample) {
                continue;
            }
            double lb = 0;
            for (Node *c : children[sorted_nodes[i]]) {
                lb = max(lb, times[position.at(c)]);
            }
            if (times[i] <= lb) {
                times[i] = lb + 1e-6;
            }
        }
    }
}

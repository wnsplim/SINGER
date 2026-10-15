//
//  Node.cpp
//  SINGER
//
//  Created by Yun Deng on 3/31/22.
//

#include "Node.hpp"

Node::Node(double t) {
    time = t;
}

void Node::set_index(int index) {
    this->index = index;
}

static inline bool site_less(const pair<double, double> &e, double p) {
    return e.first < p;
}

void Node::set_site(double pos, double s) {
    auto lo = lower_bound(mutation_sites.begin(), mutation_sites.end(), pos, site_less);
    if (lo != mutation_sites.end() and lo->first == pos) {
        lo->second = s;
        return;
    }
    size_t at = lo - mutation_sites.begin();
    size_t cur = it - mutation_sites.begin();
    mutation_sites.insert(lo, {pos, s});
    it = mutation_sites.begin() + (cur >= at ? cur + 1 : cur);
}

void Node::add_mutation(double pos) {
    set_site(pos, 1);
}

void Node::add_missing(double pos) {
    missing_sites.push_back(pos);
}

bool Node::is_missing(double pos) {
    if (binary_search(missing_sites.begin(), missing_sites.end(), pos)) {
        return true;
    }
    auto it = upper_bound(masked_intervals.begin(), masked_intervals.end(), make_pair(pos, numeric_limits<double>::infinity()));
    return it != masked_intervals.begin() and prev(it)->second > pos;
}

vector<double> Node::table_sites = {};

int Node::site_index(double pos) {
    auto it = lower_bound(table_sites.begin(), table_sites.end(), pos);
    return it != table_sites.end() and *it == pos ? (int) (it - table_sites.begin()) : -1;
}

void Node::build_missing_table() {
    missing_table.resize(table_sites.size());
    for (size_t s = 0; s < table_sites.size(); s++) {
        missing_table[s] = is_missing(table_sites[s]);
    }
}

double Node::called_width(int bin, double assayed) {
    while (width_at > 0 and width_table[width_at - 1].first >= bin) {
        width_at--;
    }
    while (width_at < width_table.size() and width_table[width_at].first < bin) {
        width_at++;
    }
    return width_at < width_table.size() and width_table[width_at].first == bin ? width_table[width_at].second : assayed;
}

bool Node::likelihood_at(double pos, double &w0, double &w1) {
    auto it = lower_bound(likelihood_sites.begin(), likelihood_sites.end(), pos,
                          [](const pair<double, pair<double, double>> &e, double p) { return e.first < p; });
    if (it == likelihood_sites.end() or it->first != pos) {
        return false;
    }
    w0 = it->second.first;
    w1 = it->second.second;
    return true;
}

void Node::write_state(double pos, double s) {
    if (s == 0) {
        auto lo = lower_bound(mutation_sites.begin(), mutation_sites.end(), pos, site_less);
        if (lo == mutation_sites.end() or lo->first != pos) {
            return;
        }
        bool on_it = (it == lo);
        size_t at = lo - mutation_sites.begin();
        size_t cur = it - mutation_sites.begin();
        mutation_sites.erase(lo);
        if (on_it) {
            it = mutation_sites.begin();
        } else {
            it = mutation_sites.begin() + (cur > at ? cur - 1 : cur);
        }
        return;
    } else if (s == 1) {
        set_site(pos, s);
    }
    return;
}

shared_ptr<Node> new_node(double t) {
    return make_shared<Node>(t);
}

void Node::move_iterator(double m) {
    if (it->first == m) {
        return;
    }
    if (it->first < m) {
        double next_pos = next(it)->first;
        if (next_pos == m) {
            ++it;
            return;
        } else if (next_pos > m) {
            return;
        }
    } else if (it->first > m) {
        double prev_pos = prev(it)->first;
        if (prev_pos <= m) {
            --it;
            return;
        }
    }
    if (abs(it->first - m) < 20) {
        while (it != mutation_sites.begin() and it->first > m) {
            --it;
        }
        while (next(it) != mutation_sites.end() and next(it)->first <= m) {
            ++it;
        }
    } else {
        it = upper_bound(mutation_sites.begin(), mutation_sites.end(), m,
                         [](double p, const pair<double, double> &e) { return p < e.first; });
        --it;
    }
    assert(it->first <= m and next(it)->first > m);
}

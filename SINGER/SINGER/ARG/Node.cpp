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
    return binary_search(missing_sites.begin(), missing_sites.end(), pos);
}

double Node::get_state(double pos) {
    move_iterator(pos);
    if (it->first == pos) {
        return it->second;
    } else {
        return 0;
    }
    return 0;
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

void Node::read_mutation(string filename) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    double x;
    while (fin >> x) {
        add_mutation(x);
    }
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

/*
Node::Node(double t) {
    time = t;
}

void Node::set_index(int index) {
    this->index = index;
}

void Node::add_mutation(double pos) {
    mutation_sites.insert(pos);
}

double Node::get_state(double pos) {
    if (mutation_sites.count(pos) > 0) {
        return 1;
    } else {
        return 0;
    }
}

void Node::write_state(double pos, double s) {
    if (s == 0) {
        mutation_sites.erase(pos);
        // ambiguous_sites.erase(pos);
        return;
    } else if (s == 1) {
        // ambiguous_sites.erase(pos);
        mutation_sites.insert(pos);
    }
    return;
}

void Node::read_mutation(string filename) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    double x;
    while (fin >> x) {
        add_mutation(x);
    }
}

shared_ptr<Node> new_node(double t) {
    return make_shared<Node>(t);
}
*/

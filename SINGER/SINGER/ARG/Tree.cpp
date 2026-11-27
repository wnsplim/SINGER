//
//  Tree.cpp
//  SINGER
//
//  Created by Yun Deng on 4/12/22.
//

#include "Tree.hpp"

Tree::Tree() {
}

double Tree::length() {
    double l = 0;
    for (auto &x : parents) {
        if (x.second->index != -1) {
            l += x.second->time - x.first->time;
        }
    }
    return l;
}

void Flat_tree::assign(const Tree &tree) {
    parents.assign(tree.parents.begin(), tree.parents.end());
    lengths.clear();
    for (auto &x : parents) {
        lengths.push_back(x.second->time - x.first->time);
    }
}

void Flat_tree::forward_update(const Recombination &r) {
    auto key_less = [](const pair<Node *, Node *> &e, Node *n) { return compare_node()(e.first, n); };
    for (const Branch &b : r.deleted_branches) {
        if (any_of(r.inserted_branches.begin(), r.inserted_branches.end(), [&b](const Branch &c) { return c.lower_node == b.lower_node; })) {
            continue;
        }
        auto it = lower_bound(parents.begin(), parents.end(), b.lower_node, key_less);
        if (it != parents.end() and it->first == b.lower_node) {
            lengths.erase(lengths.begin() + (it - parents.begin()));
            parents.erase(it);
        }
    }
    for (const Branch &b : r.inserted_branches) {
        auto it = lower_bound(parents.begin(), parents.end(), b.lower_node, key_less);
        double l = b.upper_node->time - b.lower_node->time;
        if (it != parents.end() and it->first == b.lower_node) {
            it->second = b.upper_node;
            lengths[it - parents.begin()] = l;
        } else {
            lengths.insert(lengths.begin() + (it - parents.begin()), l);
            parents.insert(it, {b.lower_node, b.upper_node});
        }
    }
}

double Flat_tree::length() const {
    double l = 0;
    for (int i = 0; i < (int) parents.size(); i++) {
        if (parents[i].second->index != -1) {
            l += lengths[i];
        }
    }
    return l;
}

static vector<map<Node *, Node *, compare_node>::node_type> spare_nodes;

void Tree::delete_branch(const Branch &b) {
    assert(b.upper_node != nullptr and b.lower_node != nullptr);
    auto nh = parents.extract(b.lower_node);
    if (nh and spare_nodes.size() < 64) {
        spare_nodes.push_back(move(nh));
    }
}

void Tree::insert_branch(const Branch &b) {
    assert(b.upper_node != nullptr and b.lower_node != nullptr);
    auto it = parents.lower_bound(b.lower_node);
    if (it != parents.end() and !parents.key_comp()(b.lower_node, it->first)) {
        it->second = b.upper_node;
    } else if (spare_nodes.empty()) {
        parents.emplace_hint(it, b.lower_node, b.upper_node);
    } else {
        auto nh = move(spare_nodes.back());
        spare_nodes.pop_back();
        nh.key() = b.lower_node;
        nh.mapped() = b.upper_node;
        parents.insert(it, move(nh));
    }
}

void Tree::forward_update(Recombination &r) {
    int prev_size = (int) parents.size();
    for (const Branch &b : r.deleted_branches) {
        delete_branch(b);
    }
    for (const Branch &b : r.inserted_branches) {
        insert_branch(b);
    }
    assert(prev_size == (int) parents.size() or r.pos == 0);
}

void Tree::backward_update(Recombination &r) {
    for (const Branch &b : r.inserted_branches) {
        delete_branch(b);
    }
    for (const Branch &b : r.deleted_branches) {
        insert_branch(b);
    }
}

void Tree::remove(Branch b, Node *n) {
    assert(b.upper_node->index >= 0);
    Branch joining_branch = find_joining_branch(b);
    Node *sibling = find_sibling(b.lower_node);
    Node *parent = parents[b.upper_node];
    delete_branch(b);
    delete_branch(Branch(sibling, b.upper_node));
    delete_branch(Branch(b.upper_node, parent));
    insert_branch(joining_branch);
    insert_branch(Branch(b.lower_node, n));
}

void Tree::add(Branch added_branch, Branch joining_branch, Node *n) {
    if (n != nullptr) {
        delete_branch(Branch(added_branch.lower_node, n));
    }
    delete_branch(joining_branch);
    insert_branch(Branch(joining_branch.lower_node, added_branch.upper_node));
    insert_branch(Branch(added_branch.upper_node, joining_branch.upper_node));
    insert_branch(added_branch);
}

Node *Tree::find_sibling(Node *n) {
    Node *p = parents[n];
    for (auto &x : parents) {
        if (x.second == p and x.first != n) {
            return x.first;
        }
    }
    return nullptr;
}

Branch Tree::find_joining_branch(Branch removed_branch) {
    if (removed_branch == Branch()) {
        return Branch();
    }
    Node *p = parents[removed_branch.upper_node];
    Node *c = find_sibling(removed_branch.lower_node);
    assert(parents[c] == removed_branch.upper_node);
    return Branch(c, p);
}

pair<Branch, double> Tree::sample_uniform_cut_point() {
    double target = uniform_random()*length();
    double acc = 0;
    Branch last;
    for (auto &x : parents) {
        if (x.second->index == -1) {
            continue;
        }
        double bl = x.second->time - x.first->time;
        if (acc + bl >= target) {
            return {Branch(x.first, x.second), x.first->time + (target - acc)};
        }
        acc += bl;
        last = Branch(x.first, x.second);
    }
    return {last, 0.5*(last.lower_node->time + last.upper_node->time)};
}

double Tree::prior_likelihood() {
    vector<pair<double, int>> events;
    for (auto &x : parents) {
        events.push_back({x.first->time, x.first->is_sample ? 1 : -1});
    }
    sort(events.begin(), events.end());
    double log_likelihood = 0;
    double prev = 0;
    int k = 0;
    for (auto &e : events) {
        log_likelihood -= 0.5*k*(k - 1)*(e.first - prev);
        k += e.second;
        prev = e.first;
    }
    return log_likelihood;
}

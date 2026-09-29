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

void Tree::delete_branch(const Branch &b) {
    assert(b.upper_node != nullptr and b.lower_node != nullptr);
    parents.erase(b.lower_node);
    unordered_set<Node *> &children_nodes = children[b.upper_node];
    if (children_nodes.size() == 1) {
        children.erase(b.upper_node);
    } else {
        children_nodes.erase(b.lower_node);  
    }
}

void Tree::insert_branch(const Branch &b) {
    assert(b.upper_node != nullptr and b.lower_node != nullptr);
    parents[b.lower_node] = b.upper_node;
    children[b.upper_node].insert(b.lower_node);
}

void Tree::forward_update(Recombination &r) {
    int prev_size = (int) parents.size();
    for (const Branch &b : r.deleted_branches) {
        delete_branch(b);
    }
    for (const Branch &b : r.inserted_branches) {
        insert_branch(b);
    }
    int after_size = (int) parents.size();
    assert(prev_size == after_size or r.pos == 0);
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
    Branch sibling_branch = Branch(sibling, b.upper_node);
    Branch parent_branch = Branch(b.upper_node, parent);
    Branch cut_branch = Branch(b.lower_node, n);
    delete_branch(b);
    delete_branch(sibling_branch);
    delete_branch(parent_branch);
    insert_branch(joining_branch);
    insert_branch(cut_branch);
}

void Tree::add(Branch added_branch, Branch joining_branch, Node *n) {
    Branch lower_branch = Branch(joining_branch.lower_node, added_branch.upper_node);
    Branch upper_branch = Branch(added_branch.upper_node, joining_branch.upper_node);
    if (n != nullptr) {
        Branch cut_branch = Branch(added_branch.lower_node, n);
        delete_branch(cut_branch);
    }
    delete_branch(joining_branch);
    insert_branch(lower_branch);
    insert_branch(upper_branch);
    insert_branch(added_branch);
}

Node *Tree::find_sibling(Node *n) {
    Node *p = parents[n];
    unordered_set<Node *> &candidates = children[p];
    Node *c;
    auto c_it = candidates.begin();
    if (*c_it != n) {
        c = *c_it;
    } else {
        c = *(next(c_it));
    }
    return c;
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

pair<Branch, double> Tree::sample_cut_point() {
    double root_time = parents.rbegin()->first->time;
    double cut_time = random()*root_time;
    vector<Branch> candidates = {};
    for (auto &x : parents) {
        if (x.second->time > cut_time and x.first->time <= cut_time) {
            candidates.push_back(Branch(x.first, x.second));
        }
    }
    int index = (int) floor(candidates.size()*uniform_random());
    index = min((int) candidates.size() - 1, index);
    return {candidates[index], cut_time};
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

double Tree::log_exp(double lambda, double x) {
    return -lambda*x + log(lambda);
}

void Tree::impute_states(double m, set<Branch> &mutation_branches) {
    map<Node *, double> states = {};
    for (const Branch &b : mutation_branches) {
        states[b.lower_node] = b.lower_node->get_state(m);
        states[b.upper_node] = b.upper_node->get_state(m);
    }
    for (auto &x : parents) {
        impute_states_helper(x.first, states);
    }
    for (auto &x : states) {
        x.first->write_state(m, x.second);
    }
}

void Tree::impute_states_helper(Node *n, map<Node *, double> &states) {
    if (n->index == -1) {
        states[n] = 0;
        return;
    }
    if (states.count(n) > 0) {
        return;
    }
    Node *p = parents[n];
    impute_states_helper(p, states);
    states[n] = states[p];
}

double Tree::random() {
    double p = uniform_random();
    return p;
}

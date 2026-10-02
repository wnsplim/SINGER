//
//  Tree.hpp
//  SINGER
//
//  Created by Yun Deng on 4/12/22.
//

#ifndef Tree_hpp
#define Tree_hpp

#include <stdio.h>
#include <math.h>
#include <map>
#include <unordered_map>
#include "random_utils.hpp"
#include "Branch.hpp"
#include "Recombination.hpp"

class Tree {

public:
    
    map<Node *, Node *, compare_node> parents = {};
    
    Tree();
    
    double length();
    
    void insert_branch(const Branch &b);
    
    void delete_branch(const Branch &b);

    void forward_update(Recombination &r);
    
    void backward_update(Recombination &r);
    
    void remove(Branch b, Node *n);
    
    void add(Branch added_branch, Branch joining_branch, Node *n);
    
    Node *find_sibling(Node *n);
    
    Branch find_joining_branch(Branch removed_branch);
    
    pair<Branch, double> sample_cut_point();

    pair<Branch, double> sample_uniform_cut_point();

    double prior_likelihood();

    double log_exp(double lambda, double x);

    void impute_states(double m, set<Branch> &mutation_branches);
    
    void impute_states_helper(Node *n, map<Node *, double> &states);
    
    double random();

};

class Flat_tree {

public:

    vector<pair<Node *, Node *>> parents = {};
    vector<double> lengths = {};

    void assign(const Tree &tree);

    void forward_update(const Recombination &r);

    double length() const;

};

#endif /* Tree_hpp */

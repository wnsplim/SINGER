//
//  ARG.hpp
//  SINGER
//
//  Created by Yun Deng on 4/14/22.
//

#ifndef ARG_hpp
#define ARG_hpp

#include <stdio.h>
#include <map>
#include "Recombination.hpp"
#include "Tree.hpp"
#include "RSP_smc.hpp"
#include "Rate_map.hpp"
#include <array>
#include <tuple>
#include <unordered_map>
#include <functional>

typedef array<double, 2> Msg;

const Msg ones = [] { Msg m; m.fill(1.0); return m; }();

inline void push_up(double x, const Msg &m, Msg &parent) {
    parent[0] *= m[0] + x*m[1];
    parent[1] *= m[1] + x*m[0];
}

inline void transmit(double x, const Msg &m, Msg &out) {
    out[0] = m[0] + x*m[1];
    out[1] = m[1] + x*m[0];
}

struct Genotype_likelihood {
    double pos;
    Node *leaves[2];
    double L[4];
};

struct Mutation_row {
    Branch branch;
    int state;
    char informed;
};

struct Cell_row {
    double pos;
    int leaf[2];
    char state[2];
};

struct Bin_gaps {
    bool built = false;
    vector<int> sample = {};
    vector<int> d = {};
    vector<double> at = {};
    vector<double> before = {};
    double tail = 0;
    bool uniform = false;
    vector<int> full = {};
    double count = 0;
};

struct Tree_order {
    vector<double> times = {};
    vector<int> up = {};
    vector<vector<int>> children = {};
    int top = -1;
    vector<int> leaf = {};

    void assign(const Flat_tree &flat, Node *root_node);

    int index_of(const Flat_tree &flat, Node *n) const;
};

struct Peel {
    vector<char> missing = {};
    vector<Msg> msg = {};
    vector<Msg> empty = {};
    int count[2] = {};
    int leaf[2] = {-1, -1};
    vector<int> path = {};
    vector<char> in_path = {};
    vector<Msg> alt = {};
    vector<Msg> outside = {};
    vector<char> in_main = {};
};

struct Site_messages {
    double unit = 0;
    vector<Msg> in = {};
    vector<Msg> above = {};
    Msg nu = ones;
    int tree = -1;
    double empty_term = 0;
    double scale = 1;
    vector<Site_messages> alt = {};
};

class ARG {

public:

    double Ne = 1;
    double time_offset = 0;
    Node_ptr root = new_node(numeric_limits<double>::infinity());
    Node *cut_node = nullptr;
    Node_ptr cut_node_owner = nullptr;
    double cut_time = 0;
    set<double> mutation_sites = {};
    map<double, vector<Mutation_row>> mutations = {};
    map<double, Recombination> recombinations = {};
    int bin_num = 0;
    double penalty = 1;
    double ancestral_prob = 0.5;
    bool any_missing = false;
    vector<double> unassayed_sites = {};
    double sequence_length = 0;
    vector<double> coordinates = {};
    vector<double> rhos = {};
    vector<double> thetas = {};
    vector<Node_ptr> node_owner = {};
    vector<Node_ptr> dead_nodes = {};
    vector<Node *> node_list = {};
    int num_samples = 0;
    void release_dead_nodes();
    set<Node *, compare_node> sample_nodes = {};
    set<Node *, compare_node> node_set = {};
    map<double, Branch> joining_branches = {};
    map<double, Branch> removed_branches = {};

    double start = 0;
    double end = 0;
    double cut_pos = 0;
    Tree cut_tree;
    Tree start_tree;
    Tree end_tree;

    ARG();

    ARG(double N, double l);

    ~ARG();

    void discretize(double s);

    void discretize(Rate_map &rm, Rate_map &mm, double start, double unit, double r);

    int get_index(double x);

    void compute_rhos_thetas(double r, double m, Rate_map &rm, Rate_map &mm, double start);

    void build_singleton_arg(const Node_ptr &n);

    void add_sample(const Node_ptr &n);

    void add_node(Node *n);

    void add_new_node(double t, bool sample);

    Tree get_tree_at(double x);

    map<double, Tree> anchors = {};
    vector<pair<double, double>> anchor_dirty = {};
    int anchor_step = 0;

    void anchors_changed(double lo, double hi);

    void rebuild_anchors();

    void window_copy_into(ARG &c, double lo_pos);

    Node *get_query_node_at(double x);

    void remove(tuple<double, Branch, double> cut_point);

    void add(map<double, Branch> &new_joining_branches, map<double, Branch> &added_branches);

    bool empties_cut_span_end_recombination(map<double, Branch> &new_joining_branches, map<double, Branch> &added_branches);

    void smc_sample_recombinations(map<double, Branch> &lineage);

    void approx_sample_recombinations();

    void read_coordinates(string filename);

    void write_coordinates(string filename);

    void read_sample(string filename);

    double get_arg_length();

    tuple<double, Branch, double> sample_internal_cut();

    tuple<double, Branch, double> sample_uniform_cut();

    vector<Genotype_likelihood> genotype_likelihoods = {};

    unordered_map<Node *, vector<int>> leaf_entries = {};

    void index_genotype_likelihoods();

    void index_missing(const vector<Node *> &leaves);

    vector<pair<double, int>> site_list = {};

    vector<pair<double, int>>::iterator site_at(double x);

    pair<const int *, const int *> entries_in(Node *leaf, double x, double y) const;

    void site_missing(const Flat_tree &flat, double pos, vector<char> &missing);

    void root_prior(double *prior) const;

    Msg leaf_message(int code) const;

    void peel_site(const Flat_tree &flat, const Tree_order &t, double unit, double pos, Peel &p, Node *partner = nullptr, int partner_state = 0, Node *leaf = nullptr, const Msg &tip = ones);

    double root_value(const Msg &top, bool seen0, bool seen1) const;

    void emission_site(const Flat_tree &flat, const Tree_order &t, double unit, double pos, Node *query, Peel &p, Site_messages &s, Node *partner = nullptr, int partner_state = 0, Node *leaf = nullptr, const Msg &tip = ones);

    void cell_sum_site(const Flat_tree &flat, const Tree_order &t, double unit, const Genotype_likelihood &g, Node *leaf, Node *query, Peel &p, Site_messages &s);

    void sparse_site(const Flat_tree &flat, const Tree_order &t, double unit, double pos, Node *query, Node *removed, Node *joined, double cut, Peel &p, Site_messages &s, const vector<Msg> &full, double empty);

    void graft_single(const Flat_tree &flat, const Tree_order &t, double unit, const vector<Msg> &in, const Branch &added, const Branch &joined, vector<Msg> &out) const;

    void entry_weights(const Flat_tree &flat, const Tree_order &t, double unit, double pos, const Genotype_likelihood &g, Peel &p, double *w);

    void set_entry(const Flat_tree &flat, const Tree_order &t, double unit, double pos, const Genotype_likelihood &g, int c, Peel &p);

    void draw_states(const Flat_tree &flat, const Tree_order &t, double unit, const Peel &p, vector<int> &state);

    void resample_site(const Flat_tree &flat, const Tree_order &t, int bin, double pos, size_t e, size_t f, Node *only, Peel &p);

    void resample_entries(Node *leaf, double x, double y);

    void resample_genotypes();

    void set_lengths(Flat_tree &flat, const function<double(Node *)> *time_of);

    void outside_messages(const Flat_tree &flat, const Tree_order &t, double unit, const Peel &p, vector<Msg> &above);

    void site_counts(const Flat_tree &flat, const Tree_order &t, double unit, const Peel &p, const vector<Msg> &above, vector<double> &differs, vector<char> &informed_tip);

    void draw_mutations(const function<double(Node *)> *time_of, map<double, vector<Mutation_row>> &out);

    vector<tuple<int, int, double, double>> branch_rows();

    void informed_edges(const function<void(Node *, Node *, double, double)> &add);

    void expected_counts(const function<double(Node *)> *time_of, const function<void(Node *, Node *, double)> &add, vector<double> *expected, const vector<tuple<int, int, double, double>> *rows);

    void mapping_counts(const function<double(Node *)> *time_of, vector<double> &expected, vector<double> &informed_bases);

    vector<double> assayed = {};
    vector<pair<double, double>> masked = {};

    void compute_assayed();

    void clear_remove_info();

    double corrected_smc_prior(map<double, Branch> &lineage, int lo, int hi, double p, const Tree &tree_p);

    Flat_tree likelihood_tree;

    void reference_peel(const Flat_tree &flat, const Tree_order &t, double unit, vector<Msg> &msg) const;

    void repeel_path(const Flat_tree &flat, const Tree_order &t, double unit, int k, vector<Msg> &msg) const;

    mutable vector<Bin_gaps> bin_gaps = {};

    const Bin_gaps &gaps_of(int bin) const;

    Branch lineage_branch_before(map<double, Branch> &lineage, double x);

    set<double> get_check_points();

    bool check_disjoint_nodes(double x, double y);

    void new_recombination(double pos, Branch prev_added_branch, Branch prev_joining_branch, Branch next_added_branch, Branch next_joining_branch);

    void remove_empty_recombinations();

    void create_node_set();

    void number_nodes();

    vector<Cell_row> cell_rows();

    void write_sample(string filename, const vector<Cell_row> &cells, const vector<double> &rescaled, const vector<map<double, vector<Mutation_row>>> &draws, const vector<double> &expected, const vector<double> &informed_bases);

    void read_nodes(const vector<double> &times);

    void read_branches(const vector<double> &v);

    void read_recombs(const vector<double> &v);

    void read_cells(string filename);

};

bool compare_edge(const tuple<int, int, double, double>& edge1, const tuple<int, int, double, double>& edge2);


#endif /* ARG_hpp */

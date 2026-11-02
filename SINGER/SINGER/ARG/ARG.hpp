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
#include "Site_model.hpp"
#include <array>
#include <tuple>
#include <unordered_map>
#include <functional>

typedef array<double, 4> Msg;

inline void push_up(const double *T, int S, const Msg &m, Msg &parent) {
    if (S == 2) {
        parent[0] *= m[0] + T[1]*m[1];
        parent[1] *= m[1] + T[2]*m[0];
        return;
    }
    for (int s = 0; s < S; s++) {
        double v = 0;
        for (int t = 0; t < S; t++) {
            v += T[S*s + t]*m[t];
        }
        parent[s] *= v;
    }
}

inline void transmit(const double *T, int S, const Msg &m, Msg &out) {
    if (S == 2) {
        out[0] = m[0] + T[1]*m[1];
        out[1] = m[1] + T[2]*m[0];
        return;
    }
    for (int s = 0; s < S; s++) {
        double v = 0;
        for (int t = 0; t < S; t++) {
            v += T[S*s + t]*m[t];
        }
        out[s] = v;
    }
}

inline void transmit_up(const double *T, int S, const Msg &a, Msg &out) {
    if (S == 2) {
        out[0] = a[0] + T[2]*a[1];
        out[1] = a[1] + T[1]*a[0];
        return;
    }
    for (int s = 0; s < S; s++) {
        double v = 0;
        for (int t = 0; t < S; t++) {
            v += T[S*t + s]*a[t];
        }
        out[s] = v;
    }
}

struct Genotype_likelihood {
    double pos;
    Node *leaves[2];
    double L[4];
};

struct Bin_gaps {
    bool built = false;
    vector<int> sample = {};
    vector<int> d = {};
    vector<double> before = {};
    double tail = 0;
    bool uniform = false;
    vector<int> full = {};
    double count = 0;
    vector<array<double, 4>> bases = {};
    array<double, 4> tail_bases = {0, 0, 0, 0};
    array<double, 4> all_bases = {0, 0, 0, 0};
    array<double, 4> count_bases = {0, 0, 0, 0};
};

struct Tree_order {
    vector<double> times = {};
    vector<int> order = {};
    vector<int> up = {};
    vector<vector<int>> kids = {};
    int top = -1;
    vector<int> leaf = {};

    void assign(const Flat_tree &flat, Node *root_node);

    int index_of(const Flat_tree &flat, Node *n) const;
};

struct Edge_transitions {
    int stamp = -1;
    int model = -1;
    double unit = -1;
    vector<double> P = {};
};

struct Peel {
    vector<char> missing = {};
    vector<Msg> msg = {};
    vector<Msg> empty = {};
    int count[4] = {0, 0, 0, 0};
    int leaf[2] = {-1, -1};
    vector<int> path = {};
    vector<char> in_path = {};
    vector<Msg> alt = {};
    vector<Msg> outside = {};
    vector<char> in_main = {};
    mutable vector<Edge_transitions> trans = {};
    mutable int trans_next = 0;
};

struct Site_messages {
    double unit = 0;
    vector<Msg> in = {};
    vector<Msg> above = {};
    Msg nu = {1.0, 1.0, 1.0, 1.0};
    int top = -1;
    int tree = -1;
    double empty_term = 0;
    double scale = 1;
    double z = 0;
    vector<Site_messages> parts = {};
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
    map<double, vector<pair<Branch, int>>> mutations = {};
    map<double, Recombination> recombinations = {};
    int bin_num = 0;
    double penalty = 1;
    double ancestral_prob = 0.5;
    Site_model model;
    vector<char> column_base = {};
    bool any_missing = false;
    vector<double> unassayed_sites = {};
    double sequence_length = 0;
    double bin_size = 0;
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
    
    void compute_rhos_thetas(double r, double m);
    
    void compute_rhos_thetas(double r, double m, Rate_map &rm, Rate_map &mm, double start);
    
    void build_singleton_arg(const Node_ptr &n);
    
    void add_sample(const Node_ptr &n);
    
    void add_node(Node *n);
    
    void add_new_node(double t, bool sample);

    void resort_after_time_change();
    
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

    void adjust_recombinations();

    void read_coordinates(string filename);

    void write_coordinates(string filename);

    void write(string node_file, string branch_file, string recomb_file, string mutation_file);

    void read(string node_file, string branch_file, string recomb_file);

    double get_arg_length();

    tuple<double, Branch, double> sample_internal_cut();

    tuple<double, Branch, double> sample_uniform_cut();

    vector<Genotype_likelihood> genotype_likelihoods = {};

    unordered_map<Node *, vector<int>> leaf_entries = {};

    const ARG *entry_owner = nullptr;

    void index_genotype_likelihoods();

    void index_missing(const vector<Node *> &leaves);

    vector<pair<double, int>> site_list = {};

    vector<pair<double, int>>::iterator site_at(double x);

    void set_collapsed(Node *leaf, double x, double y);

    const vector<Genotype_likelihood> &entries() const;

    pair<const int *, const int *> entries_in(Node *leaf, double x, double y) const;

    void site_missing(const Flat_tree &flat, double pos, vector<char> &missing);

    void root_prior(double *prior) const;

    Msg leaf_message(int code) const;

    void peel_site(const Flat_tree &flat, const Tree_order &t, double unit, double pos, Peel &p, int base = -1, const vector<char> *forced_missing = nullptr);

    double root_value(const Msg &top, bool seen0, bool seen1) const;

    double site_value(const Flat_tree &flat, const Tree_order &t, double unit, double pos, Peel &p);

    void emission_site(const Flat_tree &flat, const Tree_order &t, double unit, double pos, Node *query, Peel &p, Site_messages &s, int base = -1, const vector<char> *forced_missing = nullptr);

    void emission_single(const Flat_tree &flat, const Tree_order &t, double unit, double pos, Node *query, Peel &p, Site_messages &s, int base, const vector<char> *forced_missing);

    void draw_site_states(const Flat_tree &flat, const Tree_order &t, double unit, double pos, Peel &p, vector<int> &state);

    void mixture_counts(const Flat_tree &flat, const Tree_order &t, double unit, double pos, Peel &p, vector<Msg> &above, vector<double> &differs, vector<char> &informed_msg, vector<char> &informed_tip, int base = -1, const vector<char> *forced_missing = nullptr);

    void base_counts(double s, double t, array<double, 4> &n) const;

    void entry_weights(const Flat_tree &flat, const Tree_order &t, double unit, double pos, const Genotype_likelihood &g, Peel &p, double *w);

    void set_entry(const Flat_tree &flat, const Tree_order &t, double unit, double pos, const Genotype_likelihood &g, int c, Peel &p);

    double collapsed_site_value(const Flat_tree &flat, const Tree_order &t, double unit, double pos, const Genotype_likelihood &g, Peel &p);

    void draw_states(const Flat_tree &flat, const Tree_order &t, double unit, const Peel &p, vector<int> &state);

    void resample_site(const Flat_tree &flat, const Tree_order &t, int bin, double pos, size_t e, size_t f, Node *only, Peel &p);

    void resample_entries(Node *leaf, double x, double y);

    bool journal_on = false;

    vector<tuple<Node *, double, double>> state_journal = {};

    void undo_journal();

    void resample_genotypes();

    void draw_mutations();

    map<double, vector<char>> mutation_informed = {};

    void set_lengths(Flat_tree &flat, const function<double(Node *)> *time_of);

    void outside_messages(const Flat_tree &flat, const Tree_order &t, double unit, const Peel &p, vector<Msg> &above);

    void site_counts(const Flat_tree &flat, const Tree_order &t, double unit, const Peel &p, const vector<Msg> &above, vector<double> &differs, vector<char> &informed_msg, vector<char> &informed_tip);

    void draw_mutations(const function<double(Node *)> *time_of);

    vector<tuple<int, int, double, double>> branch_rows();

    void informed_edges(const function<void(Node *, Node *, double, double)> &add);

    void expected_counts(const function<double(Node *)> *time_of, const function<void(Node *, Node *, double)> &add, vector<double> *expected, const vector<tuple<int, int, double, double>> *rows);

    void write_counts(string filename, const vector<double> &expected, const vector<double> &informed_bases);

    void mapping_counts(const function<double(Node *)> *time_of, vector<double> &expected, vector<double> &informed_bases);

    vector<double> assayed = {};
    vector<pair<double, double>> masked = {};

    void compute_assayed();

    int num_unmapped();

    void check_incompatibility();

    int count_flipping();

    void clear_remove_info();

    double corrected_smc_prior(map<double, Branch> &lineage, int lo, int hi, double p, const Tree &tree_p);

    Flat_tree likelihood_tree;

    Tree_order likelihood_order;

    Peel likelihood_peel;

    vector<vector<Msg>> reference_msg = {}, pattern_msg = {};

    void reference_peel(const Flat_tree &flat, const Tree_order &t, double unit, vector<Msg> &msg, int base = 0) const;

    void repeel_path(const Flat_tree &flat, const Tree_order &t, double unit, int k, vector<Msg> &msg) const;

    double mutation_log_likelihood(map<double, Branch> &lineage, double x, double y, double p, const Tree &tree_p);

    mutable vector<Bin_gaps> bin_gaps = {};

    const Bin_gaps &gaps_of(int bin) const;

    Branch lineage_branch_before(map<double, Branch> &lineage, double x);

    set<double> get_check_points();

    bool check_disjoint_nodes(double x, double y);

    void new_recombination(double pos, Branch prev_added_branch, Branch prev_joining_branch, Branch next_added_branch, Branch next_joining_branch);

    void remove_empty_recombinations();

    void create_node_set();
    
    void write_nodes(string filename);
    
    void write_branches(string filename);
    
    void write_recombs(string filename);
    
    void write_mutations(string filename);
    
    void read_nodes(string filename);
    
    void read_branches(string filename);
    
    void read_recombs(string filename);
    
    void read_muts(string filename);

    void write_cells(string filename);

    void read_cells(string filename);

};

bool compare_edge(const tuple<int, int, double, double>& edge1, const tuple<int, int, double, double>& edge2);


#endif /* ARG_hpp */

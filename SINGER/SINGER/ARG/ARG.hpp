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

struct Peel {
    vector<char> missing = {};
    vector<array<double, 2>> msg = {};
    vector<array<double, 2>> empty = {};
    int count[2] = {0, 0};
    int leaf[2] = {-1, -1};
    vector<int> path = {};
    vector<char> in_path = {};
    vector<array<double, 2>> alt = {};
    vector<array<double, 2>> outside = {};
    vector<char> in_main = {};
};

struct Site_messages {
    double unit = 0;
    vector<array<double, 2>> in = {};
    vector<array<double, 2>> above = {};
    array<double, 2> nu = {1.0, 1.0};
    int top = -1;
    int tree = -1;
    double empty_term = 0;
    double scale = 1;
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
    double penalty = 0.01;
    double ancestral_prob = 0.5;
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

    void read(string node_file, string branch_file, string recomb_file, string mut_file);

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

    void peel_site(const Flat_tree &flat, const Tree_order &t, double unit, double pos, Peel &p);

    double root_value(const array<double, 2> &top, bool seen0, bool seen1) const;

    double site_value(const Flat_tree &flat, const Tree_order &t, double unit, double pos, Peel &p);

    void emission_site(const Flat_tree &flat, const Tree_order &t, double unit, double pos, Node *query, Peel &p, Site_messages &s);

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

    vector<array<double, 2>> reference_msg = {}, pattern_msg = {};

    void reference_peel(const Flat_tree &flat, const Tree_order &t, double unit, vector<array<double, 2>> &msg) const;

    void repeel_path(const Flat_tree &flat, const Tree_order &t, double unit, int k, vector<array<double, 2>> &msg) const;

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

};

bool compare_edge(const tuple<int, int, double, double>& edge1, const tuple<int, int, double, double>& edge2);


#endif /* ARG_hpp */

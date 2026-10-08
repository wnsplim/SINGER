//
//  Sampler.hpp
//  SINGER
//
//  Created by Yun Deng on 3/31/23.
//

#ifndef Sampler_hpp
#define Sampler_hpp

#include <stdio.h>
#include <chrono>
#include <sstream>
#include "ARG.hpp"
#include "Threader_smc.hpp"
#include "Binary_emission.hpp"
#include "Emission.hpp"
#include "Scaler.hpp"
#include "Rate_map.hpp"
#include "Data_reader.hpp"

class Sampler {

public:

    double rho_unit = 4e-3;
    double Ne = 1;
    Rate_map recomb_map;
    Rate_map mut_map;
    double mut_rate = 0;
    double recomb_rate = 0;
    string output_prefix = "";
    double start = 0;
    double end = 0;
    double sequence_length = 0;
    ARG arg;
    bool exact = false;
    int ploidy = 2;
    vector<double> tip_times = {};
    double tip_offset = 0;
    double bsp_c = 0.01;
    double tsp_q = 0.05;
    int random_seed = 0;
    double penalty = 0.01;
    double polar = 0.99;
    int scaling_rep = 5;
    int scaling_bin = 100;
    int sample_index = 0;
    set<Node_ptr, compare_node> sample_nodes = {};
    vector<Node_ptr> ordered_sample_nodes = {};
    Data_reader data;

    Sampler();

    Sampler(double pop_size, double r, double m);

    void set_precision(double c, double q);

    void set_input_file_prefix(string f);

    void set_output_file_prefix(string f);

    Node_ptr new_sample(int i);

    void read_tip_ages(string filename, double g);

    void order_samples();

    void load_data(double start_pos, double end_pos);

    void set_missing(vector<Node *> &leaves);

    void build_singleton_arg();

    void iterative_start();

    void sweep_samples();

    void internal_sample(int num_iters, int spacing);

    void resume_internal_sample(int num_iters, int spacing);

    void debug_resume_internal_sample(int num_iters, int spacing);

    void rescale();

    void start_log();

    void write_iterative_start();

    void write_sample();

    void load_resume_arg();

    void read_resume_point(string filename);

    void retract_log(int k);
};

#endif /* Sampler_hpp */

//
//  Threader_smc.hpp
//  SINGER
//
//  Created by Yun Deng on 4/4/23.
//

#ifndef Threader_smc_hpp
#define Threader_smc_hpp

#include <stdio.h>
#include <chrono>
#include <sstream>
#include "ARG.hpp"
#include "Polar_emission.hpp"
#include "approx_BSP.hpp"
#include "TSP.hpp"

class Threader_smc {

public:

    Threader_smc(double c, double q);
    
    ~Threader_smc();
    
    void reset();

    void thread(ARG &a, Node_ptr n);
    
    void internal_rethread(ARG &a, tuple<double, Branch, double> cut_point);

    void exact_internal_rethread(ARG &a, tuple<double, Branch, double> cut_point);

    double cut_time = 0;
    double start = 0;
    double end = 0;
    int start_index = 0;
    int end_index = 0;
    approx_BSP bsp = approx_BSP();
    TSP tsp = TSP();
    double gap;
    double cutoff;
    static bool no_data;
    shared_ptr<Binary_emission> be = make_shared<Binary_emission>();
    shared_ptr<Polar_emission> pe = make_shared<Polar_emission>();
    map<double, Branch> new_joining_branches = {};
    map<double, Branch> added_branches = {};
    double pos_lo = 0;
    Tree tree_lo;
    ARG new_arg;
    ARG old_arg;

    void get_boundary(ARG &a);
    
    void set_check_points(ARG &a);

    void run_BSP(ARG &a);

    void run_TSP(ARG &a);

    void run_TSP(ARG &a, map<double, Branch> &jb);
    
    void sample_joining_branches(ARG &a);

    void sample_joining_points(ARG &a);
    
    double acceptance_ratio(ARG &a);

    bool has_bridges(ARG &a);

    bool bridges_kept(ARG &a);

    double cut_ratio(ARG &a);

    double exact_acceptance_ratio(ARG &a);
    
    double random();

};

#endif /* Threader_smc_hpp */

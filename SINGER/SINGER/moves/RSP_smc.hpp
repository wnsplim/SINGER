//
//  RSP_smc.hpp
//  SINGER
//
//  Created by Yun Deng on 3/13/23.
//

#ifndef RSP_smc_hpp
#define RSP_smc_hpp

#include <stdio.h>
#include <math.h>
#include <map>
#include "random_utils.hpp"
#include "Branch.hpp"
#include "Tree.hpp"
#include "Recombination.hpp"

class RSP_smc {

public:

    RSP_smc();

    void set_tree(Tree &tree);

    void sample_recombination(Recombination &r, double cut_time, Tree &tree, const Branch &own);

    double log_start_density(Recombination &r);

    double log_start_marginal(Recombination &r, double cut_time, const Branch &own);

    double unchanged_recomb_length(Tree &tree);

    void approx_sample_recombination(Recombination &r, double cut_time);

    void adjust(Recombination &r, double cut_time);

    void approx_sample_recombination(Recombination &r, double cut_time, double n);

    void adjust(Recombination &r, double cut_time, double n);

private:

    vector<double> level_times = {};
    vector<double> level_rates = {};
    vector<double> level_lambda = {};

    int level_of(double s);

    double lambda(double s);

    double exp_lambda_integral(double lo, double hi);

    double draw_start_time(double lo, double hi);

    vector<Branch> source_candidates(Recombination &r);

    double start_lower_bound(const Branch &candidate, double cut_time, const Branch &own);

    double random_time(double lb, double ub);

    double random_time(double lb, double ub, double q);

    double choose_time(double lb, double ub);

    double choose_time(double lb, double ub, double n);

};

#endif /* RSP_smc_hpp */

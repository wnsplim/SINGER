//
//  approx_coalescent_calculator.hpp
//  SINGER
//
//  Created by Yun Deng on 6/25/23.
//

#ifndef approx_coalescent_calculator_hpp
#define approx_coalescent_calculator_hpp

#include <stdio.h>
#include <map>
#include <math.h>
#include "Branch.hpp"
#include "Tree.hpp"
#include "Recombination.hpp"

class approx_coalescent_calculator {

public:

    double cut_time;
    double rho = 0;
    double first_moment = 0;

    approx_coalescent_calculator(double t);

    ~approx_coalescent_calculator();

    void start(set<Branch> &branches);

    void start(Tree &tree);

    void update(Recombination &r);

    pair<double, double> compute_time_weights(double x, double y);

    void compute_first_moment();

    double prob(double x, double y);

    double find_median(double x, double y);

private:

    vector<double> t = {};
    vector<double> Lam = {};
    vector<double> G = {};
    vector<double> Q = {};
    double tail_G = 0;
    double tail_Q = 0;
    int rebuild_from = 0;

    void refresh();

    void at(double x, double &g, double &q);

};

#endif /* approx_coalescent_calculator_hpp */

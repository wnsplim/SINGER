//
//  coalescent_calculator.hpp
//  SINGER
//
//  Created by Yun Deng on 6/25/23.
//

#ifndef coalescent_calculator_hpp
#define coalescent_calculator_hpp

#include <stdio.h>
#include <map>
#include <math.h>
#include "Branch.hpp"
#include "Tree.hpp"
#include "Recombination.hpp"

class coalescent_calculator {

public:

    double cut_time;
    int extra = 0;
    double rho = 0;
    double first_moment = 0;

    coalescent_calculator(double t);

    ~coalescent_calculator();

    void start(set<Branch> &branches);

    void start(Tree &tree);

    void update(Recombination &r);

    pair<double, double> compute_time_weights(double x, double y);

    double prob(double x, double y);

    double find_median(double x, double y);

    double surv(double x);

    double rate(double x);

    double surv_inv(double p);

    double recomb_mass(double s, double t);

    size_t num_times() const { return t.size(); }

private:

    vector<double> t = {};
    vector<double> Lam = {};
    vector<double> E = {};
    vector<double> G = {};
    vector<double> Q = {};
    vector<double> W = {};
    vector<double> B = {};
    double tail_G = 0;
    double tail_Q = 0;
    int rebuild_from = 0;

    void refresh();

    void at(double x, double &g, double &q);

};

#endif /* coalescent_calculator_hpp */

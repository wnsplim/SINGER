#ifndef approx_coalescent_calculator_hpp
#define approx_coalescent_calculator_hpp

#include <stdio.h>
#include <math.h>
#include "Branch.hpp"
#include "Recombination.hpp"

class approx_coalescent_calculator {

public:

    double cut_time;
    int n0 = 0;

    approx_coalescent_calculator(double t);

    void start(set<Branch> &branches);

    void update(Recombination &r);

    double prob(double x, double y);

    double prob_integral(double x);

    double find_median(double x, double y);

};

#endif

//
//  Interval.hpp
//  SINGER
//
//  Created by Yun Deng on 4/9/22.
//

#ifndef Interval_hpp
#define Interval_hpp

#include <stdio.h>
#include <map>
#include <numeric>
#include "Recombination.hpp"

using namespace std;

class Interval {
    
public:
    
    Branch branch = Branch();
    double lb = 0;
    double ub = 0;
    double weight = 0.0;
    double time = 0.0;
    int start_pos = 0;
    int source_pos = 0;
    Node *node = nullptr;
    double s_lb = -1.0;
    double s_ub = -1.0;
    Interval *source_interval = nullptr;

    vector<double> source_weights = {};
    vector<Interval *> source_intervals = {};

    Interval(Branch b, double tl, double tu, int init_pos);

    void fill_time();

    bool full(double t);
};

class Interval_info {
    
public:
    
    Branch branch;
    double lb = 0;
    double ub = 0;
    double time = 0;
    double seed_pos = 0;
    
    Interval_info();
    
    Interval_info(Branch b, double tl, double tu);

    bool operator<(const Interval_info& other) const;
    
};

#endif /* Interval_hpp */

//
//  Binary_emission.hpp
//  SINGER
//
//  Created by Yun Deng on 4/6/23.
//

#ifndef Binary_emission_hpp
#define Binary_emission_hpp

#include <stdio.h>
#include <math.h>
#include "Emission.hpp"

using namespace std;

class Binary_emission : public Emission {
    
public:
    
    map<double, double> num_unmapped = {};
    double penalty = 0.01;
    double ancestral_prob = 0.5;
    
    vector<double> diff = vector<double>(4);
    
    Binary_emission();
    
    ~Binary_emission();
    
    double null_emit(Branch &branch, double time, double theta, Node *node) override;
    
    double mut_emit(Branch &branch, double time, double theta, double bin_size, vector<double> &mut_set, Node *node) override;
    
    double emit(Branch &branch, double time, double theta, double bin_size, vector<double> &emissions, Node *node) override;
    
    double calculate_prob(double theta, double bin_size, double ll, double lu, double l0, int sl, int su, int s0);
    
    double calculate_prob(double theta, double bin_size, int s);
    
    void get_diff(vector<double> &mut_set, const Branch &branch, Node *node);
};

#endif /* Binary_emission_hpp */

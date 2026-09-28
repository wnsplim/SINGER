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
    
    double penalty = 0.01;
    double ancestral_prob = 0.5;

    Binary_emission();

    ~Binary_emission();

    double null_emit(Branch &branch, double time, double theta, Node *node) override;

    double mut_emit(Branch &branch, double time, double theta, double bin_size, vector<double> &mut_set, Node *node) override;
};

#endif /* Binary_emission_hpp */

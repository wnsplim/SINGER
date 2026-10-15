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

    double mut_emit(Branch &branch, double time, double theta, double bin_size, vector<pair<double, int>> &mut_set, Node *node, double norm) override;
};

#endif /* Binary_emission_hpp */

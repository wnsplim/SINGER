//
//  Polar_emission.hpp
//  SINGER
//
//  Created by Yun Deng on 6/14/23.
//

#ifndef Polar_emission_hpp
#define Polar_emission_hpp

#include <stdio.h>
#include <math.h>
#include "Emission.hpp"

using namespace std;

class Polar_emission : public Emission {
    
public:
    
    double penalty = 0.01;
    double ancestral_prob = 0.5;

    Polar_emission();

    ~Polar_emission();

    double null_emit(Branch &branch, double time, double theta, Node *node) override;

    double mut_emit(Branch &branch, double time, double theta, double bin_size, vector<double> &mut_set, Node *node) override;
};


#endif /* Polar_emission_hpp */

//
//  Rate_map.cpp
//  SINGER
//
//  Created by Yun Deng on 6/4/24.
//

#include "Rate_map.hpp"

Rate_map::Rate_map() {}

void Rate_map::load_map(string mut_map_file, double start, double end) {
    ifstream fin(mut_map_file);
    if (!fin.good()) {
        cerr << "input rate map file not found" << endl;
        exit(1);
    }
    rate_distances.push_back(0);
    double left;
    double right;
    double rate;
    double mut_dist;
    while (fin >> left >> right >> rate) {
        if (coordinates.size() > 0 and left != sequence_length) {
            cerr << "Error: each segment of the rate map " << mut_map_file << " must start where the previous one ends. " << endl;
            exit(1);
        }
        coordinates.push_back(left);
        mut_dist = rate_distances.back() + rate*(right - left);
        rate_distances.push_back(mut_dist);
        sequence_length = right;
    }
    if (!fin.eof() or coordinates.empty()) {
        cerr << "Error: every line of the rate map " << mut_map_file << " must be 'left right rate'. " << endl;
        exit(1);
    }
    coordinates.push_back(sequence_length);
    if (coordinates.front() > start or sequence_length < end) {
        cerr << "Error: the rate map " << mut_map_file << " must cover -start to -end. " << endl;
        exit(1);
    }
}

int Rate_map::find_index(double x) {
    auto it = upper_bound(coordinates.begin(), coordinates.end(), x);
    it--;
    int index = (int) distance(coordinates.begin(), it);
    assert(index >= 0 and index <= coordinates.size() - 1);
    return index;
}

double Rate_map::cumulative_distance(double x) {
    int index = find_index(x);
    if (index + 1 == (int) coordinates.size()) {
        return rate_distances[index];
    }
    double prev_dist = rate_distances[index];
    double next_dist = rate_distances[index+1];
    double p = (x - coordinates[index])/(coordinates[index+1] - coordinates[index]);
    double dist = (1-p)*prev_dist + p*next_dist;
    return dist;
}

double Rate_map::segment_distance(double x, double y) {
    int i = find_index(x);
    if (i + 1 < (int) coordinates.size() and y <= coordinates[i+1]) {
        // a product, so equal-width bins inside one map segment get bit-equal values (the HMMs cache on rho)
        return (rate_distances[i+1] - rate_distances[i])/(coordinates[i+1] - coordinates[i])*(y - x);
    }
    return cumulative_distance(y) - cumulative_distance(x);
}

double Rate_map::mean_rate(double x, double y) {
    return segment_distance(x, y)/(y - x);
}

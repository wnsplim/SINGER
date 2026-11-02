//
//  Site_model.hpp
//  SINGER
//

#ifndef Site_model_hpp
#define Site_model_hpp

#include <string>
#include <vector>
using namespace std;

class Site_model {

public:

    int states = 2;
    string name = "CFN";
    vector<double> pi = {0.5, 0.5};
    vector<double> rates = {1.0};
    double alpha = 0;
    int categories = 1;
    double kappa = 1;
    vector<double> exchangeabilities = vector<double>(6, 1.0);
    vector<double> rate_matrix = {};
    int stamp = 0;

    void set(const string &model, const vector<double> &frequencies, const vector<double> &exchangeabilities, double kappa);

    void set_gamma(int k, double shape);

    void refresh() {
        set(name, pi, exchangeabilities, kappa);
        set_gamma(categories, alpha);
    }

    bool has_kappa() const { return name == "K80" or name == "HKY"; }

    bool has_frequencies() const { return name == "F81" or name == "HKY" or name == "GTR"; }

    bool has_exchangeabilities() const { return name == "GTR"; }

    void transition(double x, double *T) const {
        if (states == 2) {
            T[0] = 1;
            T[1] = x;
            T[2] = x;
            T[3] = 1;
            return;
        }
        transition_eigen(x, T);
    }

    void transition_eigen(double x, double *T) const;

    void count_matrix(double x, double *M) const;

    double row_sum(double x) const;

    static double regularized_gamma(double a, double x);

private:

    vector<double> eigenvalues = {};
    vector<double> left = {};
    vector<double> right = {};
    vector<double> coupling = {};
    vector<double> cijk = {};
    vector<double> count_tensor = {};
    vector<double> distinct = {};
    vector<double> cijg = {};

    void decompose(const vector<double> &Q);
};

#endif /* Site_model_hpp */

//
//  Site_model.hpp
//  SINGER
//

#ifndef Site_model_hpp
#define Site_model_hpp

#include <string>
#include <vector>
#include <utility>
using namespace std;

const int max_states = 4;

class Site_model {

public:

    int states = 2;
    string name = "CFN";
    // true when CFN runs the landed chain: transition [1, x; x, 1], at most one change at a uniform time, the HMM's branch-product emission
    bool branch_product = true;
    vector<double> pi = {0.5, 0.5};
    vector<double> rates = {1.0};
    double alpha = 0;
    int categories = 1;
    double kappa = 1;
    vector<double> exchangeabilities = {1.0};
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

    bool continuous_gamma() const { return categories == 0; }

    void transition(double x, double *T) const {
        if (branch_product) {
            T[0] = 1;
            T[1] = x;
            T[2] = x;
            T[3] = 1;
            return;
        }
        transition_eigen(x, T);
    }

    void transition_eigen(double x, double *T) const;

    double transition_entry(double x, int i, int j) const;

    double transition_dot(double x, int i, const double *w) const;

    double rate(int i, int j) const { return rate_matrix[states*i + j]; }

    double rate_integrated_entry(int b, double A, double l) const;

    long double log_rate_integrated_stay(int b, long double A, long double l) const;

    double h(int b, double x) const;

    vector<double> h_max = {};

    double sample_change_rate(int b, double A, double l) const;

    long double node_expectation(int b, double A, double l1, double l2, double l3, double l0) const;

    double node_excess(int b, double x1, double x2, double x3, double x0) const;

    double node_rate(int b, double A, double l1, double l2, double l3, double l0) const;

    double node_sum(double x1, double x2, double x3, int xu, int xl, const double *w) const;

    void count_matrix(double x, double *M) const;

    void sample_path(int a, int b, double x, double P_ab, vector<pair<double, int>> &path, bool change = false) const;

    static double regularized_gamma(double a, double x);

private:

    vector<double> rate_matrix = {};
    vector<double> eigenvalues = {};
    vector<double> left = {};
    vector<double> right = {};
    vector<double> coupling = {};
    vector<double> cijk = {};
    vector<double> count_tensor = {};
    vector<double> distinct = {};
    vector<double> cijg = {};
    vector<double> split_terms = {};

    void decompose(const vector<double> &Q);

    void find_h_max();
};

#endif /* Site_model_hpp */

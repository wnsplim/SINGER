//
//  Site_model.cpp
//  SINGER
//

#include "Site_model.hpp"
#include <cmath>
#include <iostream>
#include <limits>

void Site_model::set(const string &model, const vector<double> &frequencies, const vector<double> &exchangeabilities_in, double kappa_in) {
    static int stamps = 0;
    stamp = ++stamps;
    name = model;
    kappa = kappa_in;
    exchangeabilities = exchangeabilities_in.empty() ? vector<double>(6, 1.0) : exchangeabilities_in;
    if (model == "CFN") {
        states = 2;
        pi = {0.5, 0.5};
        return;
    }
    states = 4;
    pi = (model == "JC" or model == "K80") ? vector<double>{0.25, 0.25, 0.25, 0.25} : frequencies;
    vector<double> s(6, 1.0);
    if (model == "K80" or model == "HKY") {
        s[1] = s[4] = kappa;
    } else if (model == "GTR") {
        s = exchangeabilities;
    }
    vector<double> Q(16, 0.0);
    int pair_index[4][4] = {{-1, 0, 1, 2}, {0, -1, 3, 4}, {1, 3, -1, 5}, {2, 4, 5, -1}};
    double total = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i != j) {
                Q[4*i + j] = s[pair_index[i][j]]*pi[j];
                Q[4*i + i] -= Q[4*i + j];
            }
        }
        total -= pi[i]*Q[4*i + i];
    }
    for (double &q : Q) {
        q /= total;
    }
    rate_matrix = Q;
    decompose(Q);
}

void Site_model::count_matrix(double x, double *M) const {
    if (states == 2) {
        M[0] = 0;
        M[1] = x;
        M[2] = x;
        M[3] = 0;
        return;
    }
    int n = states;
    double e[4], J[16];
    for (int m = 0; m < n; m++) {
        e[m] = exp(eigenvalues[m]*x);
    }
    for (int m = 0; m < n; m++) {
        for (int q = 0; q < n; q++) {
            double y = (eigenvalues[m] - eigenvalues[q])*x;
            double f = fabs(y) < 1e-8 ? 1 + y/2 : expm1(y)/y;
            J[n*m + q] = e[q]*x*f;
        }
    }
    const double *d = count_tensor.data();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double v = 0;
            for (int k = 0; k < n*n; k++) {
                v += *(d++)*J[k];
            }
            M[n*i + j] = max(v, 0.0);
        }
    }
}

void Site_model::decompose(const vector<double> &Q) {
    int n = states;
    vector<double> B(n*n), W(n*n, 0.0);
    for (int i = 0; i < n; i++) {
        W[n*i + i] = 1.0;
        for (int j = 0; j < n; j++) {
            B[n*i + j] = sqrt(pi[i])*Q[n*i + j]/sqrt(pi[j]);
        }
    }
    int sweep = 0;
    for (; sweep < 100; sweep++) {
        double off = 0;
        for (int p = 0; p < n; p++) {
            for (int q = p + 1; q < n; q++) {
                off += B[n*p + q]*B[n*p + q];
            }
        }
        if (off < 1e-32) {
            break;
        }
        for (int p = 0; p < n; p++) {
            for (int q = p + 1; q < n; q++) {
                if (fabs(B[n*p + q]) < 1e-300) {
                    continue;
                }
                double theta = (B[n*q + q] - B[n*p + p])/(2*B[n*p + q]);
                double t = (theta >= 0 ? 1.0 : -1.0)/(fabs(theta) + sqrt(theta*theta + 1));
                double c = 1/sqrt(t*t + 1), s = t*c;
                for (int k = 0; k < n; k++) {
                    double bkp = B[n*k + p], bkq = B[n*k + q];
                    B[n*k + p] = c*bkp - s*bkq;
                    B[n*k + q] = s*bkp + c*bkq;
                }
                for (int k = 0; k < n; k++) {
                    double bpk = B[n*p + k], bqk = B[n*q + k];
                    B[n*p + k] = c*bpk - s*bqk;
                    B[n*q + k] = s*bpk + c*bqk;
                }
                for (int k = 0; k < n; k++) {
                    double wkp = W[n*k + p], wkq = W[n*k + q];
                    W[n*k + p] = c*wkp - s*wkq;
                    W[n*k + q] = s*wkp + c*wkq;
                }
            }
        }
    }
    if (sweep == 100) {
        cerr << "Error: the eigen-decomposition of the " << name << " rate matrix did not converge." << endl;
        exit(1);
    }
    eigenvalues.resize(n);
    left.resize(n*n);
    right.resize(n*n);
    for (int m = 0; m < n; m++) {
        eigenvalues[m] = B[n*m + m];
        for (int i = 0; i < n; i++) {
            left[n*i + m] = W[n*i + m]/sqrt(pi[i]);
            right[n*m + i] = W[n*i + m]*sqrt(pi[i]);
        }
    }
    coupling.assign(n*n, 0.0);
    for (int m = 0; m < n; m++) {
        for (int q = 0; q < n; q++) {
            double c = 0;
            for (int a = 0; a < n; a++) {
                for (int b = 0; b < n; b++) {
                    if (a != b) {
                        c += right[n*m + a]*Q[n*a + b]*left[n*b + q];
                    }
                }
            }
            coupling[n*m + q] = c;
        }
    }
    cijk.assign(n*n*n, 0.0);
    count_tensor.assign(n*n*n*n, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int m = 0; m < n; m++) {
                cijk[n*n*i + n*j + m] = left[n*i + m]*right[n*m + j];
                for (int q = 0; q < n; q++) {
                    count_tensor[n*n*n*i + n*n*j + n*m + q] = coupling[n*m + q]*left[n*i + m]*right[n*q + j];
                }
            }
        }
    }
    double scale = 0;
    for (int m = 0; m < n; m++) {
        scale = max(scale, fabs(eigenvalues[m]));
    }
    distinct.clear();
    vector<int> group(n, -1);
    for (int m = 0; m < n; m++) {
        for (int g = 0; g < (int) distinct.size(); g++) {
            if (fabs(eigenvalues[m] - distinct[g]) <= 1e-12*scale) {
                group[m] = g;
                break;
            }
        }
        if (group[m] < 0) {
            distinct.push_back(eigenvalues[m]);
            group[m] = (int) distinct.size() - 1;
        }
    }
    int G = (int) distinct.size();
    cijg.assign(n*n*G, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int m = 0; m < n; m++) {
                cijg[n*G*i + G*j + group[m]] += cijk[n*n*i + n*j + m];
            }
        }
    }
}

void Site_model::transition_eigen(double x, double *T) const {
    int n = states, G = (int) distinct.size();
    double e[4];
    for (int g = 0; g < G; g++) {
        e[g] = exp(distinct[g]*x);
    }
    const double *c = cijg.data();
    for (int i = 0; i < n; i++) {
        double row = 0;
        for (int j = 0; j < n; j++) {
            double v = 0;
            for (int g = 0; g < G; g++) {
                v += *(c++)*e[g];
            }
            v = max(v, 0.0);
            T[n*i + j] = v;
            row += v;
        }
        for (int j = 0; j < n; j++) {
            T[n*i + j] /= row;
        }
    }
}

double Site_model::row_sum(double x) const {
    return states == 2 ? 1 + x : 1.0;
}

double Site_model::regularized_gamma(double a, double x) {
    if (x <= 0) {
        return 0;
    }
    double lg = lgamma(a);
    if (x < a + 1) {
        double term = 1/a, sum = term;
        for (int k = 1; k < 10000; k++) {
            term *= x/(a + k);
            sum += term;
            if (fabs(term) < fabs(sum)*1e-17) {
                return sum*exp(-x + a*log(x) - lg);
            }
        }
    } else {
        double tiny = numeric_limits<double>::min();
        double b = x + 1 - a, c = 1/tiny, d = 1/b, h = d;
        for (int k = 1; k < 10000; k++) {
            double an = -k*(k - a);
            b += 2;
            d = an*d + b;
            if (fabs(d) < tiny) d = tiny;
            c = b + an/c;
            if (fabs(c) < tiny) c = tiny;
            d = 1/d;
            double delta = d*c;
            h *= delta;
            if (fabs(delta - 1) < 1e-16) {
                return 1 - exp(-x + a*log(x) - lg)*h;
            }
        }
    }
    cerr << "Error: the incomplete gamma function did not converge at shape " << a << "." << endl;
    exit(1);
}

void Site_model::set_gamma(int k, double shape) {
    alpha = shape;
    categories = k;
    rates.assign(1, 1.0);
    if (k == 1) {
        return;
    }
    vector<double> cuts(k + 1, 0.0);
    cuts[k] = numeric_limits<double>::infinity();
    for (int i = 1; i < k; i++) {
        double target = (double) i/k, lo = 0, hi = 1;
        while (regularized_gamma(shape, shape*hi) < target) {
            hi *= 2;
        }
        for (int it = 0; it < 200; it++) {
            double mid = (lo + hi)/2;
            if (regularized_gamma(shape, shape*mid) < target) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        cuts[i] = (lo + hi)/2;
    }
    rates.resize(k);
    for (int i = 0; i < k; i++) {
        double upper = i + 1 == k ? 1.0 : regularized_gamma(shape + 1, shape*cuts[i + 1]);
        double lower = i == 0 ? 0.0 : regularized_gamma(shape + 1, shape*cuts[i]);
        rates[i] = k*(upper - lower);
    }
}

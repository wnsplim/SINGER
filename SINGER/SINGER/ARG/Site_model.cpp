//
//  Site_model.cpp
//  SINGER
//

#include "Site_model.hpp"
#include "random_utils.hpp"
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
    find_h_max();
}

static double stay_phi2(double lambda, double q, double x) {
    double y = (lambda + q)*x;
    if (fabs(y) < 1e-2) {
        return exp(-q*x)*(1.0/2 + y*(1.0/6 + y*(1.0/24 + y*(1.0/120 + y*(1.0/720 + y/5040)))));
    }
    if (y > 1) {
        return (exp(lambda*x) - exp(-q*x)*(1 + y))/(y*y);
    }
    return exp(-q*x)*(expm1(y) - y)/(y*y);
}

double Site_model::h(int b, double x) const {
    int G = (int) distinct.size();
    const double *c = cijg.data() + states*G*b + G*b;
    double q = -rate_matrix[states*b + b], v = 0;
    for (int g = 0; g < G; g++) {
        double z = distinct[g] + q;
        v += c[g]*z*z*stay_phi2(distinct[g], q, x);
    }
    return v;
}

void Site_model::find_h_max() {
    double top = 0;
    for (double l : distinct) {
        top = max(top, fabs(l));
    }
    h_max.assign(states, 0.0);
    for (int b = 0; b < states; b++) {
        double q = -rate_matrix[states*b + b], h0 = h(b, 0.0);
        double end = max(2/sqrt(h0), 40/top), step = 0.05/(top + q);
        int n = (int) ceil(end/step);
        int best = 0;
        double best_v = h0;
        for (int i = 1; i <= n; i++) {
            double v = h(b, i*step);
            if (v > best_v) {
                best_v = v;
                best = i;
            }
        }
        double lo = max(0, best - 1)*step, hi = (best + 1)*step, r = (sqrt(5.0) - 1)/2;
        for (int it = 0; it < 100; it++) {
            double m1 = hi - r*(hi - lo), m2 = lo + r*(hi - lo);
            if (h(b, m1) < h(b, m2)) {
                lo = m1;
            } else {
                hi = m2;
            }
        }
        h_max[b] = max(best_v, h(b, (lo + hi)/2));
    }
}

double Site_model::transition_dot(double x, int i, const double *w) const {
    int G = (int) distinct.size();
    double e[4];
    for (int g = 0; g < G; g++) {
        e[g] = exp(distinct[g]*x);
    }
    double total = 0;
    for (int j = 0; j < states; j++) {
        if (w[j] != 0) {
            const double *c = cijg.data() + states*G*i + G*j;
            double v = 0;
            for (int g = 0; g < G; g++) {
                v += c[g]*e[g];
            }
            total += w[j]*max(v, 0.0);
        }
    }
    return total;
}

double Site_model::rate_integrated_entry(int b, double A, double l) const {
    int G = (int) distinct.size();
    const double *c = cijg.data() + states*G*b + G*b;
    double v = 0;
    for (int g = 0; g < G; g++) {
        v += c[g]*exp(-alpha*log1p(-distinct[g]*l/A));
    }
    return v;
}

long double Site_model::log_rate_integrated_stay(int b, long double A, long double l) const {
    return -(long double) alpha*log1pl(-(long double) rate_matrix[states*b + b]*l/A);
}

double Site_model::sample_change_rate(int b, double A, double l) const {
    double q = -rate_matrix[states*b + b];
    gamma_distribution<double> proposal(alpha + 1, 1/A);
    while (true) {
        double nu = proposal(random_engine);
        if (uniform_random() < -expm1(-q*nu*l)/(q*nu*l)) {
            return nu;
        }
    }
}

long double Site_model::node_expectation(int b, double A, double l1, double l2, double l3, double l0) const {
    int G = (int) distinct.size();
    const double *c = cijg.data() + states*G*b + G*b, *w = split_terms.data() + G*G*G*b;
    long double v = 0;
    if (!continuous_gamma()) {
        long double e1[4], e2[4], e3[4], p = 0;
        for (int g = 0; g < G; g++) {
            e1[g] = expl((long double) distinct[g]*l1);
            e2[g] = expl((long double) distinct[g]*l2);
            e3[g] = expl((long double) distinct[g]*l3);
            p += c[g]*expl((long double) distinct[g]*l0);
        }
        for (int i = 0; i < G; i++) {
            for (int j = 0; j < G; j++) {
                for (int k = 0; k < G; k++) {
                    v += w[G*G*i + G*j + k]*e1[i]*e2[j]*e3[k];
                }
            }
        }
        return v*p;
    }
    auto power = [&](double lambda) { return exp(-alpha*log1p(-lambda/A)); };
    double s = 0;
    for (int i = 0; i < G; i++) {
        for (int j = 0; j < G; j++) {
            for (int k = 0; k < G; k++) {
                double lambda = distinct[i]*l1 + distinct[j]*l2 + distinct[k]*l3;
                if (l0 == 0) {
                    s += w[G*G*i + G*j + k]*power(lambda);
                    continue;
                }
                for (int m = 0; m < G; m++) {
                    s += w[G*G*i + G*j + k]*c[m]*power(lambda + distinct[m]*l0);
                }
            }
        }
    }
    return s;
}

double Site_model::node_excess(int b, double x1, double x2, double x3, double x0) const {
    double x[4] = {x1, x2, x3, x0}, q = -rate_matrix[states*b + b], after[5], before = 1, v = 0, r = 0;
    after[4] = 1;
    for (int i = 3; i >= 0; i--) {
        after[i] = after[i + 1]*exp(-q*x[i]);
    }
    for (int i = 0; i < 4; i++) {
        v += x[i]*x[i]*h(b, x[i])*before*after[i + 1];
        before *= transition_entry(x[i], b, b);
    }
    for (int s = 0; s < states; s++) {
        if (s != b) {
            r += transition_entry(x1, b, s)*transition_entry(x2, s, b)*transition_entry(x3, s, b);
        }
    }
    return v + r*transition_entry(x0, b, b);
}

double Site_model::node_rate(int b, double A, double l1, double l2, double l3, double l0) const {
    double q = -rate_matrix[states*b + b], bound = h_max[b]*(l1*l1 + l2*l2 + l3*l3 + l0*l0);
    for (int s = 0; s < states; s++) {
        double qs = -rate_matrix[states*s + s];
        if (s != b) {
            bound += qs*min(min(q*l1*l2, q*l1*l3), qs*l2*l3);
        }
    }
    gamma_distribution<double> proposal(alpha + 2, 1/A);
    while (true) {
        double nu = proposal(random_engine);
        if (uniform_random()*bound*nu*nu < node_excess(b, nu*l1, nu*l2, nu*l3, nu*l0)) {
            return nu;
        }
    }
}

double Site_model::node_sum(double x1, double x2, double x3, int xu, int xl, const double *w) const {
    int G = (int) distinct.size();
    double e1[4], e2[4], e3[4], total = 0;
    for (int g = 0; g < G; g++) {
        e1[g] = exp(distinct[g]*x1);
        e2[g] = exp(distinct[g]*x2);
        e3[g] = exp(distinct[g]*x3);
    }
    for (int y = 0; y < states; y++) {
        const double *cu = cijg.data() + states*G*xu + G*y, *cl = cijg.data() + states*G*y + G*xl;
        double up = 0, low = 0, query = 0;
        for (int g = 0; g < G; g++) {
            up += cu[g]*e1[g];
            low += cl[g]*e2[g];
        }
        for (int a = 0; a < states; a++) {
            if (w[a] != 0) {
                const double *cq = cijg.data() + states*G*y + G*a;
                double v = 0;
                for (int g = 0; g < G; g++) {
                    v += cq[g]*e3[g];
                }
                query += w[a]*max(v, 0.0);
            }
        }
        total += max(up, 0.0)*max(low, 0.0)*query;
    }
    return total;
}

void Site_model::sample_path(int a, int b, double x, double P_ab, vector<pair<double, int>> &path, bool change) const {
    path.clear();
    if (states == 2) {
        if (a != b) {
            path.push_back({uniform_random()*x, b});
        }
        return;
    }
    int n = states;
    double Pbuf[16], J[4], p[4];
    double done = 0;
    int s = a;
    bool first = true;
    while (true) {
        double T = x - done;
        double q = -rate_matrix[n*s + s];
        double stay = exp(-q*T);
        if (s == b and !(first and change)) {
            double P_sb = P_ab;
            if (!first) {
                transition_eigen(T, Pbuf);
                P_sb = Pbuf[n*s + b];
            }
            if (uniform_random()*P_sb < stay) {
                return;
            }
        }
        first = false;
        for (int j = 0; j < n; j++) {
            double y = (eigenvalues[j] + q)*T;
            J[j] = stay*T*(fabs(y) < 1e-8 ? 1 + y/2 : expm1(y)/y);
        }
        double total = 0;
        for (int i = 0; i < n; i++) {
            p[i] = 0;
            if (i == s) {
                continue;
            }
            double c = 0;
            for (int j = 0; j < n; j++) {
                c += left[n*i + j]*right[n*j + b]*J[j];
            }
            p[i] = max(rate_matrix[n*s + i]*c, 0.0);
            total += p[i];
        }
        double u = uniform_random()*total, acc = 0;
        int next = -1;
        for (int i = 0; i < n; i++) {
            if (p[i] == 0) {
                continue;
            }
            acc += p[i];
            next = i;
            if (u < acc) {
                break;
            }
        }
        auto cdf = [&](double tau) {
            double c = 0;
            for (int j = 0; j < n; j++) {
                double y = -(eigenvalues[j] + q)*tau;
                c += left[n*next + j]*right[n*j + b]*exp(eigenvalues[j]*T)*tau*(fabs(y) < 1e-8 ? 1 + y/2 : expm1(y)/y);
            }
            return rate_matrix[n*s + next]*c;
        };
        double v = uniform_random()*p[next], lo = 0, hi = T;
        for (int it = 0; it < 60; it++) {
            double mid = 0.5*(lo + hi);
            if (cdf(mid) < v) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        double tau = 0.5*(lo + hi);
        path.push_back({done + tau, next});
        done += tau;
        s = next;
    }
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
                B[n*p + q] = B[n*q + p] = 0;
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
    split_terms.assign(n*G*G*G, 0.0);
    for (int b = 0; b < n; b++) {
        for (int s = 0; s < n; s++) {
            for (int i = 0; i < G; i++) {
                for (int j = 0; j < G; j++) {
                    for (int k = 0; k < G; k++) {
                        split_terms[G*G*G*b + G*G*i + G*j + k] += cijg[n*G*b + G*s + i]*cijg[n*G*s + G*b + j]*cijg[n*G*s + G*b + k];
                    }
                }
            }
        }
    }
}

double Site_model::transition_entry(double x, int i, int j) const {
    int G = (int) distinct.size();
    const double *c = cijg.data() + states*G*i + G*j;
    double v = 0;
    for (int g = 0; g < G; g++) {
        v += c[g]*exp(distinct[g]*x);
    }
    return max(v, 0.0);
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
    if (k <= 1) {
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

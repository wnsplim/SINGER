//
//  random_utils.cpp
//  SINGER
//
//  Created by Yun Deng on 5/10/23.
//

#include "random_utils.hpp"
#include <algorithm>
#include <cmath>

std::mt19937 random_engine;
std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);

double uniform_random() {
    return uniform_distribution(random_engine);
}

// geometric waiting times in place of std::binomial_distribution, whose rejection branch (n min(p, 1 - p) >= 8) is biased in libstdc++ (land/finite/da/logs/binom6b.log)
long long binomial_random(long long n, double p) {
    double q = std::min(p, 1 - p);
    long long k = 0;
    if (q > 0) {
        double scale = 1/std::log1p(-q);
        for (double at = std::floor(std::log(1 - uniform_random())*scale); at < n; at += std::floor(std::log(1 - uniform_random())*scale) + 1) {
            k++;
        }
    }
    return p > 0.5 ? n - k : k;
}

std::string get_time() {
    using namespace std::chrono;
    auto now = system_clock::now();
    auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;
    auto timer = system_clock::to_time_t(now);
    std::tm bt = *std::localtime(&timer);
    std::ostringstream oss;
    oss << "[" << std::put_time(&bt, "%H:%M:%S");
    oss << '.' << std::setfill('0') << std::setw(3) << ms.count() << "]";
    return oss.str();
}

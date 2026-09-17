//
// Created by Valentin Mayer on 16.09.26.
//

#ifndef DQMC_BLOCKING_ANALYSIS_H
#define DQMC_BLOCKING_ANALYSIS_H
#include <numeric>
#include <vector>

struct block_analysis_result {
    bool success;
    double autocorrelation_time, mean, mean_variance;
};

struct block_variance_result {
    int B;
    double mean_variance;
};


inline double delta_variance(double blocked_variance, int M) {
    return blocked_variance * std::sqrt(2. / (M - 1));
}


#endif //DQMC_BLOCKING_ANALYSIS_H

//
// Created by Valentin Mayer on 15.08.26.
//
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>
#include <numeric>

const double Z_1 = 1 - std::exp(-5); // normalization factor for alpha=1
const double Z_2 = 2 * (1 - std::exp(-5. / 2)); // normalization factor for alpha=0.5

// better to store struct or something else?
struct D {
    double alpha, tau;
};

struct block_analysis_result {
    bool success;
    double autocorrelation_time, mean, mean_variance;
};

struct block_variance_result {
    int B;
    double mean_variance;
};

double target_pdf(D diag) {
    return std::exp(-diag.tau * diag.alpha);
}

void change_tau(D &diag, const double delta, std::mt19937 &rng) {
    std::uniform_real_distribution<double> random_uniform(0, 1);

    double tau_proposed = diag.tau + delta * (random_uniform(rng) - 0.5);

    diag.tau = (tau_proposed <= 5 && tau_proposed >= 0) ? tau_proposed : diag.tau;
}

void change_alpha(D &diag, std::mt19937 &rng) {
    std::bernoulli_distribution random_bernoulli(0.5);

    diag.alpha = random_bernoulli(rng) ? 0.5 : 1;
}

std::vector<D> sample_distribution(std::mt19937 &rng, const double delta, int N = 1'000'000) {
    std::uniform_real_distribution<> random_uniform(0, 1);
    std::bernoulli_distribution random_bernoulli(0.5);

    D d_initial = {.alpha = 1, .tau = random_uniform(rng) * 5};

    std::vector<D> D_values = {d_initial};
    D_values.reserve(N);
    for (int i = 1; i < N; i++) {
        D d_curr = D_values.back();
        D d_new = d_curr;

        // 50/50 which change will be done
        random_bernoulli(rng) ? change_tau(d_new, delta, rng) : change_alpha(d_new, rng);
        double p_acc = std::min(1., target_pdf(d_new) / target_pdf(d_curr));

        if (random_uniform(rng) < p_acc) {
            D_values.push_back(d_new);
        } else {
            D_values.push_back(d_curr);
        }
    }
    return D_values;
}

double delta_variance(double blocked_variance, int M) {
    return blocked_variance * std::sqrt(2. / (M - 1));
}

template<typename Func>
block_analysis_result block_analysis(Func observable, std::vector<D> &diags) {
    int N = diags.size();

    std::vector<double> observables;
    observables.reserve(N);

    //surely there is a better way
    for (const auto &diag: diags) {
        observables.push_back(observable(diag));
    }

    std::vector<double> cumulative_sum{0.};
    cumulative_sum.resize(diags.size() + 1);

    std::partial_sum(observables.begin(), observables.end(), cumulative_sum.begin() + 1);

    double mean = cumulative_sum.back() / N;

    std::vector<block_variance_result> block_variance_results;
    int B = 1;

    bool plateau_found = false;
    std::vector<block_variance_result> plateau_points;

    while (N / B > 50) {
        int M = N / B;

        std::vector<double> block_means;
        block_means.reserve(M);

        for (int b = 0; b < M; b++) {
            double block_mean = (cumulative_sum[(b + 1) * B] - cumulative_sum[b * B]) / B;
            block_means.push_back(block_mean);
        }

        double block_mean_avg = std::accumulate(block_means.begin(), block_means.end(), 0.) / M;
        double block_mean_variance = std::accumulate(block_means.begin(), block_means.end(), 0.,
                                                     [block_mean_avg](double sum, double value) {
                                                         double diff = block_mean_avg - value;
                                                         return sum + diff * diff;
                                                     }) / M;
        double mean_variance = block_mean_variance / M;
        block_variance_results.push_back({.B = B, .mean_variance = mean_variance});

        if (block_variance_results.size() > 1 && std::abs(
                block_variance_results.back().mean_variance - (block_variance_results.rbegin() + 1)->mean_variance) <
            delta_variance(block_variance_results.back().mean_variance, M))
            plateau_found = true;

        if (plateau_found)
            plateau_points.push_back({.B = B, .mean_variance = mean_variance});

        B *= 2;
    }

    double plateau_variance = std::accumulate(plateau_points.begin(), plateau_points.end(), 0.,
                                              [N](double sum, const block_variance_result &block_variance) {
                                                  return sum + 1 / std::pow(
                                                             delta_variance(
                                                                 block_variance.mean_variance, N / block_variance.B),
                                                             2) * block_variance.mean_variance;
                                              }) / std::accumulate(plateau_points.begin(), plateau_points.end(), 0.,
                                                                   [N](double sum,
                                                                       const block_variance_result &block_variance) {
                                                                       return sum + 1 / std::pow(
                                                                                  delta_variance(
                                                                                      block_variance.mean_variance,
                                                                                      N / block_variance.B), 2);
                                                                   });

    double autocorrelation_time = 1. / 2 * plateau_variance / block_variance_results.front().mean_variance;

    if (plateau_points.empty())
        return {.success = false, .autocorrelation_time = -1, .mean = mean, .mean_variance = -1};
    return {
        .success = true, .autocorrelation_time = autocorrelation_time, .mean = mean, .mean_variance = plateau_variance
    };
}

const double I1_exact = 1. - 6. * std::exp(-5);
const double I2_exact = 2. - 37. * std::exp(-5);
const double I3_exact = 2.850810019265417;
const double I4_exact = 7.298990145866727;

int main() {
    std::random_device rd;
    std::mt19937 rng(rd());

    std::vector<D> diags = sample_distribution(rng, 10, 1'000'000);

    std::ofstream histogram_output(std::string(PROJECT_ROOT) + "/data/task5/histogram.csv");
    histogram_output << "alpha,tau\n";

    std::vector<D> diags_alpha_1;
    std::vector<D> diags_alpha_0_5;

    std::partition_copy(diags.begin(), diags.end(), std::back_inserter(diags_alpha_1),
                        std::back_inserter(diags_alpha_0_5), [](const D &diag) {
                            return diag.alpha == 1;
                        });

    auto report_integral = [](const std::string &name, const block_analysis_result &result,
                              const double normalization, const double exact) {
        const double value = normalization * result.mean;
        const double error = normalization * std::sqrt(result.mean_variance);
        std::cout << name << ": " << value << " +/- " << error
                << ", exact: " << exact << ", difference: " << std::abs(value - exact) << "\n";
    };

    block_analysis_result result1 = block_analysis([](const D &diag) { return diag.tau; }, diags_alpha_1);
    report_integral("I1 (alpha=1)", result1, Z_1, I1_exact);

    block_analysis_result result2 = block_analysis([](const D &diag) { return diag.tau * diag.tau; }, diags_alpha_1);
    report_integral("I2 (alpha=1)", result2, Z_1, I2_exact);

    block_analysis_result result3 = block_analysis([](const D &diag) { return diag.tau; }, diags_alpha_0_5);
    report_integral("I3 (alpha=0.5)", result3, Z_2, I3_exact);

    block_analysis_result result4 = block_analysis([](const D &diag) { return diag.tau * diag.tau; }, diags_alpha_0_5);
    report_integral("I4 (alpha=0.5)", result4, Z_2, I4_exact);

    for (const auto [alpha, tau]: diags) {
        histogram_output << alpha << "," << tau << "\n";
    }
    return 0;
}

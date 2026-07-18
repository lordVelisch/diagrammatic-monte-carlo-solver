#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <numeric>

struct integration_result {double value, error;};

const double NORMALIZATION_FACTOR = 1-std::exp(-5.0);

double sample_tau(std::mt19937& rng) {
    std::uniform_real_distribution<double> dist(0,1);
    double r = dist(rng);
    ;
    double tau = -std::log(1-r*NORMALIZATION_FACTOR);
    return tau;
}

template <typename Func>
integration_result integrate(std::mt19937& rng,Func g, int N=1'000'000) {
    std::vector<double> values;
    values.reserve(N);

    for (int i=0; i<N; i++) {
        double tau = sample_tau(rng);
        values.push_back(g(tau));
    }

    double mean = std::accumulate(values.begin(), values.end(), 0.)/N;

    double g_variance = std::accumulate(values.begin(), values.end(), 0., [mean](double sum, double value) {
        double d = value-mean;
        return sum+d*d;
    }) / N;

    double sd_of_mean = std::sqrt(1./N*g_variance);
    return {mean, sd_of_mean};
}

int main() {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::ofstream cdf_inverse_out(std::string(PROJECT_ROOT) + "/data/task3/cdf_inverse_data.csv");
    cdf_inverse_out << "value" << "\n";

    for (int i=0; i<100'000; i++) {
        double tau = sample_tau(rng);
        cdf_inverse_out << tau << "\n";
    }

    integration_result result1 = integrate(rng, [](double tau) {return tau*NORMALIZATION_FACTOR;});
    integration_result result2 = integrate(rng, [](double tau) {return tau*tau*NORMALIZATION_FACTOR;});

    std::cout <<  "The result of the first intergral is: " << result1.value << "; error: " << result1.error << "\n";
    std::cout <<  "The result of the second intergral is: " << result2.value << "; error: " << result2.error <<"\n";
};

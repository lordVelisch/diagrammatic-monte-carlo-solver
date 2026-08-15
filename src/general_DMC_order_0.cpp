//
// Created by Valentin Mayer on 15.08.26.
//
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <numeric>

const double Z = 1 - std::exp(-5); // normalization factor

// better to store struct or something else?
struct D {
    double alpha, tau;
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

const double I1_exact = 1. - 6. * std::exp(-5);
const double I2_exact = 2. - 37. * std::exp(-5);

int main() {
    std::random_device rd;
    std::mt19937 rng(rd());

    std::vector<D> diags = sample_distribution(rng, 10, 1'000'000);

    std::ofstream histogram_output(std::string(PROJECT_ROOT) + "/data/task5/histogram.csv");
    histogram_output << "alpha,tau\n";

    // I could do another blocking analysis here (or insert the code to check the actual errors)

    for (const auto [alpha, tau] : diags) {
        histogram_output << alpha << "," << tau << "\n";
    }
    return 0;
}

//
// Created by Valentin Mayer on 15.08.26.
//
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <numeric>

// better to store struct or something else?
struct D {
    double alpha, tau;
};

double target_pdf(D diag) {
    return std::exp(-diag.tau * diag.alpha);
}

const double Z = 1 - std::exp(-5); // normalization factor


// reference okay here? Idea was maybe I can easier pass the diag_new into multiple change_x functions
// but maybe it is just better to create inside the function and return if there is at maximum only 1 mutation done
void change_tau(D &diag, const double delta, std::mt19937 &rng) {
    std::uniform_real_distribution<double> random_uniform(0, 1);

    double tau_proposed = diag.tau + delta * (random_uniform(rng) - 0.5);
    // this way it will still calculated and get 1 as result.
    //find better solution long term

    diag.tau = (tau_proposed <= 5 && tau_proposed >= 0) ? tau_proposed : diag.tau;
}

void change_alpha(D &diag, std::mt19937 &rng) {
    std::bernoulli_distribution random_bernoulli(0.5);

    diag.alpha = random_bernoulli(rng) ? 0.5 : 1;
}

// A is the observable to calculate
template<typename Func>
double sample_distribution(std::mt19937 &rng, const double delta, Func observable,
                                        int N = 1'000'000, std::ofstream* histogram_output = nullptr) {

    // file writing for plotting of diagrams, this is a bit bothersome to have inside
    // and if statements with write_to_file everytime also seems a bit clunky. better way available or fine?
    if (histogram_output)
        *histogram_output << "alpha,tau\n";

    std::uniform_real_distribution<double> random_uniform(0, 1);
    std::bernoulli_distribution random_bernoulli(0.5);

    D D_initial = {.alpha = 1, .tau = random_uniform(rng) * 5};

    double A_sum = 0;

    D d_curr = D_initial;

    for (int i = 1; i < N; i++) {

        D diag_new = d_curr;

        // 50/50 which change will be done
        random_bernoulli(rng) ? change_tau(diag_new, delta, rng) : change_alpha(diag_new, rng);

        if (histogram_output)
            *histogram_output << d_curr.alpha << "," << d_curr.tau << "\n";

        double p_acc = std::min(1., target_pdf(diag_new) / target_pdf(d_curr));

        if (random_uniform(rng) < p_acc) {
            A_sum+=observable(diag_new);
            d_curr = diag_new;
        } else {
            A_sum+=observable(d_curr); // this I could theoretically store in temp variable
        }
    }
    return A_sum/N;
}


const double I1_exact = 1. - 6. * std::exp(-5);
const double I2_exact = 2. - 37. * std::exp(-5);

int main() {
    std::random_device rd;
    std::mt19937 rng(rd());

    // just for histogram
    std::ofstream histogram_output(std::string(PROJECT_ROOT) + "/data/task5/histogram.csv");

    sample_distribution(rng, 10, [](D diag) {return diag.tau;}, 1'000'000, &histogram_output);
    sample_distribution(rng, 10, [](D diag) { return diag.alpha; }, 1'000'000);



    return 0;
}

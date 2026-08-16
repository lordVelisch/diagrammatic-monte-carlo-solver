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

const double P_add = 0.1;
const double P_rem = 0.1;

const double DELTA = 10;
// when to do it like this vs using define vs using it as function input? I am doing this now so that I dont have to pass it always

// I guess it would be good to do a general diagram and specific implementations but for now I will just do 1, but I am not sure inheritance is not a bit overkill here.
struct D {
    int order;
    double alpha, beta, tau, tau_1, tau_2, V;
    // should V and alpha be part of the Diagram? it only and external parameter right? seems like it could be external but the task describes the diagram like this
};

bool is_valid(D d) {
    // dont hardcode the 5
    if (d.order == 0)
        return d.tau >= 0 && d.tau <= 5;
    if (d.order == 2)
        return d.tau_1 >= 0 && d.tau_2 >= d.tau_1 && d.tau >= d.tau_2 && d.tau <= 5;
    throw std::runtime_error("Only order 0 and 2 diagrams are implemented so far!");
}

double target_pdf(D d) {
    if (!is_valid(d)) return 0;
    if (d.order == 0)
        return std::exp(-d.tau * d.alpha);
    if (d.order == 2) {
        return std::exp(-d.alpha * d.tau_1) * d.V * std::exp(
                   -d.beta * (d.tau_2 - d.tau_1)) * d.V * std::exp(-d.alpha * (d.tau - d.tau_2));
    }
    throw std::runtime_error("Order should only be 0 or 2!");
}

void change_tau(D &d_new, D &d_curr, std::mt19937 &rng, double &p_acc) {
    std::uniform_real_distribution<double> random_uniform(0, 1);

    d_new.tau = d_curr.tau + DELTA * (random_uniform(rng) - 0.5);

    p_acc = std::min(1., target_pdf(d_new) / target_pdf(d_curr));
}

void add_beta(D &d_new, D &d_curr, std::mt19937 &rng, double &p_acc) {
    std::uniform_real_distribution<> random_uniform(0, 1);
    std::bernoulli_distribution random_bernoulli(0.5);

    if (d_curr.order == 2) {
        p_acc = 0;
        return;
    }

    if (d_curr.order == 0) {
        // this already caused issues with p_acc
        d_new.order = 2;
        d_new.tau_1 = random_uniform(rng) * 5;
        d_new.tau_2 = (5 - d_new.tau_1) * random_uniform(rng) + d_new.tau_1;
        d_new.beta = random_bernoulli(rng) ? 1. / 4 : 3. / 4; // todo add beta in parameter struct
        d_new.V = 0.5; // a bit incidental, but want to keep as a variable somewhere, see above todo
    }
    p_acc = std::min(1., target_pdf(d_new) / (target_pdf(d_curr) * 1. / (5 - d_new.tau_1) * 1. / 5 * 1. / 2));
}

void remove_beta(D &d_new, D &d_curr, double &p_acc) {
    d_new.order = 0;
    p_acc = std::min(1., target_pdf(d_new) * 1. / (5 - d_curr.tau_1) * 1. / 5 * 1. / 2 / target_pdf(d_curr));
    // I could set everything to 0 but later I will separate the types anyway I think
    // also I will apply the function based on order anyway (but could be not explicit enough)
}

std::vector<D> sample_diagrams(std::mt19937 &rng, int N = 1'000'000) {
    std::uniform_real_distribution<> random_uniform(0, 1);

    D d_initial = {.order = 0, .alpha = 1, .tau = random_uniform(rng) * 5};
    // I should do this actually random or use some warmup

    std::vector D_values = {d_initial};
    D_values.reserve(N);
    for (int i = 1; i < N; i++) {
        D d_curr = D_values.back();
        D d_new = d_curr;

        double rand = random_uniform(rng);
        double p_acc;
        if (rand < P_add)
            // this can be built in a general way. Todo: build a function that takes a list of actions and percentages and calls the one it needs. Only question about this: I pass p_acc because I dont want to have to deal with separate return values and the acceptance criteria always looks differnet. is this the best way? or maybe a function call p_acc = execute_action[list of actions and probabilities]
            add_beta(d_new, d_curr, rng, p_acc);
        else if (rand > (1 - P_rem))
            remove_beta(d_new, d_curr, p_acc);
        else
            change_tau(d_new, d_curr, rng, p_acc);

        // this waits for p_acc to be set, right?
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

    int N = 10'000'000;
    std::vector<D> d_values = sample_diagrams(rng, N);
    std::cout << "count" << d_values.size() << "\n";
    std::ofstream histogram_output(std::string(PROJECT_ROOT) + "/data/task6/histogram.csv");
    histogram_output << "tau\n";

    double mean = std::accumulate(d_values.begin(), d_values.end(), 0., [](double sum, D d) { return sum + d.tau; }) /
                  N;

    std::cout << "Mean result is: " << mean << "\n";

    //todo: here I need to multiply by the normalization constant
    //todo: to calculate the standard deviation I have to do a block analysis again. I was thinking about generalizing it in a separate file and reusing it since it keep coming up, but need to finalize the datatype first I guess.

    for (const auto diag: d_values) {
        histogram_output << diag.tau << "\n";
    }
    return 0;
}

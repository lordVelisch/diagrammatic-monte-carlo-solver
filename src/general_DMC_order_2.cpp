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

#include "blocking_analysis.h"

const double Z0_analytical = 1 - std::exp(-5); // normalization factor for 0th order diagram
const double Z2_analytical = 0.782413501449566; // normalization factor for 2nd order diagram

const double P_add = 0.1;
const double P_rem = 0.1;

const double DELTA = 10;

// I guess it would be good to do a general diagram and specific implementations but for now I will just do 1, but I am not sure inheritance is not a bit overkill here.
struct D
{
    int order;
    double alpha, beta, tau, tau_1, tau_2, V;
    // should V and alpha be part of the Diagram? it only and external parameter right? seems like it could be external but the task describes the diagram like this
};

struct sampling_point
{
    int order;
    double tau;
};

bool is_valid(D d)
{
    // todo dont hardcode the 5
    if (d.order == 0)
        return d.tau >= 0 && d.tau <= 5;
    if (d.order == 2)
        return d.tau_1 >= 0 && d.tau_2 >= d.tau_1 && d.tau >= d.tau_2 && d.tau <= 5;
    throw std::runtime_error("Only order 0 and 2 diagrams are implemented so far!");
}

double target_pdf(D d)
{
    if (!is_valid(d)) return 0;
    if (d.order == 0)
        return std::exp(-d.tau * d.alpha);
    if (d.order == 2)
    {
        return std::exp(-d.alpha * d.tau_1) * d.V * std::exp(
            -d.beta * (d.tau_2 - d.tau_1)) * d.V * std::exp(-d.alpha * (d.tau - d.tau_2));
    }
    throw std::runtime_error("Order should only be 0 or 2!");
}

void change_tau(D& d_new, D& d_curr, std::mt19937& rng, double& p_acc)
{
    std::uniform_real_distribution<double> random_uniform(0, 1);

    d_new.tau = d_curr.tau + DELTA * (random_uniform(rng) - 0.5);

    p_acc = std::min(1., target_pdf(d_new) / target_pdf(d_curr));
}

void add_beta(D& d_new, D& d_curr, std::mt19937& rng, double& p_acc)
{
    std::uniform_real_distribution<> random_uniform(0, 1);
    std::bernoulli_distribution random_bernoulli(0.5);

    if (d_curr.order == 2)
    {
        p_acc = 0;
        return;
    }

    if (d_curr.order == 0)
    {
        d_new.order = 2;
        d_new.tau_1 = random_uniform(rng) * 5;
        d_new.tau_2 = (5 - d_new.tau_1) * random_uniform(rng) + d_new.tau_1;
        d_new.beta = random_bernoulli(rng) ? 1. / 4 : 3. / 4; // todo add beta in parameter struct
        d_new.V = 0.5; // a bit incidental, but want to keep as a variable somewhere, see above todo
    }
    p_acc = std::min(1., target_pdf(d_new) / (target_pdf(d_curr) * 1. / (5 - d_new.tau_1) * 1. / 5 * 1. / 2));
}

void remove_beta(D& d_new, D& d_curr, double& p_acc)
{
    d_new.order = 0;
    p_acc = std::min(1., target_pdf(d_new) * 1. / (5 - d_curr.tau_1) * 1. / 5 * 1. / 2 / target_pdf(d_curr));
    // I could set everything to 0 but later I will separate the types anyway I think
    // also I will apply the function based on order anyway (but could be not explicit enough)
}

std::vector<sampling_point> sample_diagrams(std::mt19937& rng, int N = 1'000'000)
{
    std::uniform_real_distribution<> random_uniform(0, 1);

    D d_init = {.order = 0, .alpha = 1, .tau = random_uniform(rng) * 5};
    // I should do this actually random or use some warmup

    std::vector<sampling_point> obs_values = {{d_init.order, d_init.tau}};
    D d_curr = d_init;
    obs_values.reserve(N);
    for (int i = 1; i < N; i++)
    {
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

        if (random_uniform(rng) < p_acc)
        {
            obs_values.push_back({d_new.order, d_new.tau});
            d_curr = d_new;
        }
        else
        {
            obs_values.push_back({d_curr.order, d_curr.tau});
        }
    }
    return obs_values;
}

//todo: this is currently not ideal: there might be some divisions by 0 for small B and I basically mix the results
template <typename Func>
block_analysis_result block_analysis(Func observable, std::vector<sampling_point>& sample_points)
{
    int N = sample_points.size();

    std::vector<double> observables;
    observables.reserve(N);

    std::vector<int> order_0_points = {0};
    order_0_points.reserve(N + 1);

    //surely there is a better way
    for (const auto& point : sample_points)
    {
        observables.push_back(observable(point));
        order_0_points.push_back(order_0_points[order_0_points.size() - 1] + (point.order == 0 ? 1 : 0));
    }

    std::vector cumulative_sum{0.};
    cumulative_sum.resize(sample_points.size() + 1);

    std::partial_sum(observables.begin(), observables.end(), cumulative_sum.begin() + 1);

    double mean = cumulative_sum.back() / N;

    std::vector<block_variance_result> block_variance_results;
    int B = 1;

    bool plateau_found = false;
    std::vector<block_variance_result> plateau_points;

    while (N / B > 50)
    {
        int M = N / B;

        std::vector<double> block_means;
        block_means.reserve(M);

        for (int b = 0; b < M; b++)
        {
            double block_mean = (cumulative_sum[(b + 1) * B] - cumulative_sum[b * B]) / (order_0_points[(b + 1) * B] -
                order_0_points[b * B]);
            block_means.push_back(block_mean);
        }

        double block_mean_avg = std::accumulate(block_means.begin(), block_means.end(), 0.) / M;
        double block_mean_variance = std::accumulate(block_means.begin(), block_means.end(), 0.,
                                                     [block_mean_avg](double sum, double value)
                                                     {
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
                                              [N](double sum, const block_variance_result& block_variance)
                                              {
                                                  return sum + 1 / std::pow(
                                                      delta_variance(
                                                          block_variance.mean_variance, N / block_variance.B),
                                                      2) * block_variance.mean_variance;
                                              }) / std::accumulate(plateau_points.begin(), plateau_points.end(), 0.,
                                                                   [N](double sum,
                                                                       const block_variance_result& block_variance)
                                                                   {
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

int main()
{
    std::random_device rd;
    std::mt19937 rng(rd());

    int N = 10'000'000;
    std::vector<sampling_point> sampling_points = sample_diagrams(rng, N);

    const size_t order_0_count = std::count_if(sampling_points.begin(), sampling_points.end(),
                                               [](sampling_point p) { return p.order == 0; });
    const size_t order_2_count = sampling_points.size() - order_0_count;
    std::cout << "Ratio of expected order 0 diagrams vs actual order 0 diagrams: \n";
    std::cout << "Expected: " << Z0_analytical / (Z0_analytical + Z2_analytical) << "\n";
    std::cout << "Actual: " << static_cast<double>(order_0_count) / (order_0_count + order_2_count) << "\n";

    const double Z = Z0_analytical * N / order_0_count;

    auto report_integral = [](const std::string& name, const block_analysis_result& result,
                              const double normalization, const double exact)
    {
        const double value = normalization * result.mean;
        // todo: this is currently a bit weird because I compute the mean just for the observable, but the mean_variance already takes Z into account only leaving a factor of Z_0 -> there are better solutions to this
        const double error = Z0_analytical * std::sqrt(result.mean_variance);
        std::cout << name << ": " << value << " +/- " << error
            << ", exact: " << exact << ", difference: " << std::abs(value - exact) << "\n";
    };

    block_analysis_result result1 = block_analysis([](const sampling_point& point) { return point.tau; },
                                                   sampling_points);
    const double I1_exact = 3.2043134579255224;
    const double I2_exact = 9.281685042355228;

    report_integral("I1", result1, Z, I1_exact);

    block_analysis_result result2 = block_analysis([](const sampling_point& point) { return point.tau * point.tau; },
                                                   sampling_points);
    report_integral("I2", result2, Z, I2_exact);

    // Histogram for plotting:

    std::ofstream histogram_output(std::string(PROJECT_ROOT) + "/data/task6/histogram.csv");
    histogram_output << "tau\n";

    for (const auto diag : sampling_points)
    {
        histogram_output << diag.tau << "\n";
    }
    return 0;
}

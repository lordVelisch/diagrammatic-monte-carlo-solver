#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <numeric>

struct sample_distribution_result {
    double acceptance_ratio;
    std::vector<double> values;
};

struct block_analysis_result {
    bool success;
    double autocorrelation_time, mean_variance;
};

struct block_variance_result {
    int B;
    double mean_variance;
};

double target_pdf(const double x) {
    return std::exp(-x);
}

const double Z = 1-std::exp(-5); // normalization factor


template <typename Func>
sample_distribution_result sample_distribution(std::mt19937 &rng, const double delta, Func observable,
                                               int N = 1'000'000, bool write_to_file = false) {
    std::uniform_real_distribution<double> random_uniform(0, 1);


    double initial_x = random_uniform(rng);

    std::vector<double> A_values{observable(initial_x)};
    A_values.reserve(N);

    double x = initial_x;

    size_t rejection_count = 0;

    for (int i = 1; i < N; i++) {
        double A = A_values.back();

        double x_proposed = x + delta * (random_uniform(rng) - 0.5);
        if (x_proposed < 0 || x_proposed > 5) {
            A_values.push_back(A);
            rejection_count++;
            continue;
        }

        double A_proposed = observable(x_proposed);

        double p_acc = std::min(1., target_pdf(x_proposed) / target_pdf(x));

        if (random_uniform(rng) < p_acc) {
            A_values.push_back(A_proposed);
            x = x_proposed;
        } else {
            A_values.push_back(A);
            rejection_count++;
        }
    }

    if (write_to_file) {
        std::ofstream histogram_output(std::string(PROJECT_ROOT) + "/data/task4/histogram.csv");
        histogram_output << "tau\n";

        std::copy(A_values.begin(), A_values.end(), std::ostream_iterator<double>(histogram_output, "\n"));
    }
    return {.acceptance_ratio = 1. - static_cast<double>(rejection_count) / N, .values = A_values};
}

double delta_variance(double blocked_variance, int M) {
    return blocked_variance * std::sqrt(2. / (M - 1));
}

template <typename Func>
block_analysis_result block_analysis(std::mt19937& rng, const int N, double delta, Func observable, std::ofstream* block_convergence_output = nullptr, std::ofstream* delta_output_file = nullptr) {
    sample_distribution_result result = sample_distribution(rng, delta, observable, N, false);

    std::vector<double> cumulative_sum {0.};
    cumulative_sum.resize(result.values.size()+1);

    std::partial_sum(result.values.begin(), result.values.end(), cumulative_sum.begin()+1);

    std::vector<block_variance_result> block_variance_results;
    int B = 1;

    bool plateau_found = false;
    std::vector<block_variance_result> plateau_points;

    while (N/B > 50) {
        int M = N/B;

        std::vector<double> block_means;
        block_means.reserve(M);

        for (int b=0; b<M; b++) {
            double block_mean = (cumulative_sum[(b+1)*B]-cumulative_sum[b*B])/B;
            block_means.push_back(block_mean);
        }

        double block_mean_avg = std::accumulate(block_means.begin(), block_means.end(), 0.)/M;
        double block_mean_variance = std::accumulate(block_means.begin(), block_means.end(), 0., [block_mean_avg] (double sum, double value) {
            double diff = block_mean_avg-value;
            return sum + diff*diff;
        }) / M;
        double mean_variance = block_mean_variance/M;
        block_variance_results.push_back({.B = B, .mean_variance = mean_variance});

        if (block_variance_results.size() > 1 && std::abs(block_variance_results.back().mean_variance - (block_variance_results.rbegin()+1)->mean_variance) < delta_variance(block_variance_results.back().mean_variance, M))
            plateau_found = true;

        if (plateau_found)
            plateau_points.push_back({.B = B, .mean_variance = mean_variance});

        B*=2;
    }

    double plateau_variance = std::accumulate(plateau_points.begin(), plateau_points.end(), 0., [N](double sum, const block_variance_result& block_variance) {
        return sum+1/std::pow(delta_variance(block_variance.mean_variance,N/block_variance.B), 2)*block_variance.mean_variance;
    }) / std::accumulate(plateau_points.begin(), plateau_points.end(), 0., [N](double sum, const block_variance_result& block_variance) {
        return sum + 1/std::pow(delta_variance(block_variance.mean_variance, N/block_variance.B), 2);
    });

    double autocorrelation_time = 1./2 * plateau_variance/block_variance_results.front().mean_variance;

    if (block_convergence_output)
        for (auto variance_at_b: block_variance_results) *block_convergence_output << variance_at_b.B << "," << variance_at_b.mean_variance << "," << delta << "\n";
    if (delta_output_file)
        *delta_output_file << delta << "," << autocorrelation_time << "," << plateau_variance << "\n";

    if (plateau_points.empty())
        return {.success = false, .autocorrelation_time = -1, .mean_variance = -1};
    return {.success = true, .autocorrelation_time = autocorrelation_time, .mean_variance = plateau_variance};
}

template <typename Func>
double tau_int_optimization(std::mt19937& rng, Func observable, int N=1'000'000) {
    std::vector<double> trial_deltas{0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100};

    std::ofstream block_convergence_output(std::string(PROJECT_ROOT) + "/data/task4/block_convergence.csv");
    block_convergence_output << "B,error,delta\n";

    std::ofstream delta_values_output(std::string(PROJECT_ROOT) + "/data/task4/delta_values.csv");
    delta_values_output << "delta,tao_int,var_inf\n";

    std::vector<std::pair<double, double>> results;
    for (double delta: trial_deltas) {
        block_analysis_result result = block_analysis(rng, N, delta, observable, &block_convergence_output, &delta_values_output);
        if (result.success)
            results.emplace_back(delta, result.autocorrelation_time);
    }

    double delta_opt = std::min_element(results.begin(), results.end(), [](const auto& result1, const auto& result2) {
        return result1.second < result2.second;
    })->first;

    return delta_opt;
}

const double I1_exact = 1.-6.*std::exp(-5);
const double I2_exact = 2.-37.*std::exp(-5);

int main() {
    std::random_device rd;
    std::mt19937 rng(rd());

    // just for historgram
    sample_distribution_result result = sample_distribution(rng, 10, [](double x) {return x;}, 1'000'000, true);

     double delta_opt = tau_int_optimization(rng, [](const double A){return A;});

     std::cout << "The optimal delta, that is minimizing the autocorrelation times is: " << delta_opt << "\n";

     // solving integral I1:
     int N = 1'000'000;
     auto x_func = [](double x) {return x;};
     std::vector<double> samples_1 = sample_distribution(rng, delta_opt, x_func, N).values;
     double x_mean = std::accumulate(samples_1.begin(), samples_1.end(), .0)/N;
     double x_sigma = std::sqrt(block_analysis(rng, N, delta_opt, x_func).mean_variance);

     double I1 = Z*x_mean;

     std::cout << "Integral I1 result: " << I1 << " pm " << x_sigma <<"\n";
     std::cout << "Difference to exact: " << I1 - I1_exact << "\n";
     std::cout << "Accpetance ratio: " << result.acceptance_ratio << "\n";
     std::cout << "\n";

     // solving integral I2:

     auto x2_func = [](double x) {return x*x;};

     std::vector<double> samples_2 = sample_distribution(rng, delta_opt, x2_func, N).values;
     double x2_mean = std::accumulate(samples_2.begin(), samples_2.end(), .0)/N;
     double x2_sigma = std::sqrt(block_analysis(rng, N, delta_opt, x2_func).mean_variance); // the delta_opt is used from <x>_Q because it is not this important

     double I2 = Z*x2_mean;
     std::cout << "Integral I2 result: " << I2 << " pm " << x2_sigma << "\n";
     std::cout << "Difference to exact: " << I2 - I2_exact << "\n";

    return 0;
}

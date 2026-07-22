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
    double autocorrelation_time;
    int plateau_block_size;
};

double f(const double tau) {
    return std::exp(-tau);
}

template <typename Func>
sample_distribution_result sample_distribution(std::mt19937 &rng, const double DELTA, Func estimation_value_of,
                                               int N = 1'000'000, bool write_to_file = false) {
    std::uniform_real_distribution<double> random_uniform(0, 1);


    double initial_x = random_uniform(rng);

    std::vector<double> A_values{estimation_value_of(initial_x)};
    A_values.reserve(N);

    double x = initial_x;

    size_t rejection_count = 0;

    for (int i = 1; i < N; i++) {
        double A = A_values.back();

        double x_proposed = x + DELTA * (random_uniform(rng) - 0.5);
        if (x_proposed < 0 || x_proposed > 5) {
            A_values.push_back(A);
            rejection_count++;
            continue;
        }

        double A_proposed = estimation_value_of(x_proposed);

        double p_acc = std::min(1., f(x_proposed) / f(x));

        if (random_uniform(rng) < p_acc) {
            A_values.push_back(A_proposed);
            x = x_proposed;
        } else {
            A_values.push_back(A);
            rejection_count++;
        }
    }

    if (write_to_file) {
        std::ofstream out(std::string(PROJECT_ROOT) + "/data/task4/histogram.csv");
        out << "tau\n";

        std::copy(A_values.begin(), A_values.end(), std::ostream_iterator<double>(out, "\n"));
    }
    return {1. - static_cast<double>(rejection_count) / N, A_values};
}

std::vector<double>::value_type delta_variance(double blocked_variance, int M) {
    return blocked_variance * std::sqrt(2. / (M - 1));
}

block_analysis_result block_analysis(std::mt19937& rng, std::ofstream& block_convergence_output, std::ofstream& out, const int N, double delta) {
    sample_distribution_result result = sample_distribution(rng, delta, [](double A){return A;}, N, true);

    double tol = 1e-6;

    std::vector<double> cumulative_sum {0.};
    cumulative_sum.resize(result.values.size()+1);

    std::partial_sum(result.values.begin(), result.values.end(), cumulative_sum.begin()+1);

    std::vector<std::pair<int, double>> blocked_variances;
    int B = 1;

    int M = N/B;

    bool plateau_found = false;
    std::vector<std::pair<int, double>> plateau_points;

    while (M > 50) {
        M = N/B;

        std::vector<double> block_means;
        block_means.reserve(M);

        for (int b=0; b<M; b++) {
            double block_mean = (cumulative_sum[(b+1)*B]-cumulative_sum[b*B])/B;
            block_means.push_back(block_mean);
        }

        double sample_average = std::accumulate(block_means.begin(), block_means.end(), 0.)/M;
        double sample_variance = std::accumulate(block_means.begin(), block_means.end(), 0., [sample_average] (double sum, double value) {
            double diff = sample_average-value;
            return sum + diff*diff;
        }) / M;
        double blocked_variance = sample_variance/M;
        blocked_variances.emplace_back(B, blocked_variance);

        if ( blocked_variances.size() > 1 && std::abs(blocked_variances.back().second - (blocked_variances.rbegin()+1)->second) < delta_variance(blocked_variances.back().second, M))
            plateau_found = true;

        if (plateau_found)
            plateau_points.emplace_back(B, blocked_variance);

        B*=2;
    }

    for (auto variance_at_b: blocked_variances) block_convergence_output << variance_at_b.first << "," << variance_at_b.second << "," << delta << "\n";

    double plateau_variance = std::accumulate(plateau_points.begin(), plateau_points.end(), 0., [N](double sum, std::pair<int, double> plateau_variance) {
        return sum+1/std::pow(delta_variance(plateau_variance.second,N/plateau_variance.first), 2)*plateau_variance.second;
    }) / std::accumulate(plateau_points.begin(), plateau_points.end(), 0., [N](double sum, std::pair<int, double> plateau_variance) {
        return sum + 1/std::pow(delta_variance(plateau_variance.second, N/plateau_variance.first), 2);
    });

    double autocorrelation_time = 1./2 * plateau_variance/blocked_variances.front().second;

    out << delta << "," << autocorrelation_time << "," << plateau_variance << "\n";

    if (plateau_points.empty())
        return {false, 0, 0};
    return {true, autocorrelation_time, plateau_points.front().first};
}

int main() {
    std::random_device rd;
    std::mt19937 rng(rd());

    /*
    sample_distribution_result result = sample_distribution(rng, 0.3, 1, 1'000'000, true);

    std::cout << "Acceptance ratio: " << result.acceptance_ratio << "\n";
    */

    std::vector<double> trial_deltas{0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100};
    //std::vector<double> trial_deltas{1};

    std::ofstream block_convergence_output(std::string(PROJECT_ROOT) + "/data/task4/block_convergence.csv");
    block_convergence_output << "B,error,delta\n";

    std::ofstream out(std::string(PROJECT_ROOT) + "/data/task4/delta_values.csv");
    out << "delta,tao_int,var_inf\n";

    const int N = 10'000'000;

    std::vector<std::pair<double, double>> results;
    for (double delta: trial_deltas) {
        block_analysis_result result = block_analysis(rng, block_convergence_output, out, N, delta);
        if (result.success)
            results.emplace_back(delta, result.autocorrelation_time);
    }

    double delta_opt = std::min_element(results.begin(), results.end(), [](const auto result1, const auto result2) {
        return result1.second < result2.second; // is there a nicer way like in java because i just pass the autocorrelation_time
    })->first;

    std::cout << "The optimal delta, that is minimizing the autocorrelation times is: " << delta_opt << "\n";
    return 0;
}

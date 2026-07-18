#include <fstream>
#include <iostream>
#include <random>
#include <filesystem>

template <typename Callback>
double estimate_pi_area(int N, std::mt19937& rng, Callback&& on_sample) {
    int N_circle_count = 0;
    std::uniform_real_distribution<double> dist(0,1);

    for (int i=0; i<N; i++) {
        double x = dist(rng);
        double y = dist(rng);

        if (x*x + y*y <= 1)
            N_circle_count++;

        on_sample(i, static_cast<double>(N_circle_count)/(i+1)*4);
    }

    return static_cast<double>(N_circle_count)/N*4;
}

double f(double x) {
    return std::sqrt(1-x*x);
}

template <typename Callback>
double estimate_pi_integration(int N, std::mt19937& rng, Callback&& on_sample) {
    std::uniform_real_distribution<double> dist(0,1);

    double sum=0;
    for (int i=0; i<N;i++) {
        sum+=f(dist(rng));
        on_sample(i, 1./(i+1)*sum*4);
    }

    return 1./N*sum*4;
}

int main() {
    std::mt19937 rng(std::random_device{}());

    std::ofstream out(std::string(PROJECT_ROOT) + "/data/task2/pi_estimation_convergence.csv");
    if (!out)
        std::cerr << "failed to open output file \n";


    out << "N,method,pi_estimate\n";

    estimate_pi_area(1000000, rng, [&out] (int i, double pi_estimation) {
        out << i << ",area," << pi_estimation << "\n";
    });
    estimate_pi_integration(1000000, rng, [&out] (int i, double pi_estimation) {
        out << i << ",integration," << pi_estimation << "\n";
    });

    std::ofstream dist_out(std::string(PROJECT_ROOT) + "/data/task2/pi_estimation_distribution.csv");
    if (!dist_out) std::cerr << "failed to open output file \n";

    dist_out << "M,method,pi_estimate\n";

    int M = 1000;
    for (int i=0; i<M; i++) {
        double pi_estimate_area = estimate_pi_area(1000000, rng, [](int, double){});
        dist_out << i << ",area," << pi_estimate_area << "\n";

        double pi_estimate_integration = estimate_pi_integration(1000000, rng, [](int, double){});
        dist_out << i << ",integration," << pi_estimate_integration << "\n";
    }

    return 0;
}

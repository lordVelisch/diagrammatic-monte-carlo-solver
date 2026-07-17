#include <fstream>
#include <iostream>
#include <random>
#include <filesystem>

double estimate_pi(int N, std::mt19937& rng, std::ofstream& out) {
    int N_circle_count = 0;

    std::uniform_real_distribution<double> dist(0,1);

    for (int i=0; i<N; i++) {
        double x = dist(rng);
        double y = dist(rng);

        if (x*x + y*y <= 1)
            N_circle_count++;

        out << i << "," << "area," << static_cast<double>(N_circle_count)/(i+1)*4 << "\n";
    }

    double pi_estimation = static_cast<double>(N_circle_count)/N*4;

    return pi_estimation;
}

int main() {

    std::cout << "cwd: " << std::filesystem::current_path() << "\n";

    std::mt19937 rng(std::random_device{}());

    std::ofstream out(std::string(PROJECT_ROOT) + "/data/task2/pi_estimation_convergence.csv");
    if (!out)
        std::cerr << "failed to open output file \n";


    out << "N,method,pi_estimate\n";

    double result = estimate_pi(1000000, rng, out);

    return 0;
}

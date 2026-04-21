#pragma once
#include <random>
#include <utility>
#include <vector>

class PoissonProcess {
public:
    static std::vector<std::pair<double, int>> simulate(double T, double lambda) {
        std::mt19937 gen(std::random_device{}());
        std::exponential_distribution<double> exp_dist(lambda);
        std::vector<std::pair<double, int>> events;
        double t = 0.0;
        int count = 0;
        while (true) {
            t += exp_dist(gen);
            if (t >= T)
                break;
            events.emplace_back(t, ++count);
        }
        return events;
    }
};

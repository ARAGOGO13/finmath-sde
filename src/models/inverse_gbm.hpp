#pragma once
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include "../sde_base.hpp"

class InverseGBM : public SDE {
public:
    double mu;
    double sigma;
    double S0;

    InverseGBM(double mu, double sigma, double S0) : mu(mu), sigma(sigma), S0(S0) {}

    double drift(double x, double t) const override { return (sigma * sigma - mu) * x; }

    double diffusion(double x, double t) const override { return -sigma * x; }

    double diffusion_derivative(double x, double t) const override { return -sigma; }

    double theoretical_mean(double T) const { return (1.0 / S0) * std::exp((sigma * sigma - mu) * T); }

    double theoretical_var(double T) const {
        double mu_star = sigma * sigma - mu;
        return (1.0 / (S0 * S0)) * std::exp(2.0 * mu_star * T) * (std::exp(sigma * sigma * T) - 1.0);
    }

    std::vector<double> simulate_via_inverse(double T, int steps, int N, std::mt19937 &rng) const {
        double dt = T / steps;
        double sqrt_dt = std::sqrt(dt);
        std::normal_distribution<double> norm(0.0, 1.0);

        std::vector<double> results(N);
        for (int j = 0; j < N; ++j) {
            double S = S0;
            for (int i = 0; i < steps; ++i) {
                double dW = sqrt_dt * norm(rng);
                S += mu * S * dt + sigma * S * dW;
            }
            results[j] = 1.0 / S;
        }
        return results;
    }

    std::vector<double> simulate_direct(double T, int steps, int N, std::mt19937 &rng) const {
        double dt = T / steps;
        double sqrt_dt = std::sqrt(dt);
        double mu_star = sigma * sigma - mu;
        std::normal_distribution<double> norm(0.0, 1.0);

        std::vector<double> results(N);
        for (int j = 0; j < N; ++j) {
            double U = 1.0 / S0;
            for (int i = 0; i < steps; ++i) {
                double dW = sqrt_dt * norm(rng);
                U += mu_star * U * dt + (-sigma) * U * dW;
            }
            results[j] = U;
        }
        return results;
    }

    static std::pair<double, double> moments(const std::vector<double> &v) {
        double mean = 0.0;
        for (double x: v)
            mean += x;
        mean /= v.size();

        double var = 0.0;
        for (double x: v)
            var += (x - mean) * (x - mean);
        var /= (v.size() - 1);

        return {mean, var};
    }
};

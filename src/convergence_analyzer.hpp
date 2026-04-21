#pragma once
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include "monte_carlo_engine.hpp"
#include "sde_base.hpp"

struct ConvergenceResult {
    std::vector<double> dt_values;
    std::vector<double> euler_errors;
    std::vector<double> milshtein_errors;
    std::vector<double> log_dt;
    std::vector<double> log_euler;
    std::vector<double> log_milshtein;
    double euler_slope = 0.0;
    double milshtein_slope = 0.0;
};

class ConvergenceAnalyzer {
public:
    ConvergenceAnalyzer(int num_paths = 2000, std::vector<int> steps_list = {8, 16, 32, 64, 128, 256, 512}) :
        num_paths_(num_paths), steps_list_(std::move(steps_list)) {}

    ConvergenceResult analyze_gbm(double mu, double sigma, double x0, double T) const {
        const std::string RST = "\033[0m";
        const std::string BOLD = "\033[1m";
        const std::string DIM = "\033[2m";
        const std::string CYAN = "\033[36m";
        const std::string GREEN = "\033[32m";
        const std::string YELLOW = "\033[33m";

        ConvergenceResult result;
        std::mt19937 gen(42);
        std::normal_distribution<double> nd(0.0, 1.0);

        std::cout << "\n  " << BOLD << std::left << std::setw(10) << "Steps" << std::setw(14) << "dt" << std::setw(18)
                  << "Euler error" << std::setw(18) << "Milshtein error" << RST << "\n"
                  << "  " << std::string(58, '-') << "\n";

        for (int N: steps_list_) {
            double dt = T / N;
            double sqrt_dt = std::sqrt(dt);

            double euler_err_sum = 0.0, milshtein_err_sum = 0.0;

            for (int path = 0; path < num_paths_; ++path) {
                std::vector<double> dW(N);
                double W_T = 0.0;
                for (int k = 0; k < N; ++k) {
                    dW[k] = nd(gen) * sqrt_dt;
                    W_T += dW[k];
                }
                double X_exact = x0 * std::exp((mu - 0.5 * sigma * sigma) * T + sigma * W_T);

                double X_euler = x0;
                for (int k = 0; k < N; ++k)
                    X_euler += mu * X_euler * dt + sigma * X_euler * dW[k];

                double X_milshtein = x0;
                for (int k = 0; k < N; ++k) {
                    double s = sigma * X_milshtein;
                    X_milshtein += mu * X_milshtein * dt + s * dW[k] + 0.5 * s * sigma * (dW[k] * dW[k] - dt);
                }

                euler_err_sum += std::abs(X_exact - X_euler);
                milshtein_err_sum += std::abs(X_exact - X_milshtein);
            }

            double euler_err = euler_err_sum / num_paths_;
            double milshtein_err = milshtein_err_sum / num_paths_;

            result.dt_values.push_back(dt);
            result.euler_errors.push_back(euler_err);
            result.milshtein_errors.push_back(milshtein_err);
            result.log_dt.push_back(std::log(dt));
            result.log_euler.push_back(std::log(euler_err));
            result.log_milshtein.push_back(std::log(milshtein_err));

            std::cout << "  " << DIM << std::left << std::setw(10) << N << std::scientific << std::setprecision(3)
                      << std::setw(14) << dt << RST << YELLOW << std::setw(18) << euler_err << RST << GREEN
                      << std::setw(18) << milshtein_err << RST << "\n";
        }

        result.euler_slope = linear_regression_slope(result.log_dt, result.log_euler);
        result.milshtein_slope = linear_regression_slope(result.log_dt, result.log_milshtein);

        return result;
    }

private:
    int num_paths_;
    std::vector<int> steps_list_;

    static double linear_regression_slope(const std::vector<double> &x, const std::vector<double> &y) {
        int n = static_cast<int>(x.size());
        if (n < 2)
            return 0.0;
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (int i = 0; i < n; ++i) {
            sx += x[i];
            sy += y[i];
            sxx += x[i] * x[i];
            sxy += x[i] * y[i];
        }
        return (n * sxy - sx * sy) / (n * sxx - sx * sx);
    }
};

#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <sstream>
#include "monte_carlo_engine.hpp"
#include "sde_base.hpp"
#include "console_utils.hpp"

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
        ConvergenceResult result;
        std::mt19937 gen(42);
        std::normal_distribution<double> nd(0.0, 1.0);

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
        }

        // Вывод таблицы
        std::cout << "\n";
        std::vector<std::string> headers = {"Steps", "dt", "Euler error", "Milshtein error"};
        std::vector<int> widths = {10, 14, 18, 18};
        print_table_header(headers, widths);

        for (size_t i = 0; i < steps_list_.size(); ++i) {
            std::ostringstream ss_steps, ss_dt, ss_eu, ss_mil;
            ss_steps << steps_list_[i];
            ss_dt << std::scientific << std::setprecision(3) << result.dt_values[i];
            ss_eu << result.euler_errors[i];
            ss_mil << result.milshtein_errors[i];
            print_table_row({ss_steps.str(), ss_dt.str(), ss_eu.str(), ss_mil.str()}, widths);
        }

        // Наклоны
        result.euler_slope = linear_regression_slope(result.log_dt, result.log_euler);
        result.milshtein_slope = linear_regression_slope(result.log_dt, result.log_milshtein);

        std::cout << "\n  log-log slopes\n";
        print_separator(44);
        std::cout << "  " << pad("Euler", 14) << std::fixed << std::setprecision(3)
                  << std::setw(8) << result.euler_slope << "   theory 0.5\n";
        std::cout << "  " << pad("Milshtein", 14) << std::setw(8) << result.milshtein_slope
                  << "   theory 1.0\n";

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
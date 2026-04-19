#pragma once
#include "sde_base.hpp"
#include "monte_carlo_engine.hpp"
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>

// =============================================================================
// ConvergenceAnalyzer — численное сравнение схем Эйлера и Мильштейна
//
// Метод: strong error (сильная ошибка)
//   e(dt) = E[ |X_T^exact - X_T^scheme| ]
//
// Для GBM точное решение известно:
//   X_T = X_0 * exp((mu - 0.5*sigma^2)*T + sigma*W_T)
//
// Ключевая идея: мы используем ОДНИ И ТЕ ЖЕ броуновские траектории
// для обеих схем и для точного решения. Это называется
// "pathwise convergence" и позволяет честно сравнить погрешности.
//
// Теоретические порядки сходимости:
//   Эйлер:     e(dt) ~ C * dt^{0.5}  => наклон log-log графика = 0.5
//   Мильштейн: e(dt) ~ C * dt^{1.0}  => наклон log-log графика = 1.0
// =============================================================================

struct ConvergenceResult {
    std::vector<double> dt_values;
    std::vector<double> euler_errors;
    std::vector<double> milshtein_errors;
    std::vector<double> log_dt;
    std::vector<double> log_euler;
    std::vector<double> log_milshtein;
    double euler_slope     = 0.0;  // должно быть ~0.5
    double milshtein_slope = 0.0;  // должно быть ~1.0
};

class ConvergenceAnalyzer {
public:
    // num_paths — сколько траекторий усреднять для оценки E[|error|]
    // steps_list — список количеств шагов (чем больше шагов, тем меньше dt)
    ConvergenceAnalyzer(int num_paths = 2000,
                        std::vector<int> steps_list = {8, 16, 32, 64, 128, 256, 512})
        : num_paths_(num_paths)
        , steps_list_(std::move(steps_list))
    {}

    // Анализ для GBM с известным точным решением
    // mu, sigma, x0, T — параметры модели
    ConvergenceResult analyze_gbm(double mu, double sigma, double x0, double T) const
    {
        ConvergenceResult result;
        std::mt19937 gen(42);  // фиксированный seed для воспроизводимости
        std::normal_distribution<double> nd(0.0, 1.0);

        std::cout << "\n=== Анализ сходимости (GBM, mu=" << mu
                  << ", sigma=" << sigma << ") ===\n";
        // Используем printf для корректного выравнивания с кириллицей
        std::printf("%-1s %8s %30s %36s\n",
                    "N шагов", "dt", "Ошибка Эйлер", "Ошибка Мильштейн");
        std::cout << std::string(66, '-') << "\n";

        for (int N : steps_list_) {
            double dt      = T / N;
            double sqrt_dt = std::sqrt(dt);

            double euler_err_sum     = 0.0;
            double milshtein_err_sum = 0.0;

            for (int path = 0; path < num_paths_; ++path) {
                // Генерируем броуновский путь с шагом dt
                // Одинаковый для обеих схем и для точного решения
                std::vector<double> dW(N);
                double W_T = 0.0;
                for (int k = 0; k < N; ++k) {
                    dW[k]  = nd(gen) * sqrt_dt;
                    W_T   += dW[k];
                }

                // Точное решение GBM
                double X_exact = x0 * std::exp(
                    (mu - 0.5 * sigma * sigma) * T + sigma * W_T
                );

                // Схема Эйлера
                double X_euler = x0;
                for (int k = 0; k < N; ++k) {
                    X_euler += mu * X_euler * dt + sigma * X_euler * dW[k];
                }

                // Схема Мильштейна
                // Для GBM: sigma(x) = sigma*x, sigma'(x) = sigma
                // Поправка: 0.5 * sigma*x * sigma * (dW^2 - dt)
                double X_milshtein = x0;
                for (int k = 0; k < N; ++k) {
                    double s = sigma * X_milshtein;
                    X_milshtein += mu * X_milshtein * dt
                                 + s * dW[k]
                                 + 0.5 * s * sigma * (dW[k]*dW[k] - dt);
                }

                euler_err_sum     += std::abs(X_exact - X_euler);
                milshtein_err_sum += std::abs(X_exact - X_milshtein);
            }

            double euler_err     = euler_err_sum     / num_paths_;
            double milshtein_err = milshtein_err_sum / num_paths_;

            result.dt_values.push_back(dt);
            result.euler_errors.push_back(euler_err);
            result.milshtein_errors.push_back(milshtein_err);
            result.log_dt.push_back(std::log(dt));
            result.log_euler.push_back(std::log(euler_err));
            result.log_milshtein.push_back(std::log(milshtein_err));

            std::printf("%-10d %-14.3e %-20.3e %-20.3e\n",
                        N, dt, euler_err, milshtein_err);
        }

        // Оцениваем наклон методом наименьших квадратов в log-log пространстве
        result.euler_slope     = linear_regression_slope(result.log_dt, result.log_euler);
        result.milshtein_slope = linear_regression_slope(result.log_dt, result.log_milshtein);

        std::cout << "\nНаклон log-log (теория: Эйлер=0.5, Мильштейн=1.0):\n";
        std::cout << "  Эйлер:      " << std::fixed << std::setprecision(3)
                  << result.euler_slope     << "\n";
        std::cout << "  Мильштейн:  " << result.milshtein_slope << "\n";

        return result;
    }

private:
    int              num_paths_;
    std::vector<int> steps_list_;

    // МНК для оценки наклона в log-log пространстве
    static double linear_regression_slope(const std::vector<double>& x,
                                          const std::vector<double>& y)
    {
        int n = static_cast<int>(x.size());
        if (n < 2) return 0.0;
        double sx  = 0, sy  = 0, sxx = 0, sxy = 0;
        for (int i = 0; i < n; ++i) {
            sx  += x[i]; sy  += y[i];
            sxx += x[i] * x[i];
            sxy += x[i] * y[i];
        }
        return (n * sxy - sx * sy) / (n * sxx - sx * sx);
    }
};
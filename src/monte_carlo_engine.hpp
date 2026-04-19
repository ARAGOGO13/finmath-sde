#pragma once
#include "sde_base.hpp"
#include <vector>
#include <random>
#include <cmath>
#include <numeric>
#include <algorithm>

// =============================================================================
// MonteCarloEngine — движок моделирования траекторий СДУ
//
// Реализует две схемы дискретизации (раздел C.3 PDF):
//
// ЭЙЛЕР (порядок сильной сходимости O(sqrt(dt))):
//   X_{k+1} = X_k + f(X_k)*dt + sigma(X_k)*dW
//
// МИЛЬШТЕЙН (порядок сильной сходимости O(dt)):
//   X_{k+1} = X_k + f(X_k)*dt + sigma(X_k)*dW
//            + 0.5 * sigma(X_k) * sigma'(X_k) * (dW^2 - dt)
//
// Поправка Мильштейна = 0.5 * sigma * sigma' * (dW^2 - dt)
// При sigma' = 0 (Vasicek) схемы совпадают.
// При sigma(x) = sigma*x (GBM) поправка существенна.
// =============================================================================

enum class Scheme { Euler, Milshtein };

class MonteCarloEngine {
public:
    MonteCarloEngine(size_t num_scenarios, int num_steps, double T,
                     Scheme scheme = Scheme::Euler)
        : num_scenarios_(num_scenarios)
        , num_steps_(num_steps)
        , dt_(T / num_steps)
        , scheme_(scheme)
        , gen_(std::random_device{}())
        , normal_(0.0, 1.0)
    {}

    // Основной метод: симулирует num_scenarios траекторий,
    // сохраняет первые num_paths_to_save для построения графиков.
    std::vector<std::vector<double>> simulate(
        const SDE& model, double x0, size_t num_paths_to_save = 50)
    {
        std::vector<std::vector<double>> paths(
            num_paths_to_save,
            std::vector<double>(num_steps_ + 1, 0.0)
        );

        double sum_XT = 0.0;
        double sum_XT2 = 0.0;

        for (size_t j = 0; j < num_scenarios_; ++j) {
            double X = x0;
            if (j < num_paths_to_save) paths[j][0] = X;

            for (int k = 0; k < num_steps_; ++k) {
                double t   = k * dt_;
                double dW  = normal_(gen_) * std::sqrt(dt_);
                double f   = model.drift(X, t);
                double s   = model.diffusion(X, t);

                X += f * dt_ + s * dW;

                // Поправка Мильштейна
                if (scheme_ == Scheme::Milshtein) {
                    double s_prime = model.diffusion_derivative(X, t);
                    X += 0.5 * s * s_prime * (dW * dW - dt_);
                }

                if (j < num_paths_to_save) paths[j][k + 1] = X;
            }
            sum_XT  += X;
            sum_XT2 += X * X;
        }

        mean_at_T_     = sum_XT  / static_cast<double>(num_scenarios_);
        double mean2   = sum_XT2 / static_cast<double>(num_scenarios_);
        variance_at_T_ = mean2 - mean_at_T_ * mean_at_T_;

        return paths;
    }

    double get_mean_at_T()     const { return mean_at_T_; }
    double get_variance_at_T() const { return variance_at_T_; }
    double get_std_at_T()      const { return std::sqrt(variance_at_T_); }

    Scheme get_scheme() const { return scheme_; }

private:
    size_t  num_scenarios_;
    int     num_steps_;
    double  dt_;
    Scheme  scheme_;

    std::mt19937                     gen_;
    std::normal_distribution<double> normal_;

    mutable double mean_at_T_     = 0.0;
    mutable double variance_at_T_ = 0.0;
};
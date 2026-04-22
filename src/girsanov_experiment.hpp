#pragma once
#include <vector>
#include "models/black_scholes.hpp"

// Численная проверка теоремы Гирсанова и анализ цены смены меры.
//
// Три оценщика одной цены колла при разных физических дрейфах mu:
//   1. BS-формула  (зависит только от r, не от mu)
//   2. MC под risk-neutral мерой Q
//   3. MC под физической мерой P с перевзвешиванием Радона-Никодима
//
// variance_ratio = (Var P-оценщика) / (Var Q-оценщика) ≈ exp(lambda^2*T)
// при больших |lambda| = |(mu-r)/sigma| демонстрирует экспоненциальный рост.
class GirsanovExperiment {
public:
    GirsanovExperiment(double r, double sigma, double T,
                       double K, double S0, int N_mc)
        : r_(r), sigma_(sigma), T_(T), K_(K), S0_(S0), N_mc_(N_mc) {}

    struct Point {
        double mu;           // физический дрейф
        double lambda;       // (mu - r) / sigma
        double bs_price;     // аналитическая BS-цена
        double rn_price;     // MC-цена под Q
        double rn_std_err;
        double p_price;      // MC-цена под P с перевзвешиванием
        double p_std_err;
        double variance_ratio; // (p_std_err / rn_std_err)^2
    };

    // Для каждого mu в mu_grid вычисляет все три оценки цены.
    // base_seed: для mu_i используется seed_i = base_seed + i,
    //            P-оценщик — seed_i + 1000 (чтобы иметь другой поток случайных чисел).
    std::vector<Point> run(const std::vector<double>& mu_grid,
                           unsigned base_seed = 42) const {
        BlackScholes bs(r_, sigma_, T_, K_, S0_);
        // BS-цена не зависит от mu — это концептуально важно
        const double analytical = bs.call_price();

        std::vector<Point> result;
        result.reserve(mu_grid.size());

        for (size_t i = 0; i < mu_grid.size(); ++i) {
            double   mu_i   = mu_grid[i];
            unsigned seed_i = base_seed + static_cast<unsigned>(i);

            auto rn = bs.price_call_mc_rn(N_mc_, seed_i);
            auto p  = bs.price_call_mc_girsanov(mu_i, N_mc_, seed_i + 1000u);

            double lambda     = (mu_i - r_) / sigma_;
            double var_ratio  = (rn.std_error > 1e-15)
                                ? (p.std_error / rn.std_error) * (p.std_error / rn.std_error)
                                : 1.0;

            result.push_back({mu_i, lambda,
                              analytical,
                              rn.price, rn.std_error,
                              p.price,  p.std_error,
                              var_ratio});
        }

        return result;
    }

private:
    double r_, sigma_, T_, K_, S0_;
    int    N_mc_;
};
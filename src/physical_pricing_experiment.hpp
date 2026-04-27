#pragma once

#include <cmath>
#include <vector>

#include "models/black_scholes.hpp"

// =============================================================================
// Эксперимент: ценообразование колла под физической мерой P vs риск-нейтральной Q.
//
// Наивный подход (без использования теоремы Гирсанова):
//   C0_P(μ) = e^{-rT} * E^P[max(S_T - K, 0)]
//           = e^{-rT} [S0 e^{μT} Φ(d1(μ)) - K Φ(d2(μ))]
//
// Риск-нейтральная цена (Блэк–Шоулз):
//   C0_Q = S0 Φ(d1(r)) - K e^{-rT} Φ(d2(r))
//
// Разность C0_P(μ) - C0_Q — «стоимость незнания» рыночной цены риска λ = (μ - r)/σ.
// При μ = r меры совпадают: C0_P = C0_Q.
// При μ > r наивный трейдер переоценивает опцион, при μ < r — недооценивает.
// =============================================================================

class PhysicalPricingExperiment {
public:
    struct Point {
        double mu;
        double lambda;
        double price_rn;
        double price_physical_analytic;
        double price_physical_mc;
        double mc_std_err;
        double difference;
        double difference_pct;
    };

    PhysicalPricingExperiment(double r, double sigma, double T,
                              double K, double S0, int N_mc = 200000)
        : r_(r), sigma_(sigma), T_(T), K_(K), S0_(S0), N_mc_(N_mc) {}

    // Запуск расчётов по сетке mu_grid. Возвращает вектор результатов для каждого μ.
    std::vector<Point> run(const std::vector<double>& mu_grid,
                           unsigned base_seed = 42) const {
        BlackScholes bs(r_, sigma_, T_, K_, S0_);
        const double price_rn = bs.call_price();

        std::vector<Point> result;
        result.reserve(mu_grid.size());

        for (size_t i = 0; i < mu_grid.size(); ++i) {
            double mu_i    = mu_grid[i];
            double price_an = bs.call_price_physical(mu_i);
            auto mc         = bs.price_call_mc_physical(
                mu_i, N_mc_, base_seed + static_cast<unsigned>(i));
            double diff     = price_an - price_rn;
            double diff_pct = diff / (std::abs(price_rn) + 1e-15) * 100.0;
            double lambda   = (mu_i - r_) / sigma_;
            result.push_back({mu_i, lambda, price_rn,
                              price_an, mc.price, mc.std_error,
                              diff, diff_pct});
        }
        return result;
    }

private:
    double r_, sigma_, T_, K_, S0_;
    int    N_mc_;
};

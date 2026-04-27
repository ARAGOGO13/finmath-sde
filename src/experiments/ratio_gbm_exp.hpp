// =============================================================================
// Эксперимент C_ratio: отношение двух коррелированных GBM (V_t = S_t / U_t).
//
// Параметры базовых процессов:
//   dS = μ_S S dt + σ_S S dW^(1),   dU = μ_U U dt + σ_U U dW^(2),
//   d⟨W^(1), W^(2)⟩ = ρ dt.
//
// Для каждого значения ρ из {-1.0, -0.5, 0.0, 0.5, 1.0}:
//   • вычисляется эффективный процесс V_t (дрейф μ_V, волатильность σ_V),
//   • генерируются траектории,
//   • сравниваются MC- и теоретические средние и дисперсии.
// Строятся графики траекторий с ±1σ и общей дисперсии для всех ρ.
// =============================================================================

#pragma once

#include <vector>
#include <string>
#include "../models/ratio_gbm.hpp"
#include "../plot_ratio.hpp"
#include "../console_utils.hpp"

inline void run_ratio_gbm_experiment() {
    using namespace std;

    section("V_t = S_t / U_t  (ratio of correlated GBMs)",
            "mu_V = mu_S - mu_U + sigma_U^2 - rho * sigma_S * sigma_U  |  sigma_V = sqrt(sigma_S^2 + sigma_U^2 - 2*rho*sigma_S*sigma_U)");

    const double T = 1.0;
    const double mu_S = 0.08, sigma_S = 0.20;
    const double mu_U = 0.05, sigma_U = 0.15;
    const double S0 = 100.0, U0 = 50.0;    // V0 = 2.0
    const int N_ratio_steps = 10;
    const int N_ratio_mc    = 200;
    const int N_ratio_plot  = 5;

    vector<double> rhos = {-1.0, -0.5, 0.0, 0.5, 1.0};

    auto tv2 = linspace(0.0, T, N_ratio_steps + 1);
    vector<double> tg2(tv2.begin(), tv2.end());

    RatioGBM model(mu_S, sigma_S, mu_U, sigma_U, S0, U0);

    ratio_table_header();

    vector<RatioGBMResult> results;
    for (size_t i = 0; i < rhos.size(); ++i) {
        auto res = model.simulate(rhos[i], tg2,
                                  N_ratio_mc, N_ratio_plot,
                                  42 + static_cast<unsigned>(i) * 7);
        double th_E = model.theory_mean_at(T, rhos[i]);
        ratio_table_row(res, th_E);
        results.push_back(std::move(res));
    }

    plot_ratio_paths(tg2, results, "../plots/C_ratio_paths");
    plot_ratio_all_vars(tg2, results, "../plots/C_ratio_var_all.pdf");

    print_saved({
        "C_ratio_paths_m10.pdf", "C_ratio_paths_m05.pdf",
        "C_ratio_paths_00.pdf",  "C_ratio_paths_p05.pdf",
        "C_ratio_paths_p10.pdf",
        "C_ratio_var_all.pdf"
    });
}
// =============================================================================
// Эксперимент C4 (часть 1): GBM и InverseGBM (U = 1/S).
//
// Симулируются траектории GBM с параметрами μ=0.05, σ=0.20, S0=100
// и обратного процесса U_t = 1/S_t. Сравниваются средние и дисперсии
// по Монте-Карло с теоретическими значениями. Строятся четыре графика:
//   E[S_t], Var[S_t] для GBM,
//   E[U_t], Var[U_t] для InverseGBM.
// =============================================================================

#pragma once

#include <vector>
#include <random>
#include <string>
#include <array>
#include "../plotting.hpp"
#include "../console_utils.hpp"

inline void run_gbm_inverse_experiment(const std::array<double,2>& shared_ylim) {
    using namespace std;

    section("GBM + InverseGBM (U=1/S)",
            "dS = mu S dt + sigma S dW  |  dU = (sigma^2 - mu) U dt - sigma U dW");

    const double T       = 1.0;
    const int    N_steps = 10;
    const int    N_mc    = 200;
    const int    N_paths = 5;

    const double mu = 0.05, sigma = 0.20, S0 = 100.0;
    const double U0 = 1.0 / S0, mu_star = sigma * sigma - mu;
    const double dt = T / N_steps, sqdt = sqrt(dt);

    auto tv = linspace(0.0, T, N_steps + 1);
    vector<double> time_std(tv.begin(), tv.end());

    mt19937 rng_g(42), rng_i(99);
    normal_distribution<double> nd(0.0, 1.0);

    vector<vector<double>> gbm_p(N_mc, vector<double>(N_steps + 1));
    vector<vector<double>> inv_p(N_mc, vector<double>(N_steps + 1));

    // Генерация траекторий GBM (схема Мильштейна)
    for (int j = 0; j < N_mc; ++j) {
        gbm_p[j][0] = S0;
        for (int k = 0; k < N_steps; ++k) {
            double dW = sqdt * nd(rng_g), S = gbm_p[j][k];
            gbm_p[j][k + 1] = S + mu * S * dt + sigma * S * dW
                            + 0.5 * sigma * S * sigma * (dW * dW - dt);
        }
    }
    // Генерация траекторий InverseGBM (схема Мильштейна)
    for (int j = 0; j < N_mc; ++j) {
        inv_p[j][0] = U0;
        for (int k = 0; k < N_steps; ++k) {
            double dW = sqdt * nd(rng_i), U = inv_p[j][k];
            inv_p[j][k + 1] = U + mu_star * U * dt + (-sigma) * U * dW
                            + 0.5 * sigma * sigma * U * (dW * dW - dt);
        }
    }

    // Теоретические моменты
    vector<double> gm(N_steps + 1), gv(N_steps + 1), im(N_steps + 1), iv(N_steps + 1);
    for (int i = 0; i <= N_steps; ++i) {
        double t = time_std[i];
        gm[i] = S0 * exp(mu * t);
        gv[i] = S0 * S0 * exp(2 * mu * t) * (exp(sigma * sigma * t) - 1.0);
        im[i] = U0 * exp(mu_star * t);
        iv[i] = U0 * U0 * exp(2 * mu_star * t) * (exp(sigma * sigma * t) - 1.0);
    }

    vector<vector<double>> gp(gbm_p.begin(), gbm_p.begin() + N_paths);
    vector<vector<double>> ip(inv_p.begin(), inv_p.begin() + N_paths);
    plot_gbm_and_inverse(time_std, gp, gm, gv, ip, im, iv,
        "dS=mu*S*dt+sigma*S*dW  [mu=0.05, sigma=0.20, S0=100]",
        "dU=(sigma^2-mu)*U*dt-sigma*U*dW  [mu*=-0.01, sigma=0.20, U0=0.01]",
        shared_ylim,
        "../plots/C4_gbm");

    vector<double> gmc_m, gmc_v, imc_m, imc_v;
    compute_mc_moments(gbm_p, gmc_m, gmc_v);
    compute_mc_moments(inv_p, imc_m, imc_v);

    table_header(18);
    table_row("GBM  E[S_T]",   gmc_m.back(), gm.back(), 18);
    table_row("GBM  Var[S_T]", gmc_v.back(), gv.back(), 18);
    table_row("Inv  E[U_T]",   imc_m.back(), im.back(), 18);
    table_row("Inv  Var[U_T]", imc_v.back(), iv.back(), 18);
    print_saved({"C4_gbm_mean.pdf", "C4_gbm_var.pdf",
                 "C4_gbm_inv_mean.pdf", "C4_gbm_inv_var.pdf"});
}
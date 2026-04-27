// =============================================================================
// Эксперимент C4 (часть 3): модель CIR.
//
// Симулируются траектории процесса:
//   dX = θ (μ - X) dt + σ √X dW
// Параметры: θ=2.0, μ=0.05, σ=0.02, X0=0.03.
// Проверяется условие Феллера. Сравниваются средние и дисперсии
// по Монте-Карло с теоретическими значениями.
// =============================================================================

#pragma once

#include <vector>
#include <random>
#include <string>
#include <array>
#include "../models/cir.hpp"
#include "../models/ornstein_uhlenbeck.hpp"
#include "../plotting.hpp"
#include "../console_utils.hpp"

inline void run_cir_experiment(const std::array<double,2>& shared_ylim) {
    using namespace std;

    section("Cox-Ingersoll-Ross",
            "dX = theta (mu - X) dt + sigma sqrt(X) dW  |  Feller: 2*theta*mu >= sigma^2");

    const double T       = 1.0;
    const int    N_steps = 10;
    const int    N_mc    = 200;
    const int    N_paths = 5;

    const double theta = 2.0, mu = 0.05, sigma = 0.02, x0 = 0.03;
    CoxIngersollRoss cir(theta, mu, sigma);
    bool feller = cir.feller_condition();
    cout << "  Feller condition: 2*theta*mu = " << 2*theta*mu
         << ", sigma^2 = " << sigma*sigma
         << " -> " << (feller ? "satisfied" : "violated") << "\n\n";

    auto tv = linspace(0.0, T, N_steps + 1);
    vector<double> time_std(tv.begin(), tv.end());

    mt19937 rng(321);
    normal_distribution<double> nd(0.0, 1.0);
    const double dt = T / N_steps, sqdt = sqrt(dt);

    vector<vector<double>> all(N_mc, vector<double>(N_steps + 1));
    for (int j = 0; j < N_mc; ++j) {
        all[j][0] = x0;
        for (int k = 0; k < N_steps; ++k) {
            double dW = sqdt * nd(rng);
            double X = std::max(all[j][k], 0.0);
            double s = cir.diffusion(X, 0), sp = cir.diffusion_derivative(X, 0);
            all[j][k + 1] = X + cir.drift(X, 0) * dt + s * dW + 0.5 * s * sp * (dW * dW - dt);
        }
    }

    OrnsteinUhlenbeck ou_for_mean(theta, mu, sigma);
    vector<double> mth(N_steps + 1), vth(N_steps + 1);
    for (int i = 0; i <= N_steps; ++i) {
        mth[i] = ou_for_mean.exact_mean(x0, time_std[i]);
        vth[i] = cir_exact_variance(x0, time_std[i], theta, mu, sigma);
    }

    vector<vector<double>> pp(all.begin(), all.begin() + N_paths);
    plot_mean_reversion(time_std, pp, mth, vth, "CIR",
        "dX=theta*(mu-X)*dt+sigma*sqrt(X)*dW  [theta=2.0, mu=0.05, sigma=0.02, X0=0.03]",
        shared_ylim, "../plots/C4_cir");

    vector<double> mc_m, mc_v;
    compute_mc_moments(all, mc_m, mc_v);
    table_header(18);
    table_row("E[X_T]",   mc_m.back(), mth.back(), 18);
    table_row("Var[X_T]", mc_v.back(), vth.back(), 18);
    print_saved({"C4_cir_mean.pdf", "C4_cir_var.pdf"});
}
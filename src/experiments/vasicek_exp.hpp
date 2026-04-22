#pragma once

#include <vector>
#include <random>
#include <string>
#include <array>
#include "../models/ornstein_uhlenbeck.hpp"
#include "../plotting.hpp"
#include "../console_utils.hpp"

// =============================================================================
// Эксперимент: модель Васичека (Орнштейна-Уленбека)
// =============================================================================
inline void run_vasicek_experiment(const std::array<double,2>& shared_ylim) {
    using namespace std;

    section("Vasicek / Ornstein-Uhlenbeck",
            "dX = theta (mu - X) dt + sigma dW  |  sigma' = 0 => Euler = Milshtein");

    const double T       = 1.0;
    const int    N_steps = 50;
    const int    N_mc    = 2000;
    const int    N_paths = 50;

    const double theta = 2.0, mu = 0.05, sigma = 0.02, x0 = 0.03;
    OrnsteinUhlenbeck vasicek(theta, mu, sigma);

    auto tv = linspace(0.0, T, N_steps + 1);
    vector<double> time_std(tv.begin(), tv.end());

    mt19937 rng(123);
    normal_distribution<double> nd(0.0, 1.0);
    const double dt = T / N_steps, sqdt = sqrt(dt);

    vector<vector<double>> all(N_mc, vector<double>(N_steps + 1));
    for (int j = 0; j < N_mc; ++j) {
        all[j][0] = x0;
        for (int k = 0; k < N_steps; ++k) {
            double dW = sqdt * nd(rng), X = all[j][k];
            all[j][k + 1] = X + vasicek.drift(X, 0) * dt + vasicek.diffusion(X, 0) * dW;
        }
    }

    vector<double> mth(N_steps + 1), vth(N_steps + 1);
    for (int i = 0; i <= N_steps; ++i) {
        mth[i] = vasicek.exact_mean(x0, time_std[i]);
        vth[i] = vasicek.exact_var(time_std[i]);
    }

    vector<vector<double>> pp(all.begin(), all.begin() + N_paths);
    plot_mean_reversion(time_std, pp, mth, vth, "Vasicek",
        "dX=theta*(mu-X)*dt+sigma*dW  [theta=2.0, mu=0.05, sigma=0.02, X0=0.03]",
        shared_ylim, "../plots/C4_vasicek");

    vector<double> mc_m, mc_v;
    compute_mc_moments(all, mc_m, mc_v);
    table_header(18);
    table_row("E[X_T]",   mc_m.back(), mth.back(), 18);
    table_row("Var[X_T]", mc_v.back(), vth.back(), 18);
    print_saved({"C4_vasicek_mean.pdf", "C4_vasicek_var.pdf"});
}
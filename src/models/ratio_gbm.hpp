#pragma once
#include <vector>
#include <random>
#include <cmath>
#include <string>

// =============================================================================
// RatioGBM — V_t = S_t / U_t
//
// dS = mu_S * S * dt + sigma_S * S * dW1
// dU = mu_U * U * dt + sigma_U * U * dW2
// dW1 * dW2 = rho * dt
//
// По формуле Ито (задача 3):
//   dV = mu_V * V * dt + sigma_V * V * dW
//
// где:
//   mu_V    = mu_S - mu_U + sigma_U^2 - rho * sigma_S * sigma_U
//   sigma_V = sqrt(sigma_S^2 + sigma_U^2 - 2*rho*sigma_S*sigma_U)
//
// Коррелированные BM генерируются разложением Холецкого:
//   dW1 = sqrt(dt) * Z1
//   dW2 = sqrt(dt) * (rho * Z1 + sqrt(1-rho^2) * Z2)
// =============================================================================

struct RatioGBMResult {
    double rho;
    double mu_V;
    double sigma_V;

    std::vector<std::vector<double>> paths;

    std::vector<double> mc_mean;
    std::vector<double> mc_var;
    std::vector<double> theory_mean;
    std::vector<double> theory_var;
};

class RatioGBM {
public:
    double mu_S, sigma_S;
    double mu_U, sigma_U;
    double V0;   // S0 / U0

    RatioGBM(double mu_S, double sigma_S,
             double mu_U, double sigma_U,
             double S0, double U0)
        : mu_S(mu_S), sigma_S(sigma_S)
        , mu_U(mu_U), sigma_U(sigma_U)
        , V0(S0 / U0) {}

    // Теоретический эффективный дрейф V_t
    double theory_mu_V(double rho) const {
        return mu_S - mu_U + sigma_U*sigma_U - rho*sigma_S*sigma_U;
    }

    // Теоретическая эффективная волатильность V_t
    double theory_sigma_V(double rho) const {
        double s2 = sigma_S*sigma_S + sigma_U*sigma_U - 2*rho*sigma_S*sigma_U;
        return std::sqrt(std::max(s2, 0.0));
    }

    // E[V_t] = V0 * exp(mu_V * t)
    double theory_mean_at(double t, double rho) const {
        return V0 * std::exp(theory_mu_V(rho) * t);
    }

    // Var[V_t] = V0^2 * exp(2*mu_V*t) * (exp(sigma_V^2*t) - 1)
    double theory_var_at(double t, double rho) const {
        double mu_v  = theory_mu_V(rho);
        double sig_v = theory_sigma_V(rho);
        return V0*V0 * std::exp(2*mu_v*t) * (std::exp(sig_v*sig_v*t) - 1.0);
    }

    RatioGBMResult simulate(
        double rho,
        const std::vector<double>& time_grid,
        int N_mc,
        int N_plot_paths,
        unsigned seed = 42) const
    {
        RatioGBMResult res;
        res.rho     = rho;
        res.mu_V    = theory_mu_V(rho);
        res.sigma_V = theory_sigma_V(rho);

        size_t T_steps = time_grid.size();
        double T  = time_grid.back();
        int    N  = static_cast<int>(T_steps) - 1;
        double dt = T / N;
        double sqdt = std::sqrt(dt);


        std::mt19937 rng(seed);
        std::normal_distribution<double> nd(0.0, 1.0);
        double rho_c = std::sqrt(std::max(1.0 - rho*rho, 0.0));

        std::vector<std::vector<double>> all_paths(N_mc,
            std::vector<double>(T_steps));

        for (int j = 0; j < N_mc; ++j) {
            all_paths[j][0] = V0;
            for (int k = 0; k < N; ++k) {
                double Z1 = nd(rng), Z2 = nd(rng);
                double dW = sqdt * Z1;
                double V  = all_paths[j][k];
                // Milshtein
                all_paths[j][k+1] = V
                    + res.mu_V    * V * dt
                    + res.sigma_V * V * dW
                    + 0.5 * res.sigma_V * res.sigma_V * V * (dW*dW - dt);
            }
        }

        res.mc_mean.assign(T_steps, 0.0);
        res.mc_var .assign(T_steps, 0.0);
        for (size_t t = 0; t < T_steps; ++t) {
            double s = 0.0, s2 = 0.0;
            for (int j = 0; j < N_mc; ++j) {
                s  += all_paths[j][t];
                s2 += all_paths[j][t] * all_paths[j][t];
            }
            res.mc_mean[t] = s  / N_mc;
            res.mc_var [t] = s2 / N_mc - res.mc_mean[t] * res.mc_mean[t];
        }

        res.theory_mean.resize(T_steps);
        res.theory_var .resize(T_steps);
        for (size_t i = 0; i < T_steps; ++i) {
            res.theory_mean[i] = theory_mean_at(time_grid[i], rho);
            res.theory_var [i] = theory_var_at (time_grid[i], rho);
        }




        int save = std::min(N_plot_paths, N_mc);
        res.paths.assign(all_paths.begin(), all_paths.begin() + save);

        return res;
    }
};
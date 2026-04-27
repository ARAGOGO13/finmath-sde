// =============================================================================
// Модель отношения двух коррелированных геометрических броуновских движений.
//
// Пусть S_t и U_t следуют GBM с корреляцией ρ:
//   dS_t = μ_S S_t dt + σ_S S_t dW_t^(1)
//   dU_t = μ_U U_t dt + σ_U U_t dW_t^(2)
//   d⟨W^(1), W^(2)⟩_t = ρ dt
//
// Процесс V_t = S_t / U_t также является GBM:
//   dV_t = μ_V V_t dt + σ_V V_t dW_t
// где μ_V = μ_S - μ_U + σ_U² - ρ σ_S σ_U,
//      σ_V² = σ_S² + σ_U² - 2 ρ σ_S σ_U.
//
// Класс предоставляет аналитические моменты V_t и метод Монте-Карло для
// симуляции траекторий с заданной корреляцией.
// =============================================================================

#pragma once
#include <vector>
#include <random>
#include <cmath>
#include <string>

// Результат моделирования отношения двух GBM
struct RatioGBMResult {
    double rho;                     // коэффициент корреляции
    double mu_V;                    // эффективный дрейф V
    double sigma_V;                 // эффективная волатильность V

    std::vector<std::vector<double>> paths;   // сохранённые траектории для визуализации

    std::vector<double> mc_mean;    // Монте-Карло оценка E[V_t] по шагам времени
    std::vector<double> mc_var;     // Монте-Карло оценка Var[V_t] по шагам времени
    std::vector<double> theory_mean;// теоретическое E[V_t] по шагам времени
    std::vector<double> theory_var; // теоретическое Var[V_t] по шагам времени
};

class RatioGBM {
public:
    double mu_S, sigma_S;   // параметры процесса S
    double mu_U, sigma_U;   // параметры процесса U
    double V0;              // начальное значение V0 = S0 / U0

    RatioGBM(double mu_S, double sigma_S,
             double mu_U, double sigma_U,
             double S0, double U0)
        : mu_S(mu_S), sigma_S(sigma_S)
        , mu_U(mu_U), sigma_U(sigma_U)
        , V0(S0 / U0) {}

    // Эффективный дрейф процесса V_t при заданной корреляции ρ
    double theory_mu_V(double rho) const {
        return mu_S - mu_U + sigma_U * sigma_U - rho * sigma_S * sigma_U;
    }

    // Эффективная волатильность процесса V_t при заданной корреляции ρ
    double theory_sigma_V(double rho) const {
        double s2 = sigma_S * sigma_S + sigma_U * sigma_U - 2 * rho * sigma_S * sigma_U;
        return std::sqrt(std::max(s2, 0.0));
    }

    // Теоретическое математическое ожидание V_t в момент времени t
    double theory_mean_at(double t, double rho) const {
        return V0 * std::exp(theory_mu_V(rho) * t);
    }

    // Теоретическая дисперсия V_t в момент времени t
    double theory_var_at(double t, double rho) const {
        double mu_v  = theory_mu_V(rho);
        double sig_v = theory_sigma_V(rho);
        return V0 * V0 * std::exp(2 * mu_v * t) * (std::exp(sig_v * sig_v * t) - 1.0);
    }

    // Монте-Карло симуляция траекторий V_t с использованием схемы Мильштейна.
    // rho            – корреляция между броуновскими движениями S и U
    // time_grid      – сетка моментов времени, включая 0 и T
    // N_mc           – число сценариев для оценки моментов
    // N_plot_paths   – число траекторий, сохраняемых для визуализации
    // seed           – зерно генератора псевдослучайных чисел
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
                // Схема Мильштейна для процесса V_t
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
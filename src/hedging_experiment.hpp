// =============================================================================
// Эксперимент по дискретному дельта-хеджированию короткой позиции в европейском колле.
//
// Продавец получает премию C0 (BS-цена с подразумеваемой волатильностью σ_impl),
// затем ребалансирует портфель на равномерной сетке времени, поддерживая
// дельта-нейтральность путём покупки/продажи базового актива.
//
// Основные результаты (Bertsimas, Kogan, Lo 1998):
//   - При σ_real = σ_impl стандартное отклонение P&L ~ √(Δt).
//   - При σ_real ≠ σ_impl возникает арбитраж волатильности:
//     E[P&L] ≈ 0.5 · Γ · S0² · (σ_impl² - σ_real²) · T.
// =============================================================================

#pragma once

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include "models/black_scholes.hpp"

// Параметры эксперимента по хеджированию
struct HedgingConfig {
    double r          = 0.05;   // безрисковая процентная ставка
    double T          = 1.0;    // время до экспирации (лет)
    double K          = 100.0;  // страйк опциона
    double S0         = 100.0;  // начальная цена базового актива
    double sigma_real = 0.20;   // реальная волатильность для симуляции цены
    double sigma_impl = 0.20;   // подразумеваемая волатильность для расчёта премии и дельты
};

// Результат одного эксперимента по дельта-хеджированию
struct HedgingPnLResult {
    int    N_rebal;                      // число ребалансировок
    double dt;                           // шаг времени Δt = T / N_rebal
    std::vector<double> pnl_samples;     // выборка значений P&L
    double mean_pnl;                     // среднее P&L
    double std_pnl;                      // стандартное отклонение P&L
    double var_pnl;                      // дисперсия P&L
    double q05, q95;                     // 5% и 95% квантили распределения P&L
};

class DeltaHedgingExperiment {
public:
    explicit DeltaHedgingExperiment(HedgingConfig cfg = {}) : cfg_(cfg) {}

    // Симуляция дельта-хеджирования с заданным числом ребалансировок.
    // N_rebal – число шагов; N_paths – число сценариев Монте-Карло;
    // seed – зерно генератора.
    HedgingPnLResult simulate(int N_rebal, int N_paths, unsigned seed) const {
        if (N_rebal < 1) N_rebal = 1;

        const double dt        = cfg_.T / N_rebal;
        const double exp_r_dt  = std::exp(cfg_.r * dt);
        const double drift_adj = (cfg_.r - 0.5 * cfg_.sigma_real * cfg_.sigma_real) * dt;
        const double vol_sqdt  = cfg_.sigma_real * std::sqrt(dt);

        // Вспомогательная функция для вычисления дельты Блэка–Шоулза
        auto bs_delta = [&](double S, int step_k) -> double {
            double tau = cfg_.T - step_k * dt;
            if (tau < 1e-12) return (S > cfg_.K) ? 1.0 : 0.0;
            double sqrt_tau = std::sqrt(tau);
            double d1 = (std::log(S / cfg_.K) + (cfg_.r + 0.5*cfg_.sigma_impl*cfg_.sigma_impl)*tau)
                       / (cfg_.sigma_impl * sqrt_tau);
            return 0.5 * std::erfc(-d1 / std::sqrt(2.0));
        };

        const double delta0 = bs_delta(cfg_.S0, 0);

        BlackScholes bs0(cfg_.r, cfg_.sigma_impl, cfg_.T, cfg_.K, cfg_.S0);
        const double C0 = bs0.call_price();

        std::mt19937_64 rng(seed);
        std::normal_distribution<double> nd(0.0, 1.0);

        std::vector<double> samples(N_paths);

        for (int path = 0; path < N_paths; ++path) {
            double S     = cfg_.S0;
            double delta = delta0;
            double B     = C0 - delta * S;   // начальный банковский счёт

            for (int k = 0; k < N_rebal; ++k) {
                double S_new = S * std::exp(drift_adj + vol_sqdt * nd(rng));
                double B_new = B * exp_r_dt;

                if (k == N_rebal - 1) {
                    // Последний шаг: закрытие позиции
                    samples[path] = B_new + delta * S_new - std::max(S_new - cfg_.K, 0.0);
                } else {
                    double delta_new = bs_delta(S_new, k + 1);
                    B_new -= (delta_new - delta) * S_new;  // корректировка позиции
                    S     = S_new;
                    B     = B_new;
                    delta = delta_new;
                }
            }
        }

        double sum = 0.0, sum2 = 0.0;
        for (double v : samples) { sum += v; sum2 += v * v; }
        double mean = sum / N_paths;
        double var  = std::max(sum2 / N_paths - mean * mean, 0.0);

        std::vector<double> sorted = samples;
        std::sort(sorted.begin(), sorted.end());
        int i05 = static_cast<int>(0.05 * (N_paths - 1));
        int i95 = static_cast<int>(0.95 * (N_paths - 1));

        return HedgingPnLResult{
            N_rebal, dt, std::move(samples),
            mean, std::sqrt(var), var,
            sorted[i05], sorted[i95]
        };
    }

    // Сканирование по частоте ребалансировки.
    // Оценивает зависимость std(P&L) от Δt и наклон в log-log координатах.
    struct ScanResult {
        std::vector<int>    N_rebal_list;
        std::vector<double> dt_list;
        std::vector<double> mean_pnl;
        std::vector<double> std_pnl;
        double              slope_log_std_vs_log_dt; // наклон в log-log
    };

    ScanResult scan_rebalance_frequency(const std::vector<int>& N_rebal_list,
                                        int N_paths,
                                        unsigned base_seed) const {
        ScanResult result;
        result.N_rebal_list = N_rebal_list;

        for (size_t i = 0; i < N_rebal_list.size(); ++i) {
            auto r = simulate(N_rebal_list[i], N_paths, base_seed + static_cast<unsigned>(i) * 17u);
            result.dt_list .push_back(r.dt);
            result.mean_pnl.push_back(r.mean_pnl);
            result.std_pnl .push_back(r.std_pnl);
        }

        const int n = static_cast<int>(result.dt_list.size());
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (int i = 0; i < n; ++i) {
            double xi = std::log(result.dt_list[i]);
            double yi = std::log(std::max(result.std_pnl[i], 1e-15));
            sx += xi; sy += yi; sxx += xi*xi; sxy += xi*yi;
        }
        result.slope_log_std_vs_log_dt = (n*sxy - sx*sy) / (n*sxx - sx*sx);

        return result;
    }

private:
    HedgingConfig cfg_;
};
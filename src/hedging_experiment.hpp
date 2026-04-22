#pragma once
// =============================================================================
// Дискретный delta-hedging P&L эксперимент.
//
// Продавец европейского колла получает премию C_0 = BS(S_0, sigma_impl) и
// поддерживает delta-нейтральный портфель, ребалансируя на сетке 0=t_0<...<t_N=T.
//
// Ключевые результаты (Bertsimas, Kogan, Lo 1998, JFE):
//   (A) std(P&L) ~ sqrt(dt) = sqrt(T/N) при sigma_real = sigma_impl
//       [наклон +0.5 в log-log по N]
//   (B) E[P&L] ≈ 0 при sigma_real == sigma_impl
//   (C) E[P&L] ≈ 0.5 * Γ * S_0^2 * (sigma_impl^2 - sigma_real^2) * T
//       при несовпадении волатильностей ("vol arbitrage"):
//       sigma_real > sigma_impl → продавец теряет (E[P&L] < 0).
//
// Bertsimas D., Kogan L., Lo A.W. (1998) "When Is Time Continuous?"
//   Journal of Financial Economics 55(2), 173-204.
// =============================================================================
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include "models/black_scholes.hpp"

struct HedgingConfig {
    double r          = 0.05;
    double T          = 1.0;
    double K          = 100.0;
    double S0         = 100.0;
    double sigma_real = 0.20; // реальная vol для симуляции движения
    double sigma_impl = 0.20; // implied vol для хеджирования и расчёта премии
};

struct HedgingPnLResult {
    int    N_rebal;
    double dt;
    std::vector<double> pnl_samples;
    double mean_pnl;
    double std_pnl;
    double var_pnl;
    double q05, q95; // квантили для анализа хвостов
};

class DeltaHedgingExperiment {
public:
    explicit DeltaHedgingExperiment(HedgingConfig cfg = {}) : cfg_(cfg) {}

    // Прогоняет N_paths путей при данном числе ребалансировок.
    // N_rebal — число временных шагов (N_rebal-1 промежуточных ребалансировок).
    HedgingPnLResult simulate(int N_rebal, int N_paths, unsigned seed) const {
        if (N_rebal < 1) N_rebal = 1;

        const double dt        = cfg_.T / N_rebal;
        const double exp_r_dt  = std::exp(cfg_.r * dt);
        const double drift_adj = (cfg_.r - 0.5 * cfg_.sigma_real * cfg_.sigma_real) * dt;
        const double vol_sqdt  = cfg_.sigma_real * std::sqrt(dt);

        // Предвычислим для каждого шага k: tau_k = T-(k+1)*dt и производные
        // (используются для BS-дельты на промежуточных шагах k=0..N_rebal-2)
        std::vector<double> adj_k(N_rebal - 1), inv_vol_tau_k(N_rebal - 1);
        for (int k = 0; k < N_rebal - 1; ++k) {
            double tau = cfg_.T - (k + 1) * dt;
            if (tau < 1e-12) tau = 1e-12;
            double sqrt_tau        = std::sqrt(tau);
            adj_k[k]              = (cfg_.r + 0.5 * cfg_.sigma_impl * cfg_.sigma_impl) * tau;
            inv_vol_tau_k[k]      = 1.0 / (cfg_.sigma_impl * sqrt_tau);
        }

        // BS-дельта по sigma_impl, tau = оставшееся время
        auto bs_delta = [&](double S, int step_k) -> double {
            // step_k — индекс шага ДО которого мы вычисляем дельту
            // tau = T - step_k * dt (оставшееся время на начало шага step_k)
            double tau = cfg_.T - step_k * dt;
            if (tau < 1e-12) return (S > cfg_.K) ? 1.0 : 0.0;
            double sqrt_tau = std::sqrt(tau);
            double d1 = (std::log(S / cfg_.K) + (cfg_.r + 0.5*cfg_.sigma_impl*cfg_.sigma_impl)*tau)
                       / (cfg_.sigma_impl * sqrt_tau);
            return 0.5 * std::erfc(-d1 * k_inv_sqrt2);
        };

        // Дельта в t=0 (tau=T)
        const double delta0 = bs_delta(cfg_.S0, 0);

        // BS-цена по sigma_impl в t=0
        BlackScholes bs0(cfg_.r, cfg_.sigma_impl, cfg_.T, cfg_.K, cfg_.S0);
        const double C0 = bs0.call_price();

        std::mt19937_64 rng(seed);
        std::normal_distribution<double> nd(0.0, 1.0);

        std::vector<double> samples(N_paths);

        for (int path = 0; path < N_paths; ++path) {
            double S     = cfg_.S0;
            double delta = delta0;
            double B     = C0 - delta * S; // счёт в банке

            for (int k = 0; k < N_rebal; ++k) {
                // Шаг a: симулируем S_{k+1} под Q с sigma_real
                double S_new = S * std::exp(drift_adj + vol_sqdt * nd(rng));
                // Шаг b: банк растёт по ставке r
                double B_new = B * exp_r_dt;

                if (k == N_rebal - 1) {
                    // Шаг e: последняя итерация — закрываем позицию
                    // P&L = банк + выручка от продажи акций - выплата по коллу
                    samples[path] = B_new + delta * S_new - std::max(S_new - cfg_.K, 0.0);
                } else {
                    // Шаги c-d: ребалансировка
                    double delta_new = bs_delta(S_new, k + 1);
                    B_new -= (delta_new - delta) * S_new;
                    S     = S_new;
                    B     = B_new;
                    delta = delta_new;
                }
            }
        }

        // Статистики
        double sum = 0.0, sum2 = 0.0;
        for (double v : samples) { sum += v; sum2 += v * v; }
        double mean = sum / N_paths;
        double var  = std::max(sum2 / N_paths - mean * mean, 0.0);

        // Квантили
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

    // Сканирует набор значений N_rebal, возвращает std(P&L) vs dt для log-log анализа.
    struct ScanResult {
        std::vector<int>    N_rebal_list;
        std::vector<double> dt_list;
        std::vector<double> mean_pnl;
        std::vector<double> std_pnl;
        double              slope_log_std_vs_log_dt; // теоретически ≈ +0.5
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

        // OLS slope: log(std_pnl) ~ slope * log(dt) + const
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

    static constexpr double k_inv_sqrt2 = 0.7071067811865475244; // 1/sqrt(2)
};
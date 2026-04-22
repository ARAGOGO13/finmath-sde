#pragma once

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include "models/black_scholes.hpp"

// =============================================================================
// Эксперимент по дискретному дельта-хеджированию короткой позиции в европейском колле
//
// Класс DeltaHedgingExperiment моделирует стратегию дельта-нейтрального
// хеджирования для продавца европейского колл-опциона. Продавец получает
// премию C₀, рассчитанную по формуле Блэка–Шоулза с подразумеваемой
// волатильностью σ_impl, и затем ребалансирует портфель на равномерной
// временной сетке 0 = t₀ < t₁ < ... < t_N = T, поддерживая дельта-нейтральность
// путём покупки/продажи базового актива.
//
// Ключевые теоретические результаты (Bertsimas, Kogan, Lo 1998):
//   (A) При совпадении реальной и подразумеваемой волатильностей (σ_real = σ_impl)
//       стандартное отклонение P&L ведёт себя как √(Δt) = √(T/N), что
//       соответствует наклону +0.5 в логарифмическом масштабе.
//   (B) Математическое ожидание P&L приблизительно равно нулю при σ_real = σ_impl.
//   (C) При несовпадении волатильностей возникает арбитраж волатильности:
//       E[P&L] ≈ 0.5 · Γ · S₀² · (σ_impl² − σ_real²) · T,
//       где Γ – гамма опциона. Если σ_real > σ_impl, продавец в среднем несёт убытки.
//
// Эксперимент позволяет исследовать влияние частоты ребалансировки (N) и
// рассогласования волатильностей на распределение итогового P&L.
// =============================================================================

struct HedgingConfig {
    double r          = 0.05;   // безрисковая процентная ставка
    double T          = 1.0;    // время до экспирации (в годах)
    double K          = 100.0;  // страйк опциона
    double S0         = 100.0;  // начальная цена базового актива
    double sigma_real = 0.20;   // реальная волатильность для симуляции цены актива
    double sigma_impl = 0.20;   // подразумеваемая волатильность для расчёта премии и дельты
};

// -----------------------------------------------------------------------------
// Результат одного эксперимента по дельта-хеджированию
// -----------------------------------------------------------------------------
struct HedgingPnLResult {
    int    N_rebal;                      // число временных шагов (ребалансировок)
    double dt;                           // длина шага Δt = T / N_rebal
    std::vector<double> pnl_samples;     // вектор значений P&L для всех сценариев
    double mean_pnl;                     // выборочное среднее P&L
    double std_pnl;                      // выборочное стандартное отклонение P&L
    double var_pnl;                      // выборочная дисперсия P&L
    double q05, q95;                     // 5% и 95% квантили распределения P&L
};

// -----------------------------------------------------------------------------
// Основной класс эксперимента
// -----------------------------------------------------------------------------
class DeltaHedgingExperiment {
public:
    explicit DeltaHedgingExperiment(HedgingConfig cfg = {}) : cfg_(cfg) {}

    // -------------------------------------------------------------------------
    // Симуляция дельта-хеджирования с заданным числом ребалансировок
    // -------------------------------------------------------------------------
    // Параметры:
    //   N_rebal – число временных шагов (N_rebal ≥ 1). При N_rebal = 1
    //             хеджирование отсутствует (только начальная дельта).
    //   N_paths – число сценариев Монте-Карло
    //   seed    – зерно генератора псевдослучайных чисел
    //
    // Алгоритм на одном сценарии:
    //   1. В момент t=0 вычисляется BS-дельта (по σ_impl) и начальный банк:
    //      B₀ = C₀ − Δ₀·S₀.
    //   2. Для каждого шага k = 0, 1, …, N_rebal−1:
    //      a) Генерируется цена актива S_{k+1} под физической мерой
    //         с реальной волатильностью σ_real и дрейфом r.
    //      b) Банк растёт по безрисковой ставке: B ← B·exp(r·Δt).
    //      c) Если k = N_rebal−1 (последний шаг):
    //         – позиция закрывается: продаются оставшиеся акции Δ_{last}·S_T,
    //           выплачивается опционный платёж max(S_T − K, 0).
    //         – P&L = банк + выручка от продажи акций − выплата.
    //      d) Иначе (промежуточная ребалансировка):
    //         – вычисляется новая дельта Δ_new (по σ_impl).
    //         – корректируется банк: B ← B − (Δ_new − Δ_old)·S_new.
    //         – обновляются состояние S и дельта Δ.
    //
    // Возвращает структуру HedgingPnLResult с выборкой P&L и статистиками.
    // -------------------------------------------------------------------------
    HedgingPnLResult simulate(int N_rebal, int N_paths, unsigned seed) const {
        if (N_rebal < 1) N_rebal = 1;

        const double dt        = cfg_.T / N_rebal;
        const double exp_r_dt  = std::exp(cfg_.r * dt);
        const double drift_adj = (cfg_.r - 0.5 * cfg_.sigma_real * cfg_.sigma_real) * dt;
        const double vol_sqdt  = cfg_.sigma_real * std::sqrt(dt);

        std::vector<double> adj_k(N_rebal - 1), inv_vol_tau_k(N_rebal - 1);
        for (int k = 0; k < N_rebal - 1; ++k) {
            double tau = cfg_.T - (k + 1) * dt;
            if (tau < 1e-12) tau = 1e-12;
            double sqrt_tau        = std::sqrt(tau);
            adj_k[k]              = (cfg_.r + 0.5 * cfg_.sigma_impl * cfg_.sigma_impl) * tau;
            inv_vol_tau_k[k]      = 1.0 / (cfg_.sigma_impl * sqrt_tau);
        }

        auto bs_delta = [&](double S, int step_k) -> double {
            double tau = cfg_.T - step_k * dt;
            if (tau < 1e-12) return (S > cfg_.K) ? 1.0 : 0.0;
            double sqrt_tau = std::sqrt(tau);
            double d1 = (std::log(S / cfg_.K) + (cfg_.r + 0.5*cfg_.sigma_impl*cfg_.sigma_impl)*tau)
                       / (cfg_.sigma_impl * sqrt_tau);
            return 0.5 * std::erfc(-d1 * k_inv_sqrt2);
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
            double B     = C0 - delta * S;

            for (int k = 0; k < N_rebal; ++k) {
                double S_new = S * std::exp(drift_adj + vol_sqdt * nd(rng));
                double B_new = B * exp_r_dt;

                if (k == N_rebal - 1) {
                    samples[path] = B_new + delta * S_new - std::max(S_new - cfg_.K, 0.0);
                } else {
                    double delta_new = bs_delta(S_new, k + 1);
                    B_new -= (delta_new - delta) * S_new;
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

    // -------------------------------------------------------------------------
    // Сканирование по частоте ребалансировки
    // -------------------------------------------------------------------------
    // Для переданного набора значений N_rebal запускает симуляцию и вычисляет
    // статистики P&L. Дополнительно оценивает наклон зависимости log(std(P&L))
    // от log(Δt) с помощью линейной регрессии (теоретическое значение +0.5).
    //
    // Параметры:
    //   N_rebal_list – вектор чисел шагов для анализа
    //   N_paths      – число сценариев Монте-Карло на каждую симуляцию
    //   base_seed    – базовое зерно; для i-го элемента используется
    //                  base_seed + 17·i для обеспечения независимости
    //
    // Возвращает структуру ScanResult с временными шагами Δt, средними,
    // стандартными отклонениями и оценкой наклона.
    // -------------------------------------------------------------------------
    struct ScanResult {
        std::vector<int>    N_rebal_list;   // входной массив N
        std::vector<double> dt_list;        // Δt = T / N
        std::vector<double> mean_pnl;       // E[P&L] для каждого N
        std::vector<double> std_pnl;        // std(P&L) для каждого N
        double              slope_log_std_vs_log_dt; // наклон в log-log координатах
    };

    // -------------------------------------------------------------------------
    // Сканирование зависимости статистик P&L от частоты ребалансировки
    // -------------------------------------------------------------------------
    // Выполняет серию симуляций дельта-хеджирования для различных значений
    // числа шагов N_rebal из переданного списка N_rebal_list. Для каждого N
    // запускается simulate() с заданным числом сценариев N_paths, и
    // сохраняются среднее и стандартное отклонение P&L.
    //
    // Дополнительно вычисляется наклон линейной регрессии в логарифмических
    // координатах: log(std(P&L)) = α · log(Δt) + const. Согласно теории
    // (Bertsimas, Kogan, Lo 1998) при σ_real = σ_impl ожидается α ≈ 0.5.
    //
    // Параметры:
    //   N_rebal_list – вектор значений N (число ребалансировок) для анализа
    //   N_paths      – число сценариев Монте-Карло для каждой симуляции
    //   base_seed    – базовое зерно генератора случайных чисел. Для i-го
    //                  элемента списка используется зерно base_seed + 17*i,
    //                  что обеспечивает независимость выборок.
    //
    // Возвращает структуру ScanResult, содержащую:
    //   • исходный список N_rebal_list
    //   • соответствующие значения Δt = T / N
    //   • векторы средних и стандартных отклонений P&L для каждого N
    //   • оценку наклона slope_log_std_vs_log_dt
    // -------------------------------------------------------------------------
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

    static constexpr double k_inv_sqrt2 = 0.7071067811865475244; // 1/sqrt(2)
};
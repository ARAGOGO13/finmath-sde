// =============================================================================
// Модель Блэка–Шоулза для европейского колл-опциона на бездивидендный актив,
// следующий геометрическому броуновскому движению.
//
// Под риск-нейтральной мерой Q динамика цены актива:
//   dS_t = r S_t dt + σ S_t dW_t^Q
//
// Аналитическая формула для колла:
//   C = S0 Φ(d1) - K e^{-rT} Φ(d2)
//   d1 = [ln(S0/K) + (r + σ²/2) T] / (σ √T)
//   d2 = d1 - σ √T
//
// Также содержит методы Монте-Карло для оценки цены колла:
//   - под риск-нейтральной мерой Q (семплирование с дрейфом r),
//   - под физической мерой P с дрейфом μ (наивный оценщик).
// =============================================================================

#pragma once

#include <cmath>
#include <random>
#include <stdexcept>

// Результат оценки цены методом Монте-Карло
struct MCPriceResult {
    double price;      // оценка цены
    double std_error;  // стандартная ошибка
    double ci_lo;      // нижняя граница 95% доверительного интервала
    double ci_hi;      // верхняя граница 95% доверительного интервала
    int N;             // число сценариев
};

class BlackScholes {
public:
    // Параметры: безрисковая ставка, волатильность, срок, страйк, начальная цена.
    // Все должны быть положительными (кроме r, который может быть любым).
    BlackScholes(double r, double sigma, double T, double K, double S0)
        : r_(r), sigma_(sigma), T_(T), K_(K), S0_(S0) {
        if (T <= 0.0) throw std::invalid_argument("T must be positive");
        if (sigma <= 0.0) throw std::invalid_argument("sigma must be positive");
        if (S0 <= 0.0) throw std::invalid_argument("S0 must be positive");
        if (K <= 0.0) throw std::invalid_argument("K must be positive");
    }

    double r() const { return r_; }
    double sigma() const { return sigma_; }
    double T() const { return T_; }
    double K() const { return K_; }
    double S0() const { return S0_; }

    // Функция стандартного нормального распределения Φ(x)
    static double norm_cdf(double x) { return 0.5 * std::erfc(-x / k_sqrt2); }
    // Плотность стандартного нормального распределения φ(x)
    static double norm_pdf(double x) { return std::exp(-0.5 * x * x) / k_sqrt2pi; }

    // d1 = [ln(S0/K) + (r + σ²/2)T] / (σ√T)
    double d1() const { return (std::log(S0_ / K_) + (r_ + 0.5 * sigma_ * sigma_) * T_) / (sigma_ * std::sqrt(T_)); }
    // d2 = d1 - σ√T
    double d2() const { return d1() - sigma_ * std::sqrt(T_); }

    // Аналитическая цена колла (Блэк–Шоулз)
    double call_price() const { return S0_ * norm_cdf(d1()) - K_ * std::exp(-r_ * T_) * norm_cdf(d2()); }

    // Монте-Карло оценка цены колла под риск-нейтральной мерой Q.
    // Семплируется S_T с дрейфом r, затем дисконтированный выигрыш усредняется.
    MCPriceResult price_call_mc_rn(int N, unsigned seed = 42) const {
        std::mt19937_64 rng(seed);
        std::normal_distribution<double> nd(0.0, 1.0);
        const double drift_adj = (r_ - 0.5 * sigma_ * sigma_) * T_;
        const double vol_adj = sigma_ * std::sqrt(T_);
        const double discount = std::exp(-r_ * T_);
        double sum = 0.0, sum2 = 0.0;
        for (int i = 0; i < N; ++i) {
            double S_T = S0_ * std::exp(drift_adj + vol_adj * nd(rng));
            double dpay = discount * std::max(S_T - K_, 0.0);
            sum += dpay;
            sum2 += dpay * dpay;
        }
        double price = sum / N;
        double variance = sum2 / N - price * price;
        double std_err = std::sqrt(std::max(variance, 0.0) / N);
        return {price, std_err, price - 1.96 * std_err, price + 1.96 * std_err, N};
    }

    // Аналитическая цена колла под физической мерой P с дрейфом μ.
    // C0_P(μ) = e^{-rT} [ S0 e^{μT} Φ(d1(μ)) - K Φ(d2(μ)) ].
    // При μ = r совпадает с call_price().
    double call_price_physical(double mu) const {
        double sqT = std::sqrt(T_);
        double d1_mu = (std::log(S0_ / K_) + (mu + 0.5 * sigma_ * sigma_) * T_) / (sigma_ * sqT);
        double d2_mu = d1_mu - sigma_ * sqT;
        return std::exp(-r_ * T_) * (S0_ * std::exp(mu * T_) * norm_cdf(d1_mu) - K_ * norm_cdf(d2_mu));
    }

    // Монте-Карло оценка цены колла под физической мерой P (наивный подход).
    // Семплируется S_T с дрейфом μ, дисконтируется по r, весовая коррекция не применяется.
    MCPriceResult price_call_mc_physical(double mu, int N, unsigned seed) const {
        std::mt19937_64 rng(seed);
        std::normal_distribution<double> nd(0.0, 1.0);
        const double drift_adj = (mu - 0.5 * sigma_ * sigma_) * T_;
        const double vol_adj = sigma_ * std::sqrt(T_);
        const double discount = std::exp(-r_ * T_);
        double sum = 0.0, sum2 = 0.0;
        for (int i = 0; i < N; ++i) {
            double S_T = S0_ * std::exp(drift_adj + vol_adj * nd(rng));
            double dpay = discount * std::max(S_T - K_, 0.0);
            sum += dpay;
            sum2 += dpay * dpay;
        }
        double price = sum / N;
        double variance = sum2 / N - price * price;
        double std_err = std::sqrt(std::max(variance, 0.0) / N);
        return {price, std_err, price - 1.96 * std_err, price + 1.96 * std_err, N};
    }

private:
    double r_, sigma_, T_, K_, S0_;
    static constexpr double k_sqrt2 = 1.4142135623730950488;   // √2
    static constexpr double k_sqrt2pi = 2.5066282746310002416; // √(2π)
};
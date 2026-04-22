#pragma once

#include <cmath>
#include <random>
#include <stdexcept>

// =============================================================================
// Модель Блэка–Шоулза для европейского опциона колл/пут на актив, следующий GBM
//
// Под риск-нейтральной мерой Q динамика базового актива описывается СДУ:
//   dS_t = r S_t dt + σ S_t dW_t^Q
//
// где:
//   r     – безрисковая процентная ставка (годовая, непрерывное начисление)
//   σ     – волатильность актива (годовая)
//   W_t^Q – стандартное броуновское движение под мерой Q
//
// Аналитические формулы для европейского колла:
//   C = S0 Φ(d1) - K e^{-rT} Φ(d2)
//   P = K e^{-rT} Φ(-d2) - S0 Φ(-d1)   (из паритета пут-колл)
//   d1 = [ln(S0/K) + (r + σ²/2) T] / (σ √T)
//   d2 = d1 - σ √T
//
// Греки (чувствительности) для колла:
//   Δ = ∂C/∂S = Φ(d1)
//   Γ = ∂²C/∂S² = φ(d1) / (S0 σ √T)
//   ν = ∂C/∂σ = S0 φ(d1) √T
//
// Класс также предоставляет методы Монте-Карло для оценки цены колла:
//   • price_call_mc_rn      – семплирование терминального значения под мерой Q
//   • price_call_mc_girsanov – семплирование под физической мерой P с
//                              перевзвешиванием производной Радона–Никодима
//
// Все вероятностные вычисления используют высокоточные аппроксимации
// стандартного нормального распределения (CDF и PDF).
// =============================================================================

// Структура для возврата результатов Монте-Карло оценки цены колла
struct MCPriceResult {
    double price;      // оценка цены (дисконтированный средний payoff)
    double std_error;  // стандартная ошибка оценки
    double ci_lo;      // нижняя граница 95% доверительного интервала
    double ci_hi;      // верхняя граница 95% доверительного интервала
    int N;             // число использованных сценариев
};

class BlackScholes {
public:
    // Принимает параметры модели: безрисковую ставку r, волатильность sigma,
    // время до экспирации T (в годах), страйк K и текущую цену актива S0.
    // Проверяет, что T, sigma, S0, K положительны.
    // -------------------------------------------------------------------------
    BlackScholes(double r, double sigma, double T, double K, double S0) : r_(r), sigma_(sigma), T_(T), K_(K), S0_(S0) {
        if (T <= 0.0)
            throw std::invalid_argument("T must be positive");
        if (sigma <= 0.0)
            throw std::invalid_argument("sigma must be positive");
        if (S0 <= 0.0)
            throw std::invalid_argument("S0 must be positive");
        if (K <= 0.0)
            throw std::invalid_argument("K must be positive");
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

    // d1 = [ln(S0/K) + (r + sigma^2/2)*T] / (sigma*sqrt(T))
    double d1() const { return (std::log(S0_ / K_) + (r_ + 0.5 * sigma_ * sigma_) * T_) / (sigma_ * std::sqrt(T_)); }

    // d2 = d1 - sigma*sqrt(T)
    double d2() const { return d1() - sigma_ * std::sqrt(T_); }

    // Аналитическая цена колла: C = S0*Phi(d1) - K*exp(-rT)*Phi(d2)
    double call_price() const { return S0_ * norm_cdf(d1()) - K_ * std::exp(-r_ * T_) * norm_cdf(d2()); }
    // Аналитическая цена пута (put-call parity): P = C - S0 + K*exp(-rT)
    double put_price() const { return call_price() - S0_ + K_ * std::exp(-r_ * T_); }
    // Дельта колла: Phi(d1)
    double call_delta() const { return norm_cdf(d1()); }
    // Гамма колла: phi(d1) / (S0 * sigma * sqrt(T))
    double call_gamma() const { return norm_pdf(d1()) / (S0_ * sigma_ * std::sqrt(T_)); }
    // Вега колла: S0 * phi(d1) * sqrt(T)
    double call_vega() const { return S0_ * norm_pdf(d1()) * std::sqrt(T_); }

    // -------------------------------------------------------------------------
    // Монте-Карло оценка цены колла под физической мерой P с перевзвешиванием
    // -------------------------------------------------------------------------
    // Симулирует S_T под физическим дрейфом μ, затем корректирует payoff
    // с помощью производной Радона–Никодима dQ/dP.
    // Параметры:
    //   mu   – физический дрейф актива
    //   N    – число сценариев
    //   seed – зерно генератора (должно отличаться от используемого в Q-MC)
    // Возвращает структуру MCPriceResult.
    // Рыночная цена риска λ = (μ - r)/σ неявно участвует в весовой функции.
    // -------------------------------------------------------------------------
    MCPriceResult price_call_mc_girsanov(double mu, int N, unsigned seed) const {
        std::mt19937_64 rng(seed);
        std::normal_distribution<double> nd(0.0, 1.0);

        const double lambda   = (mu - r_) / sigma_;
        const double sqrtT    = std::sqrt(T_);
        const double discount = std::exp(-r_ * T_);

        double sum = 0.0, sum2 = 0.0;
        for (int i = 0; i < N; ++i) {
            double Z   = nd(rng);
            double W_T = sqrtT * Z;
            // Цена актива под мерой P (дрейф mu, а не r)
            double S_T = S0_ * std::exp((mu - 0.5 * sigma_ * sigma_) * T_ + sigma_ * W_T);
            // Производная Радона-Никодима dQ/dP|_{F_T}
            double weight          = std::exp(-0.5 * lambda * lambda * T_ - lambda * W_T);
            double payoff_weighted = std::max(S_T - K_, 0.0) * weight;
            double dpay            = discount * payoff_weighted;
            sum  += dpay;
            sum2 += dpay * dpay;
        }

        double price    = sum / N;
        double variance = sum2 / N - price * price;
        double std_err  = std::sqrt(std::max(variance, 0.0) / N);

        return MCPriceResult{price, std_err, price - k_z95 * std_err, price + k_z95 * std_err, N};
    }

    // -------------------------------------------------------------------------
    // Монте-Карло оценка цены колла под риск-нейтральной мерой Q
    // -------------------------------------------------------------------------
    // Использует точное семплирование терминального значения S_T.
    // Параметры:
    //   N    – число сценариев
    //   seed – зерно генератора псевдослучайных чисел (по умолчанию 42)
    // Возвращает структуру MCPriceResult с ценой, стандартной ошибкой и ДИ.
    // -------------------------------------------------------------------------
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

        return MCPriceResult{price, std_err, price - k_z95 * std_err, price + k_z95 * std_err, N};
    }

private:
    double r_, sigma_, T_, K_, S0_;

    static constexpr double k_sqrt2 = 1.4142135623730950488; // sqrt(2)
    static constexpr double k_sqrt2pi = 2.5066282746310002416; // sqrt(2*pi)
    static constexpr double k_z95 = 1.96; // квантиль 95% ДИ
};

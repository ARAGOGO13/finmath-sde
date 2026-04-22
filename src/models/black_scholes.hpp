#pragma once
// =============================================================================
// BlackScholes — аналитическая модель ценообразования Блэка-Шоулза для
// европейского опциона колл/пут на базовый актив, следующий GBM.
//
// Подлежащий актив под риск-нейтральной мерой Q (Brigo-Mercurio, прил. C):
//   dS = r * S * dt + sigma * S * dW^Q
//
// Аналитические формулы колла и грек:
//   d1  = [ln(S0/K) + (r + sigma^2/2) * T] / (sigma * sqrt(T))
//   d2  = d1 - sigma * sqrt(T)
//   C   = S0 * Phi(d1) - K * exp(-r*T) * Phi(d2)
//   P   = K * exp(-r*T) * Phi(-d2) - S0 * Phi(-d1)   (put-call parity)
//   Delta = Phi(d1),   Gamma = phi(d1)/(S0*sigma*sqrt(T)),   Vega = S0*phi(d1)*sqrt(T)
//
// MC-оценщик использует прямое семплирование терминального значения
// (exact sampling, нулевое смещение по сравнению с дискретизацией):
//   S_T = S0 * exp((r - sigma^2/2)*T + sigma*sqrt(T)*Z),  Z ~ N(0,1)
//   цена = exp(-r*T) * среднее(max(S_T - K, 0))
// =============================================================================

#include <cmath>
#include <random>
#include <stdexcept>

// Результат MC-оценки цены колла.
struct MCPriceResult {
    double price; // дисконтированный средний payoff
    double std_error; // стандартная ошибка MC
    double ci_lo; // нижняя граница 95% ДИ
    double ci_hi; // верхняя граница 95% ДИ
    int N; // число сценариев
};

class BlackScholes {
public:
    // r     — безрисковая ставка
    // sigma — волатильность базового актива
    // T     — время до экспирации (в годах)
    // K     — страйк
    // S0    — текущая цена базового актива
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

    // Стандартная нормальная CDF через erfc:  Phi(x) = erfc(-x/sqrt(2)) / 2
    static double norm_cdf(double x) { return 0.5 * std::erfc(-x / k_sqrt2); }

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

    // MC-оценка цены колла под физической мерой P с перевзвешиванием Радона-Никодима.
    // lambda = (mu - r) / sigma — рыночная цена риска.
    // При mu == r: lambda=0, weight=1, результат совпадает с price_call_mc_rn.
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

    // MC-оценка цены колла под мерой Q (exact sampling терминального значения).
    // seed — явный параметр для воспроизводимости результатов.
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

    // Стандартная нормальная плотность: phi(x) = exp(-x^2/2) / sqrt(2*pi)
    static double norm_pdf(double x) { return std::exp(-0.5 * x * x) / k_sqrt2pi; }

    static constexpr double k_sqrt2 = 1.4142135623730950488; // sqrt(2)
    static constexpr double k_sqrt2pi = 2.5066282746310002416; // sqrt(2*pi)
    static constexpr double k_z95 = 1.96; // квантиль 95% ДИ
};

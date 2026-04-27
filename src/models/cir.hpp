// =============================================================================
// Модель Кокса–Ингерсолла–Росса (CIR).
//
// Стохастическое дифференциальное уравнение:
//   dX_t = θ (μ - X_t) dt + σ √X_t dW_t
//
// Процесс с возвратом к среднему μ, волатильность пропорциональна √X_t.
// При выполнении условия Феллера (2θμ ≥ σ²) процесс остаётся строго положительным.
// =============================================================================

#pragma once
#include <cmath>
#include <stdexcept>
#include "../sde_base.hpp"

class CoxIngersollRoss : public SDE {
public:
    CoxIngersollRoss(double theta = 2.0, double mu = 0.05, double sigma = 0.1) :
        theta_(theta), mu_(mu), sigma_(sigma) {}

    // Коэффициент сноса: θ (μ - x)
    double drift(double x, double) const override { return theta_ * (mu_ - x); }

    // Коэффициент диффузии: σ √max(x,0) (защита от отрицательных значений)
    double diffusion(double x, double) const override { return sigma_ * std::sqrt(std::max(x, 0.0)); }

    // Производная коэффициента диффузии: σ / (2√x) для x > 0
    double diffusion_derivative(double x, double) const override {
        if (x <= 1e-10)
            return 0.0;
        return sigma_ / (2.0 * std::sqrt(x));
    }

    double theta() const { return theta_; }
    double mu() const { return mu_; }
    double sigma() const { return sigma_; }

    // Проверка условия Феллера (2θμ ≥ σ²)
    bool feller_condition() const { return 2.0 * theta_ * mu_ >= sigma_ * sigma_; }

private:
    double theta_, mu_, sigma_;
};

// Аналитическая дисперсия процесса CIR в момент времени t
inline double cir_exact_variance(double x0, double t, double theta, double mu, double sigma) {
    double e  = std::exp(-theta * t);
    double s2 = sigma * sigma;
    return x0 * (s2 / theta) * e * (1.0 - e)
         + (mu * s2) / (2.0 * theta) * (1.0 - e) * (1.0 - e);
}
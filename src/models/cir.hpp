#pragma once
#include <cmath>
#include <stdexcept>
#include "../sde_base.hpp"

// =============================================================================
// Модель Кокса–Ингерсолла–Росса (CIR)
//
// Стохастическое дифференциальное уравнение:
//   dX_t = θ (μ - X_t) dt + σ √X_t dW_t
//
// где:
//   θ (theta) – скорость возврата к среднему
//   μ (mu)    – долгосрочное среднее значение
//   σ (sigma) – волатильность (множитель при √X_t)
//
// Процесс X_t остаётся неотрицательным при выполнении условия Феллера:
//   2 θ μ ≥ σ².
// =============================================================================

class CoxIngersollRoss : public SDE {
public:
    CoxIngersollRoss(double theta = 2.0, double mu = 0.05, double sigma = 0.1) :
        theta_(theta), mu_(mu), sigma_(sigma) {}

    double drift(double x, double) const override { return theta_ * (mu_ - x); }

    double diffusion(double x, double) const override { return sigma_ * std::sqrt(std::max(x, 0.0)); }

    double diffusion_derivative(double x, double) const override {
        if (x <= 1e-10)
            return 0.0;
        return sigma_ / (2.0 * std::sqrt(x));
    }

    double theta() const { return theta_; }
    double mu() const { return mu_; }
    double sigma() const { return sigma_; }

    bool feller_condition() const { return 2.0 * theta_ * mu_ >= sigma_ * sigma_; }

private:
    double theta_, mu_, sigma_;
};

// =============================================================================
// Аналитическая дисперсия процесса CIR
// =============================================================================
// Вычисляет теоретическую дисперсию Var[X_t] для процесса CIR с параметрами
// theta, mu, sigma и начальным значением x0 в момент времени t.
// Формула: Var[X_t] = x0 * (σ²/θ) * e^{-θt} * (1 - e^{-θt})
//                   + (μ σ²)/(2θ) * (1 - e^{-θt})²
// =============================================================================
inline double cir_exact_variance(double x0, double t, double theta, double mu, double sigma) {
    double e  = std::exp(-theta * t);
    double s2 = sigma * sigma;
    return x0 * (s2 / theta) * e * (1.0 - e)
         + (mu * s2) / (2.0 * theta) * (1.0 - e) * (1.0 - e);
}

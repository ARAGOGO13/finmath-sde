#pragma once
#include "../sde_base.hpp"
#include <cmath>

// =============================================================================
// OrnsteinUhlenbeck / Vasicek — модель mean-reversion (C.4)
//
//   dX_t = theta * (mu - X_t) * dt + sigma * dW_t
//
// Это модель Vasicek (1978) для краткосрочной процентной ставки.
// theta — скорость возврата к среднему (mean reversion speed)
// mu    — долгосрочное среднее (long-run mean)
// sigma — волатильность (константа!)
//
// Диффузия: sigma(x) = sigma = const
// Производная: sigma'(x) = 0
// => схемы Эйлера и Мильштейна СОВПАДАЮТ (см. раздел C.3 PDF)
//
// Точное распределение: X_t | X_0 ~ N(mean(t), var(t)) где:
//   mean(t) = X_0 * exp(-theta*t) + mu * (1 - exp(-theta*t))
//   var(t)  = sigma^2 / (2*theta) * (1 - exp(-2*theta*t))
// =============================================================================

class OrnsteinUhlenbeck : public SDE {
public:
    OrnsteinUhlenbeck(double theta = 5.0, double mu = 0.0, double sigma = 1.0)
        : theta_(theta), mu_(mu), sigma_(sigma) {}

    double drift(double x, double /*t*/) const override {
        return theta_ * (mu_ - x);
    }

    double diffusion(double /*x*/, double /*t*/) const override {
        return sigma_;
    }

    // sigma'(x) = 0 — диффузия константна
    double diffusion_derivative(double /*x*/, double /*t*/) const override {
        return 0.0;
    }

    // Точное среднее E[X_t | X_0 = x0]
    double exact_mean(double x0, double t) const {
        return x0 * std::exp(-theta_ * t) + mu_ * (1.0 - std::exp(-theta_ * t));
    }

    // Точная дисперсия Var[X_t | X_0 = x0]
    double exact_var(double t) const {
        return (sigma_ * sigma_) / (2.0 * theta_) * (1.0 - std::exp(-2.0 * theta_ * t));
    }

    double theta() const { return theta_; }
    double mu()    const { return mu_; }
    double sigma() const { return sigma_; }

private:
    double theta_, mu_, sigma_;
};
#pragma once
#include "../sde_base.hpp"

// =============================================================================
// GeometricBrownianMotion — модель Блэка–Шоулса (C.4)
//
//   dX_t = mu * X_t * dt + sigma * X_t * dW_t
//
// Точное решение (по формуле Ито):
//   X_t = X_0 * exp((mu - 0.5*sigma^2)*t + sigma*W_t)
//
// Диффузия: sigma(x) = sigma * x
// Производная: sigma'(x) = sigma  => поправка Мильштейна ненулевая
// =============================================================================

class GeometricBrownianMotion : public SDE {
public:
    GeometricBrownianMotion(double mu = 0.05, double sigma = 0.20)
        : mu_(mu), sigma_(sigma) {}

    double drift(double x, double /*t*/) const override {
        return mu_ * x;
    }

    double diffusion(double x, double /*t*/) const override {
        return sigma_ * x;
    }

    // sigma'(x) = sigma — константа, не зависит от x
    double diffusion_derivative(double /*x*/, double /*t*/) const override {
        return sigma_;
    }

    // Точное решение для верификации схем
    double exact(double x0, double t, double W_t) const {
        return x0 * std::exp((mu_ - 0.5 * sigma_ * sigma_) * t + sigma_ * W_t);
    }

    double mu()    const { return mu_; }
    double sigma() const { return sigma_; }

private:
    double mu_, sigma_;
};
#pragma once
#include "../sde_base.hpp"

class GeometricBrownianMotion : public SDE {
public:
    GeometricBrownianMotion(double mu = 0.05, double sigma = 0.20) : mu_(mu), sigma_(sigma) {}

    double drift(double x, double) const override { return mu_ * x; }

    double diffusion(double x, double) const override { return sigma_ * x; }

    double diffusion_derivative(double, double) const override { return sigma_; }

    double exact(double x0, double t, double W_t) const {
        return x0 * std::exp((mu_ - 0.5 * sigma_ * sigma_) * t + sigma_ * W_t);
    }

    double mu() const { return mu_; }
    double sigma() const { return sigma_; }

private:
    double mu_, sigma_;
};

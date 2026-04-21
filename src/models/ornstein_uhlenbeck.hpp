#pragma once
#include <cmath>
#include "../sde_base.hpp"

class OrnsteinUhlenbeck : public SDE {
public:
    OrnsteinUhlenbeck(double theta = 5.0, double mu = 0.0, double sigma = 1.0) :
        theta_(theta), mu_(mu), sigma_(sigma) {}

    double drift(double x, double) const override { return theta_ * (mu_ - x); }

    double diffusion(double, double) const override { return sigma_; }

    double diffusion_derivative(double, double) const override { return 0.0; }

    double exact_mean(double x0, double t) const {
        return x0 * std::exp(-theta_ * t) + mu_ * (1.0 - std::exp(-theta_ * t));
    }

    double exact_var(double t) const {
        return (sigma_ * sigma_) / (2.0 * theta_) * (1.0 - std::exp(-2.0 * theta_ * t));
    }

    double theta() const { return theta_; }
    double mu() const { return mu_; }
    double sigma() const { return sigma_; }

private:
    double theta_, mu_, sigma_;
};

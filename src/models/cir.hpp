#pragma once
#include <cmath>
#include <stdexcept>
#include "../sde_base.hpp"

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

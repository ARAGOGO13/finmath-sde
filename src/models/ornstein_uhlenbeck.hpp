// =============================================================================
// Модель Орнштейна–Уленбека (процесс Vasicek).
//
// Стохастическое дифференциальное уравнение:
//   dX_t = θ (μ - X_t) dt + σ dW_t
//
// Процесс обладает свойством возврата к среднему μ со скоростью θ.
// Диффузия постоянна, поэтому схема Эйлера и Мильштейна совпадают.
// Используется для моделирования процентных ставок и других средневозвратных величин.
// =============================================================================

#pragma once
#include <cmath>
#include "../sde_base.hpp"

class OrnsteinUhlenbeck : public SDE {
public:
    OrnsteinUhlenbeck(double theta = 5.0, double mu = 0.0, double sigma = 1.0) :
        theta_(theta), mu_(mu), sigma_(sigma) {}

    // Коэффициент сноса: θ (μ - x)
    double drift(double x, double) const override { return theta_ * (mu_ - x); }

    // Коэффициент диффузии: σ (постоянный)
    double diffusion(double, double) const override { return sigma_; }

    // Производная коэффициента диффузии: 0 (постоянная волатильность)
    double diffusion_derivative(double, double) const override { return 0.0; }

    // Точное математическое ожидание в момент t
    double exact_mean(double x0, double t) const {
        return x0 * std::exp(-theta_ * t) + mu_ * (1.0 - std::exp(-theta_ * t));
    }

    // Точная дисперсия в момент t
    double exact_var(double t) const {
        return (sigma_ * sigma_) / (2.0 * theta_) * (1.0 - std::exp(-2.0 * theta_ * t));
    }

    double theta() const { return theta_; }
    double mu() const { return mu_; }
    double sigma() const { return sigma_; }

private:
    double theta_, mu_, sigma_;
};
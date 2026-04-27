// =============================================================================
// Модель геометрического броуновского движения (GBM).
//
// Стохастическое дифференциальное уравнение:
//   dS_t = μ S_t dt + σ S_t dW_t
//
// Аналитическое решение: S_t = S_0 exp( (μ - ½σ²)t + σ W_t ).
// Используется как базовый процесс для модели Блэка–Шоулза.
// =============================================================================

#pragma once
#include "../sde_base.hpp"

class GeometricBrownianMotion : public SDE {
public:
    GeometricBrownianMotion(double mu = 0.05, double sigma = 0.20) : mu_(mu), sigma_(sigma) {}

    // Коэффициент сноса: μ x
    double drift(double x, double) const override { return mu_ * x; }

    // Коэффициент диффузии: σ x
    double diffusion(double x, double) const override { return sigma_ * x; }

    // Производная коэффициента диффузии: σ (константа)
    double diffusion_derivative(double, double) const override { return sigma_; }

    // Точное решение в момент t при известном значении броуновского движения W_t
    double exact(double x0, double t, double W_t) const {
        return x0 * std::exp((mu_ - 0.5 * sigma_ * sigma_) * t + sigma_ * W_t);
    }

    double mu() const { return mu_; }
    double sigma() const { return sigma_; }

private:
    double mu_, sigma_;
};
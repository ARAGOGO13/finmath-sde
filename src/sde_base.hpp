#pragma once

// =============================================================================
// sde_base.hpp
// Базовый интерфейс для всех СДУ вида:
//   dX_t = f(X_t, t) dt + sigma(X_t, t) dW_t
//
// diffusion_derivative нужна для схемы Мильштейна:
//   поправка = 0.5 * sigma * sigma' * (dW^2 - dt)
// Для моделей с константной диффузией (Vasicek) sigma' = 0,
// и схема Мильштейна совпадает с Эйлером — что и доказывается численно.
// =============================================================================

struct SDE {
    virtual double drift(double x, double t)      const = 0;
    virtual double diffusion(double x, double t)  const = 0;

    // Производная диффузии по x: d(sigma)/dx
    // По умолчанию 0 — подходит для моделей с детерминированной диффузией
    virtual double diffusion_derivative(double x, double t) const { return 0.0; }

    virtual ~SDE() = default;
};
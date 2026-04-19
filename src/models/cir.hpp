#pragma once
#include "../sde_base.hpp"
#include <cmath>
#include <stdexcept>

// =============================================================================
// CoxIngersollRoss — квадратно-корневой процесс (C.4)
//
//   dX_t = theta * (mu - X_t) * dt + sigma * sqrt(X_t) * dW_t
//
// Модель CIR (1985) для краткосрочной процентной ставки.
// В отличие от Vasicek, гарантирует X_t >= 0 при условии Фeller:
//   2 * theta * mu >= sigma^2
//
// Диффузия: sigma(x) = sigma * sqrt(x)
// Производная: sigma'(x) = sigma / (2 * sqrt(x))  => поправка Мильштейна ненулевая
//
// Это НЕЛИНЕЙНОЕ СДУ — хороший пример для демонстрации разницы Эйлер/Мильштейн.
// =============================================================================

class CoxIngersollRoss : public SDE {
public:
    CoxIngersollRoss(double theta = 2.0, double mu = 0.05, double sigma = 0.1)
        : theta_(theta), mu_(mu), sigma_(sigma)
    {
        // Предупреждение если условие Феллера не выполнено
        if (2.0 * theta_ * mu_ < sigma_ * sigma_) {
            // Процесс может касаться нуля — sqrt(x) станет проблемой
            // В числовой схеме будем использовать max(X, 0)
        }
    }

    double drift(double x, double /*t*/) const override {
        return theta_ * (mu_ - x);
    }

    double diffusion(double x, double /*t*/) const override {
        return sigma_ * std::sqrt(std::max(x, 0.0));
    }

    // sigma'(x) = sigma / (2 * sqrt(x))
    double diffusion_derivative(double x, double /*t*/) const override {
        if (x <= 1e-10) return 0.0;  // защита от деления на ноль
        return sigma_ / (2.0 * std::sqrt(x));
    }

    double theta() const { return theta_; }
    double mu()    const { return mu_; }
    double sigma() const { return sigma_; }

    // Условие Феллера: если выполнено, процесс строго положителен
    bool feller_condition() const {
        return 2.0 * theta_ * mu_ >= sigma_ * sigma_;
    }

private:
    double theta_, mu_, sigma_;
};
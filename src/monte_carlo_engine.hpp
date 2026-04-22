#pragma once
#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <vector>
#include "sde_base.hpp"

// =============================================================================
// Универсальный движок Монте-Карло для симуляции траекторий стохастических ДУ
//
// Класс MonteCarloEngine предоставляет методы для генерации траекторий
// произвольного скалярного процесса, описываемого интерфейсом SDE.
// Поддерживаются две численные схемы дискретизации:
//   • Euler   – схема Эйлера–Маруямы (порядок сильной сходимости 0.5)
//   • Milshtein – схема Мильштейна (порядок сильной сходимости 1.0)
//
// Схема Эйлера:
//   X_{k+1} = X_k + f(X_k, t_k) Δt + g(X_k, t_k) ΔW_k
//
// Схема Мильштейна:
//   X_{k+1} = X_k + f(X_k, t_k) Δt + g(X_k, t_k) ΔW_k
//             + ½ g(X_k, t_k) g'(X_k, t_k) (ΔW_k² - Δt)
//
// где Δt = T / N_steps, ΔW_k ~ N(0, Δt), g' – производная коэффициента диффузии
// по переменной состояния.
//
// Движок позволяет задать число сценариев, число временных шагов, схему,
// а после симуляции предоставляет оценки математического ожидания и дисперсии
// процесса в конечный момент времени T.
// =============================================================================

enum class Scheme { Euler, Milshtein };

class MonteCarloEngine {
public:
    // -------------------------------------------------------------------------
    // Параметры:
    //   num_scenarios – число независимых траекторий Монте-Карло
    //   num_steps     – число шагов дискретизации на интервале [0, T]
    //   T             – длина временного интервала
    //   scheme        – выбор численной схемы (Euler по умолчанию)
    // -------------------------------------------------------------------------
    MonteCarloEngine(size_t num_scenarios, int num_steps, double T, Scheme scheme = Scheme::Euler) :
        num_scenarios_(num_scenarios), num_steps_(num_steps), dt_(T / num_steps), scheme_(scheme),
        gen_(std::random_device{}()), normal_(0.0, 1.0) {}

    // -------------------------------------------------------------------------
    // Симуляция траекторий
    // -------------------------------------------------------------------------
    // Параметры:
    //   model – ссылка на объект, реализующий интерфейс SDE
    //   x0    – начальное значение процесса
    //   num_paths_to_save – число траекторий, сохраняемых для последующего
    //                        анализа или визуализации (по умолчанию 50)
    //
    // Возвращает вектор сохранённых траекторий размером (num_paths_to_save × (num_steps+1)).
    // Попутно вычисляет и сохраняет во внутренние поля среднее и дисперсию
    // в конечный момент времени по всем num_scenarios траекториям.
    // -------------------------------------------------------------------------
    std::vector<std::vector<double>> simulate(const SDE &model, double x0, size_t num_paths_to_save = 50) {
        std::vector<std::vector<double>> paths(num_paths_to_save, std::vector<double>(num_steps_ + 1, 0.0));

        double sum_XT = 0.0;
        double sum_XT2 = 0.0;

        for (size_t j = 0; j < num_scenarios_; ++j) {
            double X = x0;
            if (j < num_paths_to_save)
                paths[j][0] = X;

            for (int k = 0; k < num_steps_; ++k) {
                double t = k * dt_;
                double dW = normal_(gen_) * std::sqrt(dt_);
                double f = model.drift(X, t);
                double s = model.diffusion(X, t);

                X += f * dt_ + s * dW;

                if (scheme_ == Scheme::Milshtein) {
                    double s_prime = model.diffusion_derivative(X, t);
                    X += 0.5 * s * s_prime * (dW * dW - dt_);
                }

                if (j < num_paths_to_save)
                    paths[j][k + 1] = X;
            }
            sum_XT += X;
            sum_XT2 += X * X;
        }

        mean_at_T_ = sum_XT / static_cast<double>(num_scenarios_);
        double mean2 = sum_XT2 / static_cast<double>(num_scenarios_);
        variance_at_T_ = mean2 - mean_at_T_ * mean_at_T_;

        return paths;
    }

    double get_mean_at_T() const;     // Монте-Карло оценка E[X_T]
    double get_variance_at_T() const; // Монте-Карло оценка Var[X_T]
    double get_std_at_T() const;      // Квадратный корень из дисперсии

    Scheme get_scheme() const { return scheme_; }

private:
    size_t num_scenarios_;   // число траекторий
    int num_steps_;          // число шагов дискретизации
    double dt_;              // величина шага Δt = T / num_steps
    Scheme scheme_;          // выбранная численная схема

    std::mt19937 gen_;                           // генератор псевдослучайных чисел
    std::normal_distribution<double> normal_;    // стандартное нормальное распределение N(0,1)

    mutable double mean_at_T_ = 0.0;      // оценка матожидания X_T
    mutable double variance_at_T_ = 0.0;  // оценка дисперсии X_T
};


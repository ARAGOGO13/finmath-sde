#pragma once
#include <cmath>
#include <matplot/matplot.h>
#include <string>
#include <vector>
#include "girsanov_experiment.hpp"

using namespace matplot;

// =============================================================================
// Визуализация результатов эксперимента по проверке теоремы Гирсанова
//
// Файл содержит функции для построения графиков, иллюстрирующих численную
// проверку теоремы Гирсанова и анализ «стоимости смены меры». Используется
// библиотека Matplot++ (matplot). Результаты сохраняются в формате PDF.
//
// Две основные функции:
//   • plot_girsanov_prices   – сравнение трёх оценщиков цены колла
//   • plot_girsanov_variance – отношение дисперсий оценок P-MC и Q-MC
// =============================================================================


// -----------------------------------------------------------------------------
// График сравнения трёх оценщиков цены колла
// -----------------------------------------------------------------------------
// Строит зависимость цены колла от рыночной цены риска λ = (μ - r)/σ.
// На одном полотне отображаются:
//   • аналитическая цена Блэка–Шоулза (чёрная сплошная линия) – эталон,
//     не зависящий от μ;
//   • оценка Монте-Карло под риск-нейтральной мерой Q (красный пунктир);
//   • оценка Монте-Карло под физической мерой P с перевзвешиванием
//     (синий пунктир).
//
// Для MC-оценщиков также показаны 95% доверительные интервалы
// (±1.96 × стандартная ошибка) в виде вертикальных отрезков.
//
// Параметры:
//   pts      – вектор результатов GirsanovExperiment::Point, содержащий
//              значения λ и все три оценки цены с их ошибками
//   filename – имя выходного PDF-файла
// -----------------------------------------------------------------------------
inline void plot_girsanov_prices(
        const std::vector<GirsanovExperiment::Point>& pts,
        const std::string& filename) {

    std::vector<double> lambdas, bs_vals, rn_vals, p_vals, rn_err, p_err;
    lambdas.reserve(pts.size());
    for (const auto& pt : pts) {
        lambdas.push_back(pt.lambda);
        bs_vals.push_back(pt.bs_price);
        rn_vals.push_back(pt.rn_price);
        p_vals .push_back(pt.p_price);
        rn_err .push_back(1.96 * pt.rn_std_err);
        p_err  .push_back(1.96 * pt.p_std_err);
    }

    auto f = figure(true); f->size(1200, 680);
    auto ax = gca();
    hold(ax, true);

    { auto l = plot(ax, lambdas, bs_vals, "-");  l->line_width(3.0); l->color("black"); }
    { auto l = plot(ax, lambdas, rn_vals, "--"); l->line_width(2.2); l->color("red");   }
    { auto l = plot(ax, lambdas, p_vals,  "--"); l->line_width(2.2); l->color("blue");  }

    {
        auto h = errorbar(ax, lambdas, rn_vals, rn_err,
                          error_bar::type::vertical, "--");
        h->color("red"); h->line_width(1.2); h->cap_size(4.f);
    }
    {
        auto h = errorbar(ax, lambdas, p_vals, p_err,
                          error_bar::type::vertical, "--");
        h->color("blue"); h->line_width(1.2); h->cap_size(4.f);
    }

    ax->title("Girsanov verification: three estimators of the same call price");
    ax->xlabel("lambda = (mu - r) / sigma   (market price of risk)");
    ax->ylabel("Call price");
    legend(ax, {"BS analytic (black)", "Q-MC risk-neutral (red)", "P-MC Girsanov (blue)"});
    grid(ax, true);
    save(filename);
    cla();
}

// -----------------------------------------------------------------------------
// Вспомогательная функция: второй момент выплаты колла (E[(S_T - K)^+²])
// -----------------------------------------------------------------------------
// Вычисляет теоретическое значение E[(S_T - K)^+²] для логнормального
// распределения S_T с заданным дрейфом mu_eff. Используется для получения
// точного выражения дисперсии оценки Монте-Карло при заданном λ.
// -----------------------------------------------------------------------------
inline double second_moment_call_(double S0, double K, double mu_eff,
                                  double sigma, double T) {
    const double k_inv_sqrt2 = 0.7071067811865475244;
    auto Phi = [&](double x) { return 0.5 * std::erfc(-x * k_inv_sqrt2); };
    double sqT = std::sqrt(T);
    double d1  = (std::log(S0/K) + (mu_eff + 0.5*sigma*sigma)*T) / (sigma*sqT);
    double d2  = d1 - sigma*sqT;
    double d1p = d1 + sigma*sqT;
    return S0*S0 * std::exp((2*mu_eff + sigma*sigma)*T) * Phi(d1p)
         - 2*K*S0 * std::exp(mu_eff*T) * Phi(d1)
         + K*K * Phi(d2);
}

// -----------------------------------------------------------------------------
// График отношения дисперсий оценок P-MC и Q-MC
// -----------------------------------------------------------------------------
// Строит зависимость Var(P-MC) / Var(Q-MC) от λ в логарифмическом масштабе.
// Показывает:
//   • теоретическую кривую для колл-опциона (чёрная сплошная линия),
//     учитывающую корреляцию между выплатой и весом Радона–Никодима;
//   • Монте-Карло оценку отношения дисперсий (синие кружки с пунктирной
//     линией);
//   • для справки – экспоненту exp(λ² T) (серый пунктир), которая является
//     отношением дисперсий для выплат, некоррелированных с весом.
//
// График наглядно демонстрирует, что при λ > 0 вес гасит высокие выплаты,
// уменьшая дисперсию (ratio < 1), а при λ < 0 – усиливает, вызывая рост
// дисперсии (ratio ≫ 1).
//
// Параметры:
//   pts      – результаты эксперимента (используются λ и variance_ratio)
//   filename – имя выходного PDF-файла
//   T, r, sigma, K, S0 – параметры модели, необходимые для теоретической
//                        кривой
// -----------------------------------------------------------------------------
inline void plot_girsanov_variance(
        const std::vector<GirsanovExperiment::Point>& pts,
        const std::string& filename,
        double T = 1.0,
        double r = 0.05, double sigma = 0.20,
        double K = 100.0, double S0 = 100.0) {

    const double k_inv_sqrt2 = 0.7071067811865475244;
    auto Phi = [&](double x) { return 0.5 * std::erfc(-x * k_inv_sqrt2); };
    double sqT = std::sqrt(T);
    double d1  = (std::log(S0/K) + (r + 0.5*sigma*sigma)*T) / (sigma*sqT);
    double d2  = d1 - sigma*sqT;
    double C   = S0 * Phi(d1) - K * std::exp(-r*T) * Phi(d2);

    double varQ = std::exp(-2*r*T) * second_moment_call_(S0, K, r, sigma, T) - C*C;

    std::vector<double> lambdas, ratios, theory_call, theory_rn;
    lambdas.reserve(pts.size());
    for (const auto& pt : pts) {
        double lam    = pt.lambda;
        double mu     = r + lam*sigma;
        double mu_eff = 2*r - mu;
        double M2pp   = second_moment_call_(S0, K, mu_eff, sigma, T);
        double varP   = std::exp(-2*r*T) * std::exp(lam*lam*T) * M2pp - C*C;
        double ratio_theory_call = std::max(varP / varQ, 1e-5);

        lambdas    .push_back(lam);
        ratios     .push_back(std::max(pt.variance_ratio, 1e-5));
        theory_call.push_back(ratio_theory_call);
        theory_rn  .push_back(std::exp(lam*lam*T));
    }

    auto f = figure(true); f->size(1200, 680);
    auto ax = gca();
    hold(ax, true);

    { auto l = plot(ax, lambdas, theory_call, "-");  l->line_width(2.8); l->color("black"); }
    { auto l = plot(ax, lambdas, ratios, "o--"); l->line_width(2.2); l->color("blue"); l->marker_size(7); }
    { auto l = plot(ax, lambdas, theory_rn, ":");    l->line_width(1.8); l->color({0.4f, 0.4f, 0.4f, 1.0f}); }

    ax->y_axis().scale(axis_type::axis_scale::log);

    ax->title("Cost of changing measure: Var(P-MC) / Var(Q-MC) vs market price of risk lambda");
    ax->xlabel("lambda = (mu - r) / sigma");
    ax->ylabel("Variance ratio  [log scale]");
    legend(ax, {"theory (call payoff, exact)",
                "MC: (std_P / std_Q)^2",
                "exp(lambda^2 * T)  (var of dQ/dP, not of X_P)"});
    grid(ax, true);
    save(filename);
    cla();
}
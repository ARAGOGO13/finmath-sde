#pragma once
#include <cmath>
#include <matplot/matplot.h>
#include <string>
#include <vector>
#include "girsanov_experiment.hpp"

using namespace matplot;

// G1: три оценщика цены колла vs lambda = (mu-r)/sigma.
// Стиль: аналитика — чёрный solid, Q-MC — красный dashed + error bars,
//        P-MC — синий dashed + error bars (±1.96 * std_err).
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

    // Прокси-линии для легенды (3 слота)
    { auto l = plot(ax, lambdas, bs_vals, "-");  l->line_width(3.0); l->color("black"); }
    { auto l = plot(ax, lambdas, rn_vals, "--"); l->line_width(2.2); l->color("red");   }
    { auto l = plot(ax, lambdas, p_vals,  "--"); l->line_width(2.2); l->color("blue");  }

    // Error bars ±1.96σ для Q-MC (красный)
    {
        auto h = errorbar(ax, lambdas, rn_vals, rn_err,
                          error_bar::type::vertical, "--");
        h->color("red"); h->line_width(1.2); h->cap_size(4.f);
    }
    // Error bars ±1.96σ для P-MC (синий)
    {
        auto h = errorbar(ax, lambdas, p_vals, p_err,
                          error_bar::type::vertical, "--");
        h->color("blue"); h->line_width(1.2); h->cap_size(4.f);
    }

    // Параметры в заголовке (K и T из первой точки достаточно взять глобально)
    ax->title("Girsanov verification: three estimators of the same call price");
    ax->xlabel("lambda = (mu - r) / sigma   (market price of risk)");
    ax->ylabel("Call price");
    legend(ax, {"BS analytic (black)", "Q-MC risk-neutral (red)", "P-MC Girsanov (blue)"});
    grid(ax, true);
    save(filename);
    cla();
}

// Вспомогательная: E[(S_T - K)^+^2] под лог-нормальной S_T с эффективным дрейфом mu_eff.
// (Используется для точной теоретической кривой variance ratio.)
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

// G1: стоимость смены меры — отношение дисперсий vs lambda (лог. шкала по Y).
//
// Теоретическая кривая (для call-опциона) — АСИММЕТРИЧНАЯ:
//   ratio(lambda) = [exp(lambda^2 * T) * M2(2r - mu, sigma, T) - exp(2rT)*C^2]
//                 / [M2(r, sigma, T) - exp(2rT)*C^2]
//
// NB: exp(lambda^2*T) — дисперсия самой dQ/dP, но для взвешенного payoff'а колла
// формула другая из-за корреляции payoff(S_T) и weight(W_T).
// Для call при lambda>0 вес гасит высокие payoff'ы → variance падает (ratio<1).
// Для call при lambda<0 вес усиливает высокие payoff'ы → variance растёт (ratio>>1).
inline void plot_girsanov_variance(
        const std::vector<GirsanovExperiment::Point>& pts,
        const std::string& filename,
        double T = 1.0,
        double r = 0.05, double sigma = 0.20,
        double K = 100.0, double S0 = 100.0) {

    // BS price (та же для всех mu)
    const double k_inv_sqrt2 = 0.7071067811865475244;
    auto Phi = [&](double x) { return 0.5 * std::erfc(-x * k_inv_sqrt2); };
    double sqT = std::sqrt(T);
    double d1  = (std::log(S0/K) + (r + 0.5*sigma*sigma)*T) / (sigma*sqT);
    double d2  = d1 - sigma*sqT;
    double C   = S0 * Phi(d1) - K * std::exp(-r*T) * Phi(d2);

    // Var^Q(X_Q) — знаменатель
    double varQ = std::exp(-2*r*T) * second_moment_call_(S0, K, r, sigma, T) - C*C;

    std::vector<double> lambdas, ratios, theory_call, theory_rn;
    lambdas.reserve(pts.size());
    for (const auto& pt : pts) {
        double lam    = pt.lambda;
        double mu     = r + lam*sigma;
        double mu_eff = 2*r - mu;                                  // дрейф под Q' (для E[X_P^2])
        double M2pp   = second_moment_call_(S0, K, mu_eff, sigma, T);
        double varP   = std::exp(-2*r*T) * std::exp(lam*lam*T) * M2pp - C*C;
        double ratio_theory_call = std::max(varP / varQ, 1e-5);

        lambdas    .push_back(lam);
        ratios     .push_back(std::max(pt.variance_ratio, 1e-5));
        theory_call.push_back(ratio_theory_call);
        theory_rn  .push_back(std::exp(lam*lam*T));  // только для справки
    }

    auto f = figure(true); f->size(1200, 680);
    auto ax = gca();
    hold(ax, true);

    // Корректная теория для call-опциона (чёрный solid)
    { auto l = plot(ax, lambdas, theory_call, "-");  l->line_width(2.8); l->color("black"); }
    // MC-результат (синий o--)
    { auto l = plot(ax, lambdas, ratios, "o--"); l->line_width(2.2); l->color("blue"); l->marker_size(7); }
    // Референс: exp(lambda^2 * T) — дисперсия самой dQ/dP (серый пунктир)
    { auto l = plot(ax, lambdas, theory_rn, ":");    l->line_width(1.8); l->color({0.4f, 0.4f, 0.4f, 1.0f}); }

    // Логарифмическая шкала по Y
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
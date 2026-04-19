#pragma once
#include <matplot/matplot.h>
#include <string>
#include <utility>
#include <vector>
#include "convergence_analyzer.hpp"

// =============================================================================
// plotting.hpp — все функции визуализации для презентации
// =============================================================================

// ---- Траектории СДУ ---------------------------------------------------------
inline void plot_sde_paths(const std::vector<double> &time, const std::vector<std::vector<double>> &paths,
                           const std::string &title_str, const std::string &color = "blue") {
    using namespace matplot;
    figure(true);
    auto ax = gca();
    hold(ax, true);
    for (const auto &p: paths) {
        auto line = plot(ax, time, p);
        line->line_width(0.6);
        line->color(color);
    }
    ax->title(title_str);
    ax->xlabel("Время t");
    ax->ylabel("X_t");
    grid(ax, on);
}

// ---- Траектории с теоретическим средним ------------------------------------
// ИСПРАВЛЕНО:
// • траектория МК (представитель) и E[X_t] теперь на переднем плане (рисуются последними)
// • легенда по-прежнему точно соответствует цветам
inline void plot_paths_with_mean(const std::vector<double> &time, const std::vector<std::vector<double>> &paths,
                                 const std::vector<double> &mean_theory, const std::string &title_str) {
    using namespace matplot;
    figure(true);
    auto ax = gca();
    hold(ax, true);

    // 1. Dummy-линии для корректной легенды (рисуются первыми)
    if (!paths.empty()) {
        auto dummy_mk = plot(ax, time, paths[0]);
        dummy_mk->line_width(2.5);
        dummy_mk->color({0.2, 0.4, 0.8, 1.0});
    }
    auto dummy_mean = plot(ax, time, mean_theory, "-r");
    dummy_mean->line_width(2.5);

    // 2. Все полупрозрачные траектории МК — фон (будут под важными линиями)
    for (const auto &p: paths) {
        auto line = plot(ax, time, p);
        line->line_width(0.8);
        line->color({0.2, 0.4, 0.8, 0.5});
    }

    // 3. Реальные важные линии — НА ПЕРЕДНЕМ ПЛАНЕ (поверх всех)
    if (!paths.empty()) {
        auto real_mk = plot(ax, time, paths[0]);
        real_mk->line_width(2.5);
        real_mk->color({0.2, 0.4, 0.8, 1.0});
    }
    auto real_mean = plot(ax, time, mean_theory, "-r");
    real_mean->line_width(2.5);

    ax->title(title_str);
    ax->xlabel("Время t");
    ax->ylabel("X_t");
    legend({"Траектории МК", "E[X_t] теория"});
    grid(ax, on);
}

// ---- Процесс Пуассона -------------------------------------------------------
inline void plot_poisson(const std::vector<std::pair<double, int>> &events, double T) {
    using namespace matplot;
    figure(true);
    auto ax = gca();
    std::vector<double> t_vals = {0.0};
    std::vector<double> count_vals = {0.0};
    for (const auto &e: events) {
        t_vals.push_back(e.first);
        count_vals.push_back(static_cast<double>(e.second));
        t_vals.push_back(e.first);
        count_vals.push_back(static_cast<double>(e.second) - 1.0);
    }
    t_vals.push_back(T);
    count_vals.push_back(events.empty() ? 0.0 : static_cast<double>(events.back().second));

    auto st = stairs(ax, t_vals, count_vals);
    st->line_width(2.0);
    st->color("blue");
    ax->title("Процесс Пуассона — скачкообразный аналог броуновского движения (C.6)");
    ax->xlabel("Время t");
    ax->ylabel("N_t (число скачков)");
    grid(ax, on);
}

// ---- График сходимости Эйлер vs Мильштейн ----------------------------------
inline void plot_convergence(const ConvergenceResult &res, const std::string &filename) {
    using namespace matplot;
    figure(true);
    auto ax = gca();
    hold(ax, true);

    auto euler_line = loglog(ax, res.dt_values, res.euler_errors, "-ob");
    euler_line->line_width(2.0);
    euler_line->marker_size(8);

    auto milshtein_line = loglog(ax, res.dt_values, res.milshtein_errors, "-sr");
    milshtein_line->line_width(2.0);
    milshtein_line->marker_size(8);

    double dt0 = res.dt_values.front();
    double err0_e = res.euler_errors.front();
    double err0_m = res.milshtein_errors.front();

    std::vector<double> ref_x, ref_half, ref_one;
    for (double dt: res.dt_values) {
        ref_x.push_back(dt);
        ref_half.push_back(err0_e * std::pow(dt / dt0, 0.5));
        ref_one.push_back(err0_m * std::pow(dt / dt0, 1.0));
    }

    auto ref_half_line = loglog(ax, ref_x, ref_half, "--");
    ref_half_line->color("cyan");
    ref_half_line->line_width(1.2);

    auto ref_one_line = loglog(ax, ref_x, ref_one, "--");
    ref_one_line->color("magenta");
    ref_one_line->line_width(1.2);

    ax->title("Сходимость схем: Эйлер vs Мильштейн (GBM)");
    ax->xlabel("Шаг dt (log scale)");
    ax->ylabel("E[|X_T^exact - X_T^scheme|] (log scale)");
    legend({"Эйлер (наклон≈0.5)", "Мильштейн (наклон≈1.0)", "ref: O(dt^0.5)", "ref: O(dt^1.0)"});
    grid(ax, on);

    save(filename);
    cla();
}

// ---- Сравнение траекторий двух схем на одном графике -----------------------
// ИСПРАВЛЕНО: Эйлер, Мильштейн и E[X_t] теперь тоже на переднем плане
inline void plot_scheme_comparison(const std::vector<double> &time, const std::vector<std::vector<double>> &euler_paths,
                                   const std::vector<std::vector<double>> &milshtein_paths,
                                   const std::vector<double> &exact_mean, const std::string &title_str) {
    using namespace matplot;
    figure(true);
    auto ax = gca();
    hold(ax, true);

    // 1. Dummy-линии для легенды (будут скрыты)
    if (!euler_paths.empty()) {
        auto dummy_e = plot(ax, time, euler_paths[0]);
        dummy_e->line_width(2.5);
        dummy_e->color({0.2, 0.4, 0.9, 1.0});
    }
    if (!milshtein_paths.empty()) {
        auto dummy_m = plot(ax, time, milshtein_paths[0]);
        dummy_m->line_width(2.5);
        dummy_m->color({0.9, 0.3, 0.2, 1.0});
    }
    auto dummy_mean = plot(ax, time, exact_mean, "-k");
    dummy_mean->line_width(2.5);

    // 2. Все полупрозрачные траектории — фон
    for (const auto &p: euler_paths) {
        auto line = plot(ax, time, p);
        line->line_width(0.7);
        line->color({0.2, 0.4, 0.9, 0.4});
    }
    for (const auto &p: milshtein_paths) {
        auto line = plot(ax, time, p);
        line->line_width(0.7);
        line->color({0.9, 0.3, 0.2, 0.4});
    }

    // 3. Реальные важные линии — НА ПЕРЕДНЕМ ПЛАНЕ
    if (!euler_paths.empty()) {
        auto real_e = plot(ax, time, euler_paths[0]);
        real_e->line_width(2.5);
        real_e->color({0.2, 0.4, 0.9, 1.0});
    }
    if (!milshtein_paths.empty()) {
        auto real_m = plot(ax, time, milshtein_paths[0]);
        real_m->line_width(2.5);
        real_m->color({0.9, 0.3, 0.2, 1.0});
    }
    auto real_mean = plot(ax, time, exact_mean, "-k");
    real_mean->line_width(2.5);

    ax->title(title_str);
    ax->xlabel("Время t");
    ax->ylabel("X_t");
    legend({"Эйлер", "Мильштейн", "E[X_t] теория"});
    grid(ax, on);
}

// ---- Гистограмма распределения X_T ------------------------------------------
inline void plot_distribution(const std::vector<double> &samples, double theory_mean, double theory_std,
                              const std::string &title_str) {
    using namespace matplot;
    figure(true);
    auto ax = gca();

    hist(ax, samples, 50);
    ax->title(title_str);
    ax->xlabel("X_T");
    ax->ylabel("Частота");

    hold(ax, true);
    std::string full_title = title_str + "  |  E[X_T] теория=" + std::to_string(theory_mean).substr(0, 6) +
                             ", std=" + std::to_string(theory_std).substr(0, 5);
    ax->title(full_title);
    grid(ax, on);
}

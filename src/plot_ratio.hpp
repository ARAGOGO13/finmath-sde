#pragma once
#include <array>
#include <cmath>
#include <limits>
#include <matplot/matplot.h>
#include <string>
#include <vector>
#include "models/ratio_gbm.hpp"

using namespace matplot;

// Цветовая шкала cool-to-warm для rho = {-1, -0.5, 0, +0.5, +1}
static const std::vector<std::string> RHO_COLORS = {
    "blue", "green", "magenta", "orange", "red"
};
static const std::vector<std::string> RHO_LABELS = {
    "rho=-1.0", "rho=-0.5", "rho=0.0", "rho=+0.5", "rho=+1.0"
};

// Суффикс для имени выходного PDF-файла по значению rho
static std::string rho_tag(double rho) {
    if (rho < -0.9) return "m10";
    if (rho < -0.3) return "m05";
    if (rho <  0.3) return "00";
    if (rho <  0.8) return "p05";
    return "p10";
}

// Сохраняет по одному PDF на каждый rho: пути V_t (оранжевый, полупрозрачный),
// E[V_t] теория (чёрный), E[V_t] МК (красный пунктир), полоса ±1σ (чёрный пунктир).
//
// Все пять графиков имеют одинаковый диапазон Y — он вычисляется заранее по
// глобальному охвату теор. среднего ±1.5σ по всем rho. Благодаря этому разница
// в σ_V между rho=-1 и rho=+1 становится визуально очевидной.
//
// Паттерн двойной отрисовки: прокси-линии для легенды рисуются первыми (слоты 1–4),
// затем все пути, затем средние и границы повторно поверх путей для акцента.
inline void plot_ratio_paths(
        const std::vector<double>& time_grid,
        const std::vector<RatioGBMResult>& results,
        const std::string& base) {
    // Оранжевый, alpha=0.45: единый цвет по всем rho — плотность путей как единственная переменная
    static const std::array<float,4> path_col = {0.90f, 0.50f, 0.05f, 0.45f};
    size_t Ts = time_grid.size();

    // Глобальный диапазон Y по теор. ±1.5σ (не по сырым путям — у них возможны выбросы)
    double y_lo =  std::numeric_limits<double>::max();
    double y_hi = -std::numeric_limits<double>::max();
    for (const auto& r : results) {
        for (size_t t = 0; t < Ts; ++t) {
            double sd = std::sqrt(std::max(r.theory_var[t], 0.0));
            y_hi = std::max(y_hi, r.theory_mean[t] + 1.5*sd);
            y_lo = std::min(y_lo, std::max(r.theory_mean[t] - 1.5*sd, 0.0));
        }
    }
    double margin = 0.06*(y_hi - y_lo);
    y_lo = std::max(0.0, y_lo - margin);
    y_hi += margin;

    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];

        // Вычисляем границы ±1σ для текущего rho
        std::vector<double> upper(Ts), lower(Ts);
        for (size_t t = 0; t < Ts; ++t) {
            double sd = std::sqrt(std::max(r.theory_var[t], 0.0));
            upper[t] = r.theory_mean[t] + sd;
            lower[t] = std::max(r.theory_mean[t] - sd, 0.0);
        }

        auto f = figure(true); f->size(1100, 620);
        auto ax = gca();
        hold(ax, true);

        // Прокси-линии — слоты 1 (теория), 2 (МК), 3 (±1σ), 4 (пути)
        { auto l = plot(ax, time_grid, r.theory_mean, "-");  l->line_width(3.0); l->color("black"); }
        { auto l = plot(ax, time_grid, r.mc_mean,     "--"); l->line_width(2.0); l->color("red");   }
        { auto l = plot(ax, time_grid, upper,         "--"); l->line_width(1.8); l->color("black"); }
        { auto l = plot(ax, time_grid, r.paths[0]);          l->line_width(0.5); l->color(path_col); }

        // Остальные пути — слотов легенды не занимают
        for (size_t j = 1; j < r.paths.size(); ++j) {
            auto l = plot(ax, time_grid, r.paths[j]); l->line_width(0.5); l->color(path_col);
        }

        // Повторная отрисовка средних и границ поверх путей — новых слотов не добавляет
        { auto l = plot(ax, time_grid, lower,         "--"); l->line_width(1.8); l->color("black"); }
        { auto l = plot(ax, time_grid, upper,         "--"); l->line_width(1.8); l->color("black"); }
        { auto l = plot(ax, time_grid, r.theory_mean, "-");  l->line_width(3.0); l->color("black"); }
        { auto l = plot(ax, time_grid, r.mc_mean,     "--"); l->line_width(2.0); l->color("red");   }

        std::string title = RHO_LABELS[i]
            + "   mu_V="    + std::to_string(r.mu_V   ).substr(0,6)
            + "   sigma_V=" + std::to_string(r.sigma_V).substr(0,5);
        ax->title(title);
        ax->xlabel("t"); ax->ylabel("V_t");
        ylim(ax, {y_lo, y_hi});
        legend(ax, {"E[V_t] theory", "E[V_t] MC", "+/- 1 sigma", "paths"});
        grid(ax, true);

        save(base + "_" + rho_tag(r.rho) + ".pdf");
        cla();
    }
}

// Сохраняет один PDF: Var[V_t] для всех rho на одних осях.
// Теория (сплошная, цвет по rho) рисуется первой — получает 5 слотов легенды.
// МК (пунктир, тот же цвет) рисуется после — слотов легенды не получает,
// но идентифицируется по цвету совместно с теоретической кривой.
inline void plot_ratio_all_vars(
        const std::vector<double>& time_grid,
        const std::vector<RatioGBMResult>& results,
        const std::string& filename) {
    auto f = figure(true); f->size(1200, 650);
    auto ax = gca();
    hold(ax, true);

    // Теоретические кривые первыми — слоты легенды 1–5
    for (size_t i = 0; i < results.size(); ++i) {
        auto lt = plot(ax, time_grid, results[i].theory_var, "-");
        lt->line_width(2.8); lt->color(RHO_COLORS[i]);
    }
    // МК после — слотов легенды не получают (цвет совпадает с теорией)
    for (size_t i = 0; i < results.size(); ++i) {
        auto lm = plot(ax, time_grid, results[i].mc_var, "--");
        lm->line_width(1.5); lm->color(RHO_COLORS[i]);
    }

    ax->title("Var[V_t] for all rho   (solid = theory,  dashed = MC)");
    ax->xlabel("t"); ax->ylabel("Var[V_t]");
    legend(ax, {"rho = -1.0", "rho = -0.5", "rho =  0.0", "rho = +0.5", "rho = +1.0"});
    grid(ax, true);
    save(filename);
    cla();
}
#pragma once
#include <matplot/matplot.h>
#include <vector>
#include <string>
#include <cmath>
#include "models/ratio_gbm.hpp"

using namespace matplot;

// Cool-to-warm: blue→green→magenta→orange→red  (rho: -1→0→+1)
// Избегаем "cyan" (слишком светлый) и "black" (неразличим в полупрозрачных путях)
static const std::vector<std::string> RHO_COLORS = {
    "blue", "green", "magenta", "orange", "red"
};
static const std::vector<std::string> RHO_LABELS = {
    "rho=-1.0", "rho=-0.5", "rho=0.0", "rho=+0.5", "rho=+1.0"
};

static std::string rho_tag(double rho) {
    if (rho < -0.9) return "m10";
    if (rho < -0.3) return "m05";
    if (rho <  0.3) return "00";
    if (rho <  0.8) return "p05";
    return "p10";
}

// =============================================================================
// Plot A — пути + E[V_t] + ±1σ bounds: ОДИН файл на каждый rho
// =============================================================================
inline void plot_ratio_paths(
    const std::vector<double>& time_grid,
    const std::vector<RatioGBMResult>& results,
    const std::string& base)
{
    // Цвет как у rho=+0.5 (orange) для всех rho — наглядно, не конкурирует с mean-линиями
    // alpha=0.45 — меньше полупрозрачности (более плотные пути)
    static const std::vector<std::array<float,4>> path_cols = {
        {0.90f, 0.50f, 0.05f, 0.45f},
        {0.90f, 0.50f, 0.05f, 0.45f},
        {0.90f, 0.50f, 0.05f, 0.45f},
        {0.90f, 0.50f, 0.05f, 0.45f},
        {0.90f, 0.50f, 0.05f, 0.45f},
    };

    size_t Ts = time_grid.size();

    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];

        // Compute ±1σ theory bounds
        std::vector<double> upper(Ts), lower(Ts);
        for (size_t t = 0; t < Ts; ++t) {
            double sd = std::sqrt(std::max(r.theory_var[t], 0.0));
            upper[t] = r.theory_mean[t] + sd;
            lower[t] = std::max(r.theory_mean[t] - sd, 0.0);
        }

        auto f = figure(true); f->size(1100, 620);
        auto ax = gca();
        hold(ax, true);

        // ── Legend proxy lines (slots 1, 2, 3, 4) ──────────────────────────────
        // slot 1: theory mean
        auto lth_leg = plot(ax, time_grid, r.theory_mean, "-");
        lth_leg->line_width(3.0); lth_leg->color(RHO_COLORS[i]);
        // slot 2: MC mean
        auto lmc_leg = plot(ax, time_grid, r.mc_mean, "--");
        lmc_leg->line_width(2.0); lmc_leg->color("red");
        // slot 3: ±1σ bound
        auto lu_leg = plot(ax, time_grid, upper, "--");
        lu_leg->line_width(1.8); lu_leg->color(RHO_COLORS[i]);
        // slot 4: representative path
        {
            auto l = plot(ax, time_grid, r.paths[0]);
            l->line_width(0.5); l->color(path_cols[i]);
        }
        // ── Remaining paths (slots 5..N+4, no label) ────────────────────────
        for (size_t j = 1; j < r.paths.size(); ++j) {
            auto l = plot(ax, time_grid, r.paths[j]);
            l->line_width(0.5); l->color(path_cols[i]);
        }
        // ── Visual lines redrawn ON TOP (no new legend slots) ────────────────
        {   // lower bound
            auto ll = plot(ax, time_grid, lower, "--");
            ll->line_width(1.8); ll->color(RHO_COLORS[i]);
        }
        {   // upper bound
            auto lu = plot(ax, time_grid, upper, "--");
            lu->line_width(1.8); lu->color(RHO_COLORS[i]);
        }
        {   // theory mean
            auto lth = plot(ax, time_grid, r.theory_mean, "-");
            lth->line_width(3.0); lth->color(RHO_COLORS[i]);
        }
        {   // MC mean
            auto lmc = plot(ax, time_grid, r.mc_mean, "--");
            lmc->line_width(2.0); lmc->color("red");
        }

        std::string title = RHO_LABELS[i]
            + "   mu_V=" + std::to_string(r.mu_V   ).substr(0,6)
            + "   sigma_V=" + std::to_string(r.sigma_V).substr(0,5);
        ax->title(title);
        ax->xlabel("t"); ax->ylabel("V_t");
        // 4 labels → slots 1=theory(color), 2=MC(red), 3=±1σ(dashed), 4=paths(thin)
        legend(ax, {"E[V_t] theory", "E[V_t] MC", "+/- 1 sigma", "paths"});
        grid(ax, true);

        save(base + "_" + rho_tag(r.rho) + ".pdf");
        cla();
    }
}

// =============================================================================
// Plot B — Var[V_t]: все rho на одних осях
//   sigma_V: rho=-1 -> 0.35,  rho=0 -> 0.25,  rho=+1 -> 0.05  (разница ~58x!)
// =============================================================================
inline void plot_ratio_all_vars(
    const std::vector<double>& time_grid,
    const std::vector<RatioGBMResult>& results,
    const std::string& filename)
{
    auto f = figure(true); f->size(1200, 650);
    auto ax = gca();
    hold(ax, true);

    // Theory lines FIRST → legend entries 1-5 (correct colors)
    for (size_t i = 0; i < results.size(); ++i) {
        auto lt = plot(ax, time_grid, results[i].theory_var, "-");
        lt->line_width(2.8); lt->color(RHO_COLORS[i]);
    }
    // MC dashed AFTER → unlabeled (only 5 legend entries assigned above)
    for (size_t i = 0; i < results.size(); ++i) {
        auto lm = plot(ax, time_grid, results[i].mc_var, "--");
        lm->line_width(1.5); lm->color(RHO_COLORS[i]);
    }

    ax->title("Var[V_t] for all rho   (solid = theory,  dashed = MC)");
    ax->xlabel("t"); ax->ylabel("Var[V_t]");
    legend(ax, {"rho = -1.0", "rho = -0.5", "rho =  0.0",
                "rho = +0.5", "rho = +1.0"});
    grid(ax, true);
    save(filename);
    cla();
}
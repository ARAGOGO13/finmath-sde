// =============================================================================
// Визуализация модели отношения двух коррелированных GBM (RatioGBM).
//
// Для каждого значения корреляции ρ строятся:
//   • траектории процесса V_t с наложением теоретического среднего и ±1σ;
//   • сводный график дисперсии Var[V_t] для всех ρ.
//
// Также содержит вспомогательные функции для вывода таблицы результатов
// и форматирования тегов имён файлов.
// =============================================================================

#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <matplot/matplot.h>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include "models/ratio_gbm.hpp"
#include "console_utils.hpp"

using namespace matplot;

// Цвета для различных значений корреляции
static const std::vector<std::string> RHO_COLORS = {
    "blue", "green", "magenta", "orange", "red"
};
static const std::vector<std::string> RHO_LABELS = {
    "rho=-1.0", "rho=-0.5", "rho=0.0", "rho=+0.5", "rho=+1.0"
};

// Строковый тег для значения ρ (используется в именах файлов)
static std::string rho_tag(double rho) {
    if (rho < -0.9) return "m10";
    if (rho < -0.3) return "m05";
    if (rho <  0.3) return "00";
    if (rho <  0.8) return "p05";
    return "p10";
}

// Заголовок таблицы результатов для отношения GBM
inline void ratio_table_header() {
    std::vector<std::string> headers = {"rho", "mu_V", "sig_V", "E[V_T] MC", "E[V_T] th", "Err%"};
    std::vector<int> widths = {7, 12, 12, 14, 14, 8};
    print_table_header(headers, widths);
}

// Строка таблицы с результатами для одного ρ
inline void ratio_table_row(const RatioGBMResult &r, double theory_E) {
    double err = std::abs(r.mc_mean.back() - theory_E) / (std::abs(theory_E) + 1e-15) * 100.0;
    std::ostringstream ss_rho, ss_mu, ss_sig, ss_mc, ss_th, ss_err;
    ss_rho << std::fixed << std::setprecision(2) << r.rho;
    ss_mu  << std::setprecision(5) << r.mu_V;
    ss_sig << std::setprecision(5) << r.sigma_V;
    ss_mc  << r.mc_mean.back();
    ss_th  << theory_E;
    ss_err << std::setprecision(3) << err;
    print_table_row_with_error({ss_rho.str(), ss_mu.str(), ss_sig.str(),
                                ss_mc.str(), ss_th.str(), ss_err.str()},
                               {7, 12, 12, 14, 14, 8}, err);
}

// Графики траекторий V_t для каждого значения ρ (отдельные PDF)
inline void plot_ratio_paths(
        const std::vector<double>& time_grid,
        const std::vector<RatioGBMResult>& results,
        const std::string& base) {
    static const std::array<float,4> path_col = {0.90f, 0.50f, 0.05f, 0.45f};
    size_t Ts = time_grid.size();

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

        std::vector<double> upper(Ts), lower(Ts);
        for (size_t t = 0; t < Ts; ++t) {
            double sd = std::sqrt(std::max(r.theory_var[t], 0.0));
            upper[t] = r.theory_mean[t] + sd;
            lower[t] = std::max(r.theory_mean[t] - sd, 0.0);
        }

        auto f = figure(true); f->size(1100, 620);
        auto ax = gca();
        hold(ax, true);

        { auto l = plot(ax, time_grid, r.theory_mean, "-");  l->line_width(3.0); l->color("black"); }
        { auto l = plot(ax, time_grid, r.mc_mean,     "--"); l->line_width(2.0); l->color("red");   }
        { auto l = plot(ax, time_grid, upper,         "--"); l->line_width(1.8); l->color("black"); }
        { auto l = plot(ax, time_grid, r.paths[0]);          l->line_width(0.5); l->color(path_col); }

        for (size_t j = 1; j < r.paths.size(); ++j) {
            auto l = plot(ax, time_grid, r.paths[j]); l->line_width(0.5); l->color(path_col);
        }

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

// Сводный график дисперсии Var[V_t] для всех ρ на одном полотне
inline void plot_ratio_all_vars(
        const std::vector<double>& time_grid,
        const std::vector<RatioGBMResult>& results,
        const std::string& filename) {
    auto f = figure(true); f->size(1200, 650);
    auto ax = gca();
    hold(ax, true);

    for (size_t i = 0; i < results.size(); ++i) {
        auto lt = plot(ax, time_grid, results[i].theory_var, "-");
        lt->line_width(2.8); lt->color(RHO_COLORS[i]);
    }
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
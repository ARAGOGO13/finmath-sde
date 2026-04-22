#pragma once
#include <cmath>
#include <matplot/matplot.h>
#include <string>
#include <vector>
#include "hedging_experiment.hpp"

using namespace matplot;

// Цвета для гистограмм (по числу N_rebal вариантов)
static const std::vector<std::string> HIST_COLORS = {"blue", "red", "green", "orange"};

// G2a: Наложенные гистограммы P&L для нескольких N_rebal.
// Нормировка: probability (сумма столбцов ≤ 1).
inline void plot_hedging_pnl_histograms(
        const std::vector<HedgingPnLResult>& results,
        const std::vector<std::string>& labels,
        const std::string& title_str,
        const std::string& filename) {

    auto f = figure(true); f->size(1200, 680);
    auto ax = gca();
    hold(ax, true);

    for (size_t i = 0; i < results.size(); ++i) {
        const std::string& col = HIST_COLORS[i % HIST_COLORS.size()];
        auto h = hist(ax, results[i].pnl_samples, 80UL);
        h->normalization(histogram::normalization::probability);
        h->face_color(col);
        h->face_alpha(0.45f);
        h->edge_color("none");
    }

    ax->title(title_str);
    ax->xlabel("Hedging P&L");
    ax->ylabel("Probability");
    legend(ax, labels);
    grid(ax, true);
    save(filename);
    cla();
}

// G2b: std(P&L) vs dt в log-log, наклон ≈ +0.5 (теорема Bertsimas-Kogan-Lo).
// Аналог plot_convergence в plotting.hpp.
inline void plot_hedging_std_vs_dt(
        const DeltaHedgingExperiment::ScanResult& scan,
        const std::string& filename) {

    const int n = static_cast<int>(scan.dt_list.size());
    std::vector<double> log_dt(n), log_std(n);
    for (int i = 0; i < n; ++i) {
        log_dt[i]  = std::log(scan.dt_list[i]);
        log_std[i] = std::log(std::max(scan.std_pnl[i], 1e-15));
    }

    // Сдвиг b для теор. прямой slope=0.5: min MSE через mean
    double b_ref = 0.0, b_fit = 0.0;
    for (int i = 0; i < n; ++i) {
        b_ref += log_std[i] - 0.5     * log_dt[i];
        b_fit += log_std[i] - scan.slope_log_std_vs_log_dt * log_dt[i];
    }
    b_ref /= n; b_fit /= n;

    double x0 = log_dt.front(), x1 = log_dt.back();
    std::vector<double> tx = {x0, x1};
    std::vector<double> ref_line = {0.5*x0 + b_ref, 0.5*x1 + b_ref};
    std::vector<double> fit_line = {scan.slope_log_std_vs_log_dt*x0 + b_fit,
                                    scan.slope_log_std_vs_log_dt*x1 + b_fit};

    // Форматируем slope для заголовка
    double sl = scan.slope_log_std_vs_log_dt;
    int sl_int = static_cast<int>(sl * 1000);
    std::string slope_str = std::to_string(sl_int / 1000) + "."
                          + std::to_string(std::abs(sl_int) % 1000);

    auto f = figure(true); f->size(1100, 680);
    auto ax = gca();
    hold(ax, true);

    // Данные — оранжевые точки
    auto ld = plot(ax, log_dt, log_std, "o");
    ld->line_width(2.5); ld->color("orange"); ld->marker_size(9);

    // Линейная регрессия — оранжевый пунктир
    auto lf = plot(ax, tx, fit_line, "--");
    lf->line_width(2.0); lf->color("orange");

    // Эталонная прямая slope=+0.5 — чёрный пунктир
    auto lr = plot(ax, tx, ref_line, "--");
    lr->line_width(2.0); lr->color("black");

    ax->title("Discrete hedging: std(P&L) ~ sqrt(dt)   [slope=" + slope_str
              + ",  theory 0.5]");
    ax->xlabel("log(dt)   [dt = T / N_rebal]");
    ax->ylabel("log( std(P&L) )");
    legend(ax, {"MC std(P&L)", "OLS fit (slope=" + slope_str + ")",
                "slope +0.5 (theory)"});
    grid(ax, true);
    save(filename);
    cla();
}

// G2c: Эффект несоответствия волатильностей.
// Два графика (subplot 1x2): mean(P&L) и std(P&L) vs sigma_real.
// Теоретическая кривая: E[P&L] ≈ 0.5 * Γ * S0² * (sigma_impl² - sigma_real²) * T
// (Γ = phi(d1) / (S0 * sigma_impl * sqrt(T)), вычислено при sigma=sigma_impl, K=S0).
inline void plot_hedging_vol_mismatch(
        const std::vector<HedgingPnLResult>& results,
        const std::vector<double>& sigma_real_list,
        double sigma_impl,
        double S0, double K, double T, double r,
        const std::string& filename) {

    // Вычисляем Γ по sigma_impl при S0=K (ATM)
    const double sqrt_T = std::sqrt(T);
    const double d1_atm = ((r + 0.5*sigma_impl*sigma_impl)*T)
                         / (sigma_impl * sqrt_T);
    const double k_sqrt2pi = 2.5066282746310002416;
    const double phi_d1    = std::exp(-0.5 * d1_atm * d1_atm) / k_sqrt2pi;
    const double gamma_atm = phi_d1 / (S0 * sigma_impl * sqrt_T);

    // Теоретическая поправка: 0.5*Γ*S0²*(sigma_impl² - sigma_real²)*T
    const int n = static_cast<int>(sigma_real_list.size());
    std::vector<double> mean_vals(n), std_vals(n), theory_mean(n);
    for (int i = 0; i < n; ++i) {
        mean_vals[i]  = results[i].mean_pnl;
        std_vals[i]   = results[i].std_pnl;
        theory_mean[i] = 0.5 * gamma_atm * S0 * S0
                       * (sigma_impl*sigma_impl - sigma_real_list[i]*sigma_real_list[i]) * T;
    }

    auto f = figure(true); f->size(2200, 700);

    // ---- Левая панель: mean(P&L) ----
    {
        auto ax = f->add_subplot(1, 2, 0);
        hold(ax, true);

        // MC mean — красный пунктир с маркерами
        auto lm = plot(ax, sigma_real_list, mean_vals, "o--");
        lm->line_width(2.2); lm->color("red"); lm->marker_size(8);

        // Теория — чёрный solid
        auto lt = plot(ax, sigma_real_list, theory_mean, "-");
        lt->line_width(2.8); lt->color("black");

        // Вертикальная линия на sigma_impl
        std::vector<double> vx = {sigma_impl, sigma_impl};
        double ylo = *std::min_element(theory_mean.begin(), theory_mean.end());
        double yhi = *std::max_element(theory_mean.begin(), theory_mean.end());
        double margin = 0.1 * std::abs(yhi - ylo) + 0.5;
        std::vector<double> vy = {ylo - margin, yhi + margin};
        auto lv = plot(ax, vx, vy, "--");
        lv->line_width(1.5); lv->color("blue");

        ax->title("Volatility mis-specification: E[P&L] vs sigma_real  (sigma_impl="
                  + std::to_string(sigma_impl).substr(0,4) + ")");
        ax->xlabel("sigma_real");
        ax->ylabel("E[ P&L ]");
        legend(ax, {"MC mean", "theory 0.5*Gamma*S0^2*(sig_impl^2-sig_r^2)*T",
                    "sigma_impl (hedge vol)"});
        grid(ax, true);
    }

    // ---- Правая панель: std(P&L) ----
    {
        auto ax = f->add_subplot(1, 2, 1);
        hold(ax, true);

        auto lm = plot(ax, sigma_real_list, std_vals, "o--");
        lm->line_width(2.2); lm->color("red"); lm->marker_size(8);

        // Вертикальная линия
        std::vector<double> vx = {sigma_impl, sigma_impl};
        double ylo = *std::min_element(std_vals.begin(), std_vals.end());
        double yhi = *std::max_element(std_vals.begin(), std_vals.end());
        double margin = 0.05 * (yhi - ylo) + 0.02;
        std::vector<double> vy = {std::max(0.0, ylo - margin), yhi + margin};
        auto lv = plot(ax, vx, vy, "--");
        lv->line_width(1.5); lv->color("blue");

        ax->title("Volatility mis-specification: std(P&L) vs sigma_real  (sigma_impl="
                  + std::to_string(sigma_impl).substr(0,4) + ")");
        ax->xlabel("sigma_real");
        ax->ylabel("std( P&L )");
        legend(ax, {"MC std", "sigma_impl (hedge vol)"});
        grid(ax, true);
    }

    f->save(filename);
    cla();
}
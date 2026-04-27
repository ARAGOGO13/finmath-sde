// =============================================================================
// Универсальные функции для визуализации результатов Монте-Карло.
//
// Включает:
//   • построение траекторий с наложением среднего и теоретической кривой,
//   • графики дисперсии,
//   • анализ сильной сходимости численных схем,
//   • комплексную визуализацию для GBM, InverseGBM, Vasicek, CIR, Poisson.
//
// Все графики сохраняются в формате PDF.
// =============================================================================

#pragma once
#include <array>
#include <cmath>
#include <matplot/matplot.h>
#include <string>
#include <vector>
#include "convergence_analyzer.hpp"

using namespace matplot;

// Вычисление выборочного среднего и дисперсии по набору траекторий.
// paths – матрица траекторий (M × N); mc_mean, mc_var – выходные векторы длины N.
inline void compute_mc_moments(const std::vector<std::vector<double>> &paths,
                               std::vector<double> &mc_mean,
                               std::vector<double> &mc_var) {
    if (paths.empty()) return;
    size_t N = paths[0].size(), M = paths.size();
    mc_mean.assign(N, 0.0);
    mc_var.assign(N, 0.0);
    for (size_t t = 0; t < N; ++t) {
        double s = 0.0, s2 = 0.0;
        for (size_t j = 0; j < M; ++j) { s += paths[j][t]; s2 += paths[j][t]*paths[j][t]; }
        mc_mean[t] = s / (double)M;
        mc_var[t]  = s2 / (double)M - mc_mean[t]*mc_mean[t];
    }
}

// Рисование траекторий, MC-среднего и теоретического среднего на одном графике
static void draw_paths_mean(axes_handle ax,
                            const std::vector<double> &time,
                            const std::vector<std::vector<double>> &paths,
                            const std::vector<double> &mc_mean,
                            const std::vector<double> &th_mean,
                            const std::string &title_str,
                            const std::string &ylabel_str,
                            std::array<float, 4> path_color,
                            std::array<double, 2> y_lim = {0.0, 0.0}) {
    hold(ax, true);

    // Прокси-линии для легенды
    auto lmc_p = plot(ax, time, mc_mean, "--"); lmc_p->line_width(2.8); lmc_p->color("red");
    auto lth_p = plot(ax, time, th_mean, "-");  lth_p->line_width(2.8); lth_p->color("black");

    { auto l = plot(ax, time, paths[0]); l->line_width(0.4); l->color(path_color); }
    for (size_t i = 1; i < paths.size(); ++i) {
        auto l = plot(ax, time, paths[i]); l->line_width(0.4); l->color(path_color);
    }

    // Повтор поверх путей для легенды
    { auto l = plot(ax, time, mc_mean, "--"); l->line_width(2.8); l->color("red"); }
    { auto l = plot(ax, time, th_mean, "-");  l->line_width(2.8); l->color("black"); }

    ax->title(title_str); ax->xlabel("t"); ax->ylabel(ylabel_str);
    legend(ax, {"MC mean", "theory", "paths"});
    grid(ax, true);
    if (y_lim[1] > y_lim[0]) ylim(ax, {y_lim[0], y_lim[1]});
}

// Рисование дисперсии: MC и теория
static void draw_var(axes_handle ax,
                     const std::vector<double> &time,
                     const std::vector<double> &mc_var,
                     const std::vector<double> &th_var,
                     const std::string &title_str,
                     const std::string &ylabel_str) {
    hold(ax, true);
    auto lmc = plot(ax, time, mc_var, "--"); lmc->line_width(2.8); lmc->color("red");
    auto lth = plot(ax, time, th_var, "-");  lth->line_width(2.8); lth->color("black");
    ax->title(title_str); ax->xlabel("t"); ax->ylabel(ylabel_str);
    legend(ax, {"MC", "theory"});
    grid(ax, true);
}

// Комплексная визуализация GBM и InverseGBM (четыре графика)
inline void plot_gbm_and_inverse(
        const std::vector<double> &time,
        const std::vector<std::vector<double>> &gbm_paths,
        const std::vector<double> &gbm_mean_th, const std::vector<double> &gbm_var_th,
        const std::vector<std::vector<double>> &inv_paths,
        const std::vector<double> &inv_mean_th, const std::vector<double> &inv_var_th,
        const std::string &sde_gbm, const std::string &sde_inv,
        std::array<double,2> inv_mean_ylim,
        const std::string &base) {
    std::vector<double> gmc_m, gmc_v, imc_m, imc_v;
    compute_mc_moments(gbm_paths, gmc_m, gmc_v);
    compute_mc_moments(inv_paths, imc_m, imc_v);

    { auto f = figure(true); f->size(1000, 600);
      draw_paths_mean(gca(), time, gbm_paths, gmc_m, gbm_mean_th,
                      "GBM   E[S_t]   " + sde_gbm, "S_t",
                      {0.9f, 0.35f, 0.1f, 0.35f});
      save(base + "_mean.pdf"); cla(); }

    { auto f = figure(true); f->size(1000, 550);
      draw_var(gca(), time, gmc_v, gbm_var_th,
               "GBM   Var[S_t]   " + sde_gbm, "Var[S_t]");
      save(base + "_var.pdf"); cla(); }

    { auto f = figure(true); f->size(1000, 600);
      draw_paths_mean(gca(), time, inv_paths, imc_m, inv_mean_th,
                      "InvGBM (U=1/S)   E[U_t]   " + sde_inv, "U_t",
                      {0.9f, 0.35f, 0.1f, 0.35f}, inv_mean_ylim);
      save(base + "_inv_mean.pdf"); cla(); }

    { auto f = figure(true); f->size(1000, 550);
      draw_var(gca(), time, imc_v, inv_var_th,
               "InvGBM (U=1/S)   Var[U_t]   " + sde_inv, "Var[U_t]");
      save(base + "_inv_var.pdf"); cla(); }
}

// Визуализация сильной сходимости: log ошибки от log N
inline void plot_convergence(const ConvergenceResult &res, const std::string &filename) {
    size_t n = res.log_dt.size();
    std::vector<double> log_N(n);
    for (size_t i = 0; i < n; ++i) log_N[i] = -res.log_dt[i];

    auto f = figure(true); f->size(1000, 680);
    auto ax = gca();
    hold(ax, true);

    auto le = plot(ax, log_N, res.log_euler,     "o-");
    le->line_width(2.5); le->color("orange"); le->marker_size(9);
    auto lm = plot(ax, log_N, res.log_milshtein, "s-");
    lm->line_width(2.5); lm->color("green");  lm->marker_size(9);

    // Опорные прямые с наклонами -0.5 и -1.0
    double b_e = 0.0, b_m = 0.0;
    for (size_t i = 0; i < n; ++i) {
        b_e += res.log_euler[i]     + 0.5 * log_N[i];
        b_m += res.log_milshtein[i] + 1.0 * log_N[i];
    }
    b_e /= (double)n; b_m /= (double)n;
    double x0 = log_N.front(), x1 = log_N.back();
    std::vector<double> tx = {x0, x1};
    auto lte = plot(ax, tx, std::vector<double>{-0.5*x0+b_e, -0.5*x1+b_e}, "--");
    lte->line_width(1.8); lte->color("orange");
    auto ltm = plot(ax, tx, std::vector<double>{-1.0*x0+b_m, -1.0*x1+b_m}, "--");
    ltm->line_width(1.8); ltm->color("green");

    ax->title("Strong convergence:  Euler vs Milshtein   (x=log N,  y=log error)");
    ax->xlabel("log(N)   [ N = 8, 16, 32, ..., 4096 ]");
    ax->ylabel("log( E[ |X_exact - X_scheme| ] )");
    legend(ax, {"Euler", "Milshtein",
                "slope -0.5  (Euler theory)", "slope -1.0  (Milshtein theory)"});
    grid(ax, true);
    save(filename);
    cla();
}

// Визуализация процессов с возвратом к среднему (Vasicek, CIR)
inline void plot_mean_reversion(
        const std::vector<double> &time,
        const std::vector<std::vector<double>> &paths,
        const std::vector<double> &mean_th, const std::vector<double> &var_th,
        const std::string &model_name, const std::string &sde_str,
        std::array<double,2> mean_ylim,
        const std::string &base) {
    std::vector<double> mc_mean, mc_var;
    compute_mc_moments(paths, mc_mean, mc_var);

    { auto f = figure(true); f->size(1000, 620);
      draw_paths_mean(gca(), time, paths, mc_mean, mean_th,
                      model_name + "   E[X_t]   " + sde_str, "X_t",
                      {0.9f, 0.35f, 0.1f, 0.35f}, mean_ylim);
      save(base + "_mean.pdf"); cla(); }

    { auto f = figure(true); f->size(1000, 560);
      draw_var(gca(), time, mc_var, var_th,
               model_name + "   Var[X_t]   " + sde_str, "Var[X_t]");
      save(base + "_var.pdf"); cla(); }
}

// Визуализация пуассоновского процесса (три графика: пути, среднее, дисперсия)
inline void plot_poisson_full(const std::vector<double> &time_grid,
                              const std::vector<std::vector<std::vector<double>>> &all_paths,
                              const std::vector<double> &lambdas,
                              const std::string &base) {
    std::vector<std::string> cols   = {"blue", "red", "green"};
    std::vector<std::string> labels = {"lambda=2", "lambda=5", "lambda=15"};
    size_t T_steps = time_grid.size();

    { auto f = figure(true); f->size(1100, 500);
      auto ax = gca(); hold(ax, true);
      for (size_t li = 0; li < lambdas.size(); ++li) {
          if (all_paths[li].empty()) continue;
          auto l = stairs(ax, time_grid, all_paths[li][0]);
          l->line_width(2.2); l->color(cols[li]);
      }
      ax->title("Poisson sample paths");
      ax->xlabel("t"); ax->ylabel("N_t");
      legend(ax, labels); grid(ax, true);
      save(base + "_paths.pdf"); cla(); }

    std::vector<std::vector<double>> mc_means(lambdas.size()), mc_vars(lambdas.size());
    for (size_t li = 0; li < lambdas.size(); ++li) {
        size_t M = all_paths[li].size();
        mc_means[li].assign(T_steps, 0.0);
        mc_vars[li].assign(T_steps, 0.0);
        for (size_t j = 0; j < M; ++j)
            for (size_t t = 0; t < T_steps; ++t)
                mc_means[li][t] += all_paths[li][j][t];
        for (auto &v : mc_means[li]) v /= (double)M;
        for (size_t j = 0; j < M; ++j)
            for (size_t t = 0; t < T_steps; ++t) {
                double d = all_paths[li][j][t] - mc_means[li][t];
                mc_vars[li][t] += d*d;
            }
        for (auto &v : mc_vars[li]) v /= (double)M;
    }

    auto draw_moment = [&](const std::string &title_str, const std::string &ylabel_str,
                           const std::string &filename,
                           const std::vector<std::vector<double>> &mc_vals) {
        auto f = figure(true); f->size(1100, 520);
        auto ax = gca(); hold(ax, true);
        for (size_t li = 0; li < lambdas.size(); ++li) {
            std::vector<double> theory(T_steps);
            for (size_t t = 0; t < T_steps; ++t) theory[t] = lambdas[li]*time_grid[t];
            auto lm = plot(ax, time_grid, mc_vals[li], "--"); lm->line_width(2.2); lm->color(cols[li]);
            auto lt = plot(ax, time_grid, theory,      "-");  lt->line_width(1.6); lt->color(cols[li]);
        }
        ax->title(title_str); ax->xlabel("t"); ax->ylabel(ylabel_str);
        legend(ax, {"MC l=2","th l=2","MC l=5","th l=5","MC l=15","th l=15"});
        grid(ax, true); save(filename); cla();
    };

    draw_moment("E[N_t]:  MC (dashed)  vs  lambda*t (solid)",
                "E[N_t]", base + "_mean.pdf", mc_means);
    draw_moment("Var[N_t]:  MC (dashed)  vs  lambda*t (solid)    [E=Var for Poisson]",
                "Var[N_t]", base + "_var.pdf", mc_vars);
}

// Вычисление общего диапазона Y для графиков средних (InvGBM, Vasicek, CIR)
inline std::array<double,2> compute_shared_mean_ylim(double T) {
    const double U0_a = 1.0/100.0;
    const double sig_a = 0.20;
    const double mu_star_a = sig_a * sig_a - 0.05;
    double sd_inv = U0_a * std::sqrt(std::exp(sig_a * sig_a * T) - 1.0);
    double inv_hi = U0_a * std::exp(mu_star_a * T) + 3.0 * sd_inv;
    double inv_lo = std::max(0.0, U0_a * std::exp(mu_star_a * T) - 3.0 * sd_inv);

    const double theta_a = 2.0;
    const double mu_a = 0.05;
    const double sig_vc = 0.02;
    const double x0_a = 0.03;
    double sigma_stat = sig_vc / std::sqrt(2.0 * theta_a);
    double vc_hi = mu_a + 3.0 * sigma_stat;
    double vc_lo = std::max(0.0, x0_a - 3.0 * sigma_stat);

    double y_lo = std::min(inv_lo, vc_lo);
    double y_hi = std::max(inv_hi, vc_hi);
    double margin = 0.08 * (y_hi - y_lo);
    return {std::max(0.0, y_lo - margin), y_hi + margin};
}
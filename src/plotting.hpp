#pragma once
#include <cmath>
#include <matplot/matplot.h>
#include <string>
#include <vector>
#include "convergence_analyzer.hpp"

using namespace matplot;

inline void compute_mc_moments(const std::vector<std::vector<double>> &paths, std::vector<double> &mc_mean,
                               std::vector<double> &mc_var) {
    if (paths.empty())
        return;
    size_t N = paths[0].size();
    size_t M = paths.size();
    mc_mean.assign(N, 0.0);
    mc_var.assign(N, 0.0);
    for (size_t t = 0; t < N; ++t) {
        double s = 0.0, s2 = 0.0;
        for (size_t j = 0; j < M; ++j) {
            s += paths[j][t];
            s2 += paths[j][t] * paths[j][t];
        }
        mc_mean[t] = s / M;
        mc_var[t] = s2 / M - mc_mean[t] * mc_mean[t];
    }
}

static void draw_paths_mean(axes_handle ax, const std::vector<double> &time,
                            const std::vector<std::vector<double>> &paths, const std::vector<double> &mc_mean,
                            const std::vector<double> &th_mean, const std::string &title_str,
                            const std::string &ylabel_str, std::array<float, 4> path_color) {
    hold(ax, true);
    for (const auto &p: paths) {
        auto l = plot(ax, time, p);
        l->line_width(0.4);
        l->color(path_color);
    }
    auto lmc = plot(ax, time, mc_mean, "--");
    lmc->line_width(2.8);
    lmc->color("red");
    auto lth = plot(ax, time, th_mean, "-");
    lth->line_width(2.8);
    lth->color("black");
    ax->title(title_str);
    ax->xlabel("t");
    ax->ylabel(ylabel_str);
    legend(ax, {"paths", "MC mean", "theory"});
    grid(ax, true);
}

static void draw_var(axes_handle ax, const std::vector<double> &time, const std::vector<double> &mc_var,
                     const std::vector<double> &th_var, const std::string &title_str, const std::string &ylabel_str) {
    hold(ax, true);
    auto lmc = plot(ax, time, mc_var, "--");
    lmc->line_width(2.8);
    lmc->color("red");
    auto lth = plot(ax, time, th_var, "-");
    lth->line_width(2.8);
    lth->color("black");
    ax->title(title_str);
    ax->xlabel("t");
    ax->ylabel(ylabel_str);
    legend(ax, {"MC", "theory"});
    grid(ax, true);
}


inline void plot_gbm_and_inverse(const std::vector<double> &time, const std::vector<std::vector<double>> &gbm_paths,
                                 const std::vector<double> &gbm_mean_th, const std::vector<double> &gbm_var_th,
                                 const std::vector<std::vector<double>> &inv_paths,
                                 const std::vector<double> &inv_mean_th, const std::vector<double> &inv_var_th,
                                 const std::string &base) {
    std::vector<double> gmc_m, gmc_v, imc_m, imc_v;
    compute_mc_moments(gbm_paths, gmc_m, gmc_v);
    compute_mc_moments(inv_paths, imc_m, imc_v);

    {
        auto f = figure(true);
        f->size(1000, 600);
        draw_paths_mean(gca(), time, gbm_paths, gmc_m, gbm_mean_th, "GBM   E[S_t]   MC vs Theory", "S_t",
                        {0.2f, 0.45f, 0.9f, 0.15f});
        save(base + "_mean.pdf");
        cla();
    }
    {
        auto f = figure(true);
        f->size(1000, 550);
        draw_var(gca(), time, gmc_v, gbm_var_th, "GBM   Var[S_t]   MC vs Theory", "Var[S_t]");
        save(base + "_var.pdf");
        cla();
    }
    {
        auto f = figure(true);
        f->size(1000, 600);
        draw_paths_mean(gca(), time, inv_paths, imc_m, inv_mean_th, "InverseGBM (U=1/S)   E[U_t]   MC vs Theory", "U_t",
                        {0.9f, 0.35f, 0.1f, 0.15f});
        save(base + "_inv_mean.pdf");
        cla();
    }
    {
        auto f = figure(true);
        f->size(1000, 550);
        draw_var(gca(), time, imc_v, inv_var_th, "InverseGBM (U=1/S)   Var[U_t]   MC vs Theory", "Var[U_t]");
        save(base + "_inv_var.pdf");
        cla();
    }
}

inline void plot_convergence(const ConvergenceResult &res, const std::string &filename) {
    auto f = figure(true);
    f->size(1000, 680);
    auto ax = gca();
    hold(ax, true);

    auto le = plot(ax, res.log_dt, res.log_euler, "o-");
    le->line_width(2.5);
    le->color("orange");
    le->marker_size(9);

    auto lm = plot(ax, res.log_dt, res.log_milshtein, "s-");
    lm->line_width(2.5);
    lm->color("green");
    lm->marker_size(9);

    size_t n = res.log_dt.size();
    double b_e = 0.0, b_m = 0.0;
    for (size_t i = 0; i < n; ++i) {
        b_e += res.log_euler[i] - 0.5 * res.log_dt[i];
        b_m += res.log_milshtein[i] - 1.0 * res.log_dt[i];
    }
    b_e /= n;
    b_m /= n;
    double x0 = res.log_dt.front(), x1 = res.log_dt.back();
    std::vector<double> tx = {x0, x1};
    std::vector<double> tye = {0.5 * x0 + b_e, 0.5 * x1 + b_e};
    std::vector<double> tym = {1.0 * x0 + b_m, 1.0 * x1 + b_m};

    auto lte = plot(ax, tx, tye, "--");
    lte->line_width(1.8);
    lte->color("orange");
    auto ltm = plot(ax, tx, tym, "--");
    ltm->line_width(1.8);
    ltm->color("green");

    ax->title("Strong convergence:  Euler vs Milshtein  (log-log)");
    ax->xlabel("log(dt)");
    ax->ylabel("log( E[|error|] )");
    legend(ax, {"Euler (measured)", "Milshtein (measured)", "slope 0.5 (theory)", "slope 1.0 (theory)"});
    grid(ax, true);
    save(filename);
    cla();
}


inline void plot_mean_reversion(const std::vector<double> &time, const std::vector<std::vector<double>> &paths,
                                const std::vector<double> &mean_th, const std::vector<double> &var_th,
                                const std::string &model_name, const std::string &base) {
    std::vector<double> mc_mean, mc_var;
    compute_mc_moments(paths, mc_mean, mc_var);

    {
        auto f = figure(true);
        f->size(1000, 620);
        draw_paths_mean(gca(), time, paths, mc_mean, mean_th, model_name + "   E[X_t]   MC vs Theory", "X_t",
                        {0.2f, 0.5f, 0.85f, 0.15f});
        save(base + "_mean.pdf");
        cla();
    }
    {
        auto f = figure(true);
        f->size(1000, 560);
        draw_var(gca(), time, mc_var, var_th, model_name + "   Var[X_t]   MC vs Theory", "Var[X_t]");
        save(base + "_var.pdf");
        cla();
    }
}

inline void plot_poisson_full(const std::vector<double> &time_grid,
                              const std::vector<std::vector<std::vector<double>>> &all_paths,
                              const std::vector<double> &lambdas, const std::string &base) {
    std::vector<std::string> cols = {"blue", "red", "green"};
    std::vector<std::string> labels = {"lambda=2", "lambda=5", "lambda=15"};

    {
        auto f = figure(true);
        f->size(1100, 500);
        auto ax = gca();
        hold(ax, true);
        for (size_t li = 0; li < lambdas.size(); ++li) {
            if (all_paths[li].empty())
                continue;
            auto l = stairs(ax, time_grid, all_paths[li][0]);
            l->line_width(2.2);
            l->color(cols[li]);
        }
        ax->title("Poisson sample paths");
        ax->xlabel("t");
        ax->ylabel("N_t");
        legend(ax, labels);
        grid(ax, true);
        save(base + "_paths.pdf");
        cla();
    }

    size_t T_steps = time_grid.size();
    std::vector<std::vector<double>> mc_means(lambdas.size());
    std::vector<std::vector<double>> mc_vars(lambdas.size());
    for (size_t li = 0; li < lambdas.size(); ++li) {
        size_t M = all_paths[li].size();
        mc_means[li].assign(T_steps, 0.0);
        mc_vars[li].assign(T_steps, 0.0);
        for (size_t j = 0; j < M; ++j)
            for (size_t t = 0; t < T_steps; ++t)
                mc_means[li][t] += all_paths[li][j][t];
        for (auto &v: mc_means[li])
            v /= M;
        for (size_t j = 0; j < M; ++j)
            for (size_t t = 0; t < T_steps; ++t) {
                double d = all_paths[li][j][t] - mc_means[li][t];
                mc_vars[li][t] += d * d;
            }
        for (auto &v: mc_vars[li])
            v /= M;
    }

    {
        auto f = figure(true);
        f->size(1100, 520);
        auto ax = gca();
        hold(ax, true);
        for (size_t li = 0; li < lambdas.size(); ++li) {
            std::vector<double> theory(T_steps);
            for (size_t t = 0; t < T_steps; ++t)
                theory[t] = lambdas[li] * time_grid[t];
            auto lm = plot(ax, time_grid, mc_means[li], "--");
            lm->line_width(2.2);
            lm->color(cols[li]);
            auto lt = plot(ax, time_grid, theory, "-");
            lt->line_width(1.6);
            lt->color(cols[li]);
        }
        ax->title("E[N_t]:  MC (dashed)  vs  lambda*t (solid)");
        ax->xlabel("t");
        ax->ylabel("E[N_t]");
        legend(ax, {"MC  l=2", "th  l=2", "MC  l=5", "th  l=5", "MC  l=15", "th  l=15"});
        grid(ax, true);
        save(base + "_mean.pdf");
        cla();
    }

    {
        auto f = figure(true);
        f->size(1100, 520);
        auto ax = gca();
        hold(ax, true);
        for (size_t li = 0; li < lambdas.size(); ++li) {
            std::vector<double> theory(T_steps);
            for (size_t t = 0; t < T_steps; ++t)
                theory[t] = lambdas[li] * time_grid[t];
            auto lm = plot(ax, time_grid, mc_vars[li], "--");
            lm->line_width(2.2);
            lm->color(cols[li]);
            auto lt = plot(ax, time_grid, theory, "-");
            lt->line_width(1.6);
            lt->color(cols[li]);
        }
        ax->title("Var[N_t]:  MC (dashed)  vs  lambda*t (solid)    [E=Var for Poisson]");
        ax->xlabel("t");
        ax->ylabel("Var[N_t]");
        legend(ax, {"MC  l=2", "th  l=2", "MC  l=5", "th  l=5", "MC  l=15", "th  l=15"});
        grid(ax, true);
        save(base + "_var.pdf");
        cla();
    }
}

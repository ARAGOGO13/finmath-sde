#pragma once

#include <vector>
#include <string>
#include "../models/poisson_process.hpp"
#include "../plotting.hpp"
#include "../console_utils.hpp"

// =============================================================================
// Эксперимент: процессы Пуассона
// =============================================================================
inline void run_poisson_experiment() {
    using namespace std;

    section("Poisson processes",
            "E[N_t] = Var[N_t] = lambda * t  |  inter-arrival ~ Exp(lambda)");

    const double T = 1.0;
    const int Np = 300;
    const int Nmc_p = 3000;

    vector<double> lambdas = {2.0, 5.0, 15.0};

    auto tg = linspace(0.0, T, Np + 1);
    vector<double> tgrid(tg.begin(), tg.end());

    vector<vector<vector<double>>> all_p(lambdas.size());
    for (size_t li = 0; li < lambdas.size(); ++li) {
        all_p[li].resize(Nmc_p);
        for (int j = 0; j < Nmc_p; ++j)
            all_p[li][j] = simulate_poisson_path(tgrid, lambdas[li]);
    }

    plot_poisson_full(tgrid, all_p, lambdas, "../plots/C6_poisson");

    table_header(18);
    for (size_t li = 0; li < lambdas.size(); ++li) {
        double sn = 0, sn2 = 0;
        for (int j = 0; j < Nmc_p; ++j) {
            double n = all_p[li][j].back();
            sn += n;
            sn2 += n * n;
        }
        double mn = sn / Nmc_p;
        double vn = sn2 / Nmc_p - mn * mn;
        double lth = lambdas[li] * T;
        string L = to_string(static_cast<int>(lambdas[li]));
        table_row("E[N_T]   lambda=" + L, mn, lth, 18);
        table_row("Var[N_T] lambda=" + L, vn, lth, 18);
    }
    print_saved({"C6_poisson_paths.pdf", "C6_poisson_mean.pdf", "C6_poisson_var.pdf"});
}
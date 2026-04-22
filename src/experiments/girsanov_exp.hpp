#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include "../girsanov_experiment.hpp"
#include "../plot_girsanov.hpp"
#include "../console_utils.hpp"

// =============================================================================
// Эксперимент по проверке теоремы Гирсанова
// =============================================================================
inline void run_girsanov_experiment() {
    using namespace std;

    section("Girsanov theorem — three estimators of the same call price",
            "BS analytic  |  MC under Q (risk-neutral)  |  MC under P (weighted)");

    const double r = 0.05, sigma = 0.20, T = 1.0, K = 100.0, S0 = 100.0;
    const int N_mc = 100000;

    vector<double> mu_grid;
    for (int i = 0; i < 15; ++i)
        mu_grid.push_back(-0.10 + i * (0.50 / 14.0));

    GirsanovExperiment exp_g(r, sigma, T, K, S0, N_mc);
    auto pts = exp_g.run(mu_grid, 42);

    plot_girsanov_prices(pts, "../plots/G1_girsanov_prices.pdf");
    plot_girsanov_variance(pts, "../plots/G1_girsanov_variance.pdf",
                           T, r, sigma, K, S0);

    // Вывод таблицы результатов
    vector<string> headers = {"lambda", "BS", "Q-MC", "P-MC", "Q-err", "P-err", "Var ratio"};
    vector<int> widths = {9, 10, 10, 10, 9, 9, 10};
    print_table_header(headers, widths);

    for (const auto& p : pts) {
        ostringstream ss_lam, ss_bs, ss_q, ss_p, ss_qe, ss_pe, ss_vr;
        ss_lam << fixed << setprecision(3) << p.lambda;
        ss_bs  << p.bs_price;
        ss_q   << p.rn_price;
        ss_p   << p.p_price;
        ss_qe  << p.rn_std_err;
        ss_pe  << p.p_std_err;
        ss_vr  << p.variance_ratio;
        print_table_row({ss_lam.str(), ss_bs.str(), ss_q.str(), ss_p.str(),
                         ss_qe.str(), ss_pe.str(), ss_vr.str()}, widths);
    }

    print_saved({"G1_girsanov_prices.pdf", "G1_girsanov_variance.pdf"});
}
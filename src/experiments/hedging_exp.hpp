#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include "../hedging_experiment.hpp"
#include "../plot_hedging.hpp"
#include "../console_utils.hpp"

// =============================================================================
// Эксперимент по дискретному дельта-хеджированию короткого колла
// =============================================================================
inline void run_hedging_experiment() {
    using namespace std;

    section("Discrete delta-hedging of a short call",
            "premium received by BS(sigma_impl), rebalanced N times, P&L at T");

    HedgingConfig cfg;  // r=0.05, T=1, K=100, S0=100, sigma=0.20
    DeltaHedgingExperiment exp_h(cfg);

    // ---- A: Гистограммы P&L для разных N_rebal ----
    vector<int> N_rebal_hist = {10, 50, 250, 1000};
    vector<HedgingPnLResult> hist_results;
    vector<string> hist_labels;
    for (int Nr : N_rebal_hist) {
        auto r = exp_h.simulate(Nr, 20000, 42u + static_cast<unsigned>(Nr));
        hist_results.push_back(std::move(r));
        hist_labels.push_back("N=" + std::to_string(Nr));
    }
    plot_hedging_pnl_histograms(hist_results, hist_labels,
        "Hedging P&L distribution  (sigma_real = sigma_impl = 0.20)",
        "../plots/G2_hedging_pnl_hist.pdf");

    vector<string> headersA = {"N_rebal", "mean(P&L)", "std(P&L)", "q05", "q95"};
    vector<int> widthsA = {10, 14, 14, 10, 10};
    print_table_header(headersA, widthsA);

    for (size_t i = 0; i < hist_results.size(); ++i) {
        const auto& r = hist_results[i];
        ostringstream ss_n, ss_m, ss_s, ss_q05, ss_q95;
        ss_n << r.N_rebal;
        ss_m << r.mean_pnl;
        ss_s << r.std_pnl;
        ss_q05 << r.q05;
        ss_q95 << r.q95;
        print_table_row({ss_n.str(), ss_m.str(), ss_s.str(), ss_q05.str(), ss_q95.str()}, widthsA);
    }

    // ---- B: std(P&L) ~ dt^alpha в log-log ----
    vector<int> N_rebal_scan = {4, 8, 16, 32, 64, 128, 256, 512, 1024};
    auto scan = exp_h.scan_rebalance_frequency(N_rebal_scan, 10000, 1000u);
    plot_hedging_std_vs_dt(scan, "../plots/G2_hedging_std_vs_dt.pdf");

    cout << "\n  std(P&L) ~ dt^alpha   measured alpha = "
         << fixed << setprecision(3) << scan.slope_log_std_vs_log_dt
         << "   theory 0.5\n";

    // ---- C: Volatility mismatch ----
    vector<double> sigma_real_list = {0.10, 0.14, 0.17, 0.20, 0.23, 0.26, 0.30};
    vector<HedgingPnLResult> mismatch_results;
    for (double sr : sigma_real_list) {
        HedgingConfig cfg2 = cfg;
        cfg2.sigma_real = sr;
        cfg2.sigma_impl = 0.20;
        DeltaHedgingExperiment exp_c(cfg2);
        auto r = exp_c.simulate(250, 20000, 42u + static_cast<unsigned>(sr * 1000));
        mismatch_results.push_back(std::move(r));
    }
    plot_hedging_vol_mismatch(mismatch_results, sigma_real_list, 0.20,
                              100.0, 100.0, 1.0, 0.05,
                              "../plots/G2_hedging_vol_mismatch.pdf");

    const double d1_atm    = ((0.05 + 0.5 * 0.20 * 0.20) * 1.0) / (0.20 * std::sqrt(1.0));
    const double k_sqrt2pi = 2.5066282746310002416;
    const double phi_d1    = std::exp(-0.5 * d1_atm * d1_atm) / k_sqrt2pi;
    const double gamma_atm = phi_d1 / (100.0 * 0.20 * std::sqrt(1.0));

    vector<string> headersC = {"sigma_real", "mean(P&L)", "std(P&L)", "theory mean"};
    vector<int> widthsC = {13, 14, 14, 14};
    print_table_header(headersC, widthsC);

    for (size_t i = 0; i < mismatch_results.size(); ++i) {
        double sr = sigma_real_list[i];
        double theory_mean = 0.5 * gamma_atm * 100.0 * 100.0 * (0.20*0.20 - sr*sr) * 1.0;
        const auto& r = mismatch_results[i];
        ostringstream ss_sr, ss_m, ss_s, ss_th;
        ss_sr << sr;
        ss_m  << r.mean_pnl;
        ss_s  << r.std_pnl;
        ss_th << theory_mean;
        print_table_row({ss_sr.str(), ss_m.str(), ss_s.str(), ss_th.str()}, widthsC);
    }

    print_saved({"G2_hedging_pnl_hist.pdf",
                 "G2_hedging_std_vs_dt.pdf",
                 "G2_hedging_vol_mismatch.pdf"});
}
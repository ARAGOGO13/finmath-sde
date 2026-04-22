#pragma once

#include "../convergence_analyzer.hpp"
#include "../models/geometric_brownian_motion.hpp"
#include "../plotting.hpp"
#include "../console_utils.hpp"

// =============================================================================
// Эксперимент: сильная сходимость схем Эйлера и Мильштейна для GBM
// =============================================================================
inline void run_convergence_experiment() {
    using namespace std;
    section("Euler vs Milshtein", "strong error e(dt) = E[|X_exact - X_scheme|]");

    const double T = 1.0;
    GeometricBrownianMotion gbm(0.05, 0.20);
    ConvergenceAnalyzer analyzer(3000, {8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096});
    auto conv = analyzer.analyze_gbm(gbm.mu(), gbm.sigma(), 100.0, T);
    plot_convergence(conv, "../plots/C3_convergence.pdf");

    cout << "\n  log-log slopes\n";
    cout << "  " << string(44, '-') << "\n";
    auto srow = [&](const string &name, double slope, double th) {
        cout << "  " << pad(name, 14)
             << fixed << setprecision(3) << setw(8) << slope
             << "   theory " << th << "\n";
    };
    srow("Euler",     conv.euler_slope,     0.5);
    srow("Milshtein", conv.milshtein_slope, 1.0);
    print_saved({"C3_convergence.pdf"});
}
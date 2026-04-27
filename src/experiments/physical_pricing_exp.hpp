#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "console_utils.hpp"
#include "physical_pricing_experiment.hpp"
#include "plot_physical_pricing.hpp"

// =============================================================================
// Эксперимент G1: ценообразование колла под физической мерой P vs риск-нейтральной Q.
//
// Строятся следующие графики:
//   G1a — зависимость C0_P(μ) и C0_Q от физического дрейфа μ;
//   G1b — распределение дисконтированного выигрыша для ITM-траекторий (4 значения μ);
//   G1c — распределение P&L покупателя = payoff - C0_P  (те же 4 значения μ);
//   G1d — структура исходов (OTM / ITM-убыток / прибыль) и ожидаемый P&L во всём диапазоне μ.
//
// После построения отдельные PDF объединяются скриптом merge_plots.sh в итоговый отчёт.
// =============================================================================

namespace g1_detail {

// Функция стандартного нормального распределения Φ(x)
static inline double Phi(double x) {
    return 0.5 * std::erfc(-x / 1.4142135623730950488);
}

// Аналитическая цена колла под физической мерой P с дрейфом μ
static inline double c0p(double mu, double r, double sigma,
                         double T, double K, double S0) {
    double sqT   = std::sqrt(T);
    double d1    = (std::log(S0 / K) + (mu + 0.5 * sigma * sigma) * T)
                   / (sigma * sqT);
    double d2    = d1 - sigma * sqT;
    return std::exp(-r * T)
           * (S0 * std::exp(mu * T) * Phi(d1) - K * Phi(d2));
}

// Риск-нейтральная цена колла (BS)
static inline double c0q(double r, double sigma, double T,
                         double K, double S0) {
    double sqT = std::sqrt(T);
    double d1  = (std::log(S0 / K) + (r + 0.5 * sigma * sigma) * T)
                 / (sigma * sqT);
    double d2  = d1 - sigma * sqT;
    return S0 * Phi(d1) - K * std::exp(-r * T) * Phi(d2);
}

// Вероятностная структура P&L покупателя при заданном μ
struct PnLProbs {
    double p_profit;      // P(P&L > 0) — безусловная вероятность прибыли
    double p_loss;        // P(P&L < 0) — безусловная вероятность убытка (включая OTM)
    double p_otm;         // P(OTM) — опцион не исполнился (убыток всей премии)

    double p_profit_itm;  // P(прибыль | ITM)
    double p_loss_itm;    // P(убыток   | ITM)

    double expected_pnl;  // Ожидаемый P&L = C0_P - C0_Q
};

// Аналитический расчёт вероятностей P&L для заданного μ
static inline PnLProbs compute_probs(double mu, double r, double sigma,
                                      double T, double K, double S0) {
    const double sqT  = std::sqrt(T);
    const double cp   = c0p(mu, r, sigma, T, K, S0);
    const double cq   = c0q(r, sigma, T, K, S0);

    // Порог безубыточности для цены актива: e^{-rT}*(S_T - K) = C0_P
    const double K_eff = K + cp * std::exp(r * T);

    const double drift_adj = (mu - 0.5 * sigma * sigma) * T;
    const double zK  = (std::log(K     / S0) - drift_adj) / (sigma * sqT);
    const double zp  = (std::log(K_eff / S0) - drift_adj) / (sigma * sqT);

    const double p_itm    = 1.0 - Phi(zK);   // P(S_T > K)
    const double p_profit = 1.0 - Phi(zp);   // P(S_T > K_eff)
    const double p_otm    = Phi(zK);          // P(S_T <= K)
    const double p_loss   = 1.0 - p_profit;  // P(P&L ≤ 0) — включает OTM

    const double p_profit_itm = (p_itm > 1e-12) ? p_profit / p_itm : 0.0;
    const double p_loss_itm   = 1.0 - p_profit_itm;

    return PnLProbs{
        p_profit, p_loss, p_otm,
        p_profit_itm, p_loss_itm,
        cp - cq
    };
}

// Форматирование вероятности в проценты
static inline std::string pct(double v, int prec = 1) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(prec) << v * 100.0 << "%";
    return ss.str();
}

// Форматирование числа со знаком
static inline std::string sgn_fmt(double v, int prec = 2) {
    std::ostringstream ss;
    ss << (v >= 0 ? "+" : "") << std::fixed << std::setprecision(prec) << v;
    return ss.str();
}

} // namespace g1_detail

// =============================================================================
// Запуск полного эксперимента G1: цены, вероятности, распределения
// =============================================================================
inline void run_physical_pricing_experiment() {
    using namespace std;
    using namespace g1_detail;

    section("Physical measure pricing vs Risk-Neutral",
            "Analytic: C0_P(mu) = e^{-rT}*[S0*e^{mu*T}*Phi(d1) - K*Phi(d2)]");

    const double r = 0.05, sigma = 0.20, T = 1.0;
    const double K = 100.0, S0 = 100.0;
    const int    N_mc = 200000;

    // Сетка значений физического дрейфа μ
    vector<double> mu_grid;
    const double mu_min = -0.10;
    const double mu_max =  0.25;
    const int    n_mu   = 20;
    for (int i = 0; i < n_mu; ++i)
        mu_grid.push_back(mu_min + i * (mu_max - mu_min) / (n_mu - 1));

    PhysicalPricingExperiment exp_p(r, sigma, T, K, S0, N_mc);
    auto pts = exp_p.run(mu_grid, 42u);

    // ── Таблица 1: аналитические и MC-цены ──
    vector<string> h1 = {
        "mu", "lambda", "RN price",
        "Phys analytic", "MC check", "Diff", "Diff%"
    };
    vector<int> w1 = {8, 9, 12, 16, 12, 10, 8};
    print_table_header(h1, w1);

    for (const auto& p : pts) {
        ostringstream ss_mu, ss_lam, ss_rn, ss_an, ss_mc, ss_diff, ss_dp;
        ss_mu << fixed << setprecision(3) << p.mu;
        ss_lam << fixed << setprecision(3) << p.lambda;
        ss_rn << fixed << setprecision(2) << p.price_rn;
        ss_an << fixed << setprecision(2) << p.price_physical_analytic;
        ss_mc << fixed << setprecision(2) << p.price_physical_mc;
        ss_diff << fixed << setprecision(4) << p.difference;
        ss_dp << fixed << setprecision(2) << p.difference_pct;
        print_table_row_with_error({ss_mu.str(), ss_lam.str(), ss_rn.str(), ss_an.str(), ss_mc.str(), ss_diff.str(), ss_dp.str()},
                                   w1, std::abs(p.price_physical_mc - p.price_physical_analytic));
    }

    // ── Таблица 2: вероятности P&L ──
    vector<string> h2 = {
        "mu", "lambda", "C0_P", "C0_Q", "P(P&L>0)", "P(P&L<0)", "P(OTM)",
        "P(profit|ITM)", "P(loss|ITM)", "E[P&L]"
    };
    vector<int> w2 = {8, 8, 9, 9, 10, 10, 9, 13, 13, 10};
    print_table_header(h2, w2);

    for (const auto& p : pts) {
        auto pb = compute_probs(p.mu, r, sigma, T, K, S0);
        const double err_val = std::abs(pb.expected_pnl);

        ostringstream mu_s, lam_s, cp_s, cq_s, pp_s, pl_s, po_s, ppi_s, pli_s, ep_s;
        mu_s << fixed << setprecision(3) << p.mu;
        lam_s << fixed << setprecision(3) << p.lambda;
        cp_s << fixed << setprecision(2) << p.price_physical_analytic;
        cq_s << fixed << setprecision(2) << p.price_rn;
        pp_s << pct(pb.p_profit, 1);
        pl_s << pct(pb.p_loss, 1);
        po_s << pct(pb.p_otm, 1);
        ppi_s << pct(pb.p_profit_itm, 1);
        pli_s << pct(pb.p_loss_itm, 1);
        ep_s << sgn_fmt(pb.expected_pnl, 2);

        print_table_row_with_error(
            {mu_s.str(), lam_s.str(), cp_s.str(), cq_s.str(), pp_s.str(), pl_s.str(), po_s.str(),
             ppi_s.str(), pli_s.str(), ep_s.str()},
            w2, err_val);
    }

    // ── Характерные значения μ для текстовой верификации ──
    const vector<double> mu_hist = {-0.05, 0.05, 0.15, 0.22};
    cout << "\n";
    for (double mu_i : mu_hist) {
        auto pb = compute_probs(mu_i, r, sigma, T, K, S0);
        double cp = c0p(mu_i, r, sigma, T, K, S0);
        double cq = c0q(r, sigma, T, K, S0);
        double K_eff = K + cp * std::exp(r * T);

        cout << "mu = " << fixed << setprecision(2) << mu_i
             << "  |  lambda = " << fixed << setprecision(3) << (mu_i - r) / sigma
             << "  |  C0_P = " << fixed << setprecision(2) << cp
             << "  |  C0_Q = " << fixed << setprecision(2) << cq << "\n";
        cout << "K_eff = " << fixed << setprecision(2) << K_eff << "\n";
        cout << "P(P&L > 0) = " << pct(pb.p_profit, 2) << "\n";
        cout << "P(P&L < 0) = " << pct(pb.p_loss, 2) << "\n";
        cout << "P(P&L > 0 | ITM) = " << pct(pb.p_profit_itm, 2) << "\n";
        cout << "P(P&L < 0 | ITM) = " << pct(pb.p_loss_itm, 2) << "\n";
        cout << "E[P&L] = " << sgn_fmt(pb.expected_pnl, 4) << "\n\n";
    }

    // ── График 1 ──
    plot_physical_vs_rn_prices(pts, "../plots/G1_physical_vs_rn_prices.pdf", r);

    // ── График 2 ──
    std::vector<PnLProbPoint> prob_pts;
    prob_pts.reserve(mu_grid.size());
    for (double mu_i : mu_grid) {
        auto pb = compute_probs(mu_i, r, sigma, T, K, S0);
        double cp = c0p(mu_i, r, sigma, T, K, S0);
        double cq = c0q(r, sigma, T, K, S0);
        prob_pts.push_back({mu_i, (mu_i - r) / sigma, cp, cq, pb.p_profit, pb.p_loss, pb.p_otm, pb.p_profit_itm, pb.p_loss_itm, pb.expected_pnl});
    }
    plot_pnl_stacked(prob_pts, "../plots/G1_pnl_stacked.pdf", r);

    // ── График 3 ──
    plot_expected_pnl(prob_pts, "../plots/G1_expected_pnl.pdf", r);

    // ── Графики распределений ──
    const int N_dist = 500000;
    const vector<double> mu_hist_vals = {-0.05, 0.05, 0.15, 0.22};
    vector<string> payoff_files;
    for (size_t i = 0; i < mu_hist_vals.size(); ++i) {
        string fn = "../plots/G1_payoff_mu_" + to_string(i) + ".pdf";
        plot_single_payoff_distribution(r, sigma, T, K, S0, mu_hist_vals[i], N_dist, fn, 100);
        payoff_files.push_back(fn);
    }

    vector<string> pnl_files;
    for (size_t i = 0; i < mu_hist_vals.size(); ++i) {
        string fn = "../plots/G1_pnl_mu_" + to_string(i) + ".pdf";
        plot_single_pnl_distribution(r, sigma, T, K, S0, mu_hist_vals[i], N_dist, fn, 100);
        pnl_files.push_back(fn);
    }

    vector<string> saved = {"G1_physical_vs_rn_prices.pdf", "G1_pnl_stacked.pdf", "G1_expected_pnl.pdf"};
    for (size_t i = 0; i < mu_hist_vals.size(); ++i) saved.push_back("G1_payoff_mu_" + to_string(i) + ".pdf");
    for (size_t i = 0; i < mu_hist_vals.size(); ++i) saved.push_back("G1_pnl_mu_" + to_string(i) + ".pdf");
    saved.push_back("G1_payoff_2x2.pdf");
    saved.push_back("G1_pnl_2x2.pdf");
    saved.push_back("G1_physical_pricing_merged.pdf");
    print_saved(saved);
}

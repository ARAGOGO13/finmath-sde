#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <matplot/matplot.h>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "physical_pricing_experiment.hpp"

using namespace matplot;

// =============================================================================
// Вспомогательные функции
// =============================================================================

// Функция стандартного нормального распределения Φ(x)
static inline double phi_cdf(double x) {
    return 0.5 * std::erfc(-x / 1.4142135623730950488);
}

// Физическая цена колла C0_P(μ)
static inline double c0_physical(double mu, double r, double sigma, double T, double K, double S0) {
    double sqT = std::sqrt(T);
    double d1_mu = (std::log(S0 / K) + (mu + 0.5 * sigma * sigma) * T) / (sigma * sqT);
    double d2_mu = d1_mu - sigma * sqT;
    return std::exp(-r * T) * (S0 * std::exp(mu * T) * phi_cdf(d1_mu) - K * phi_cdf(d2_mu));
}

// Риск-нейтральная цена колла (Блэк–Шоулз)
static inline double c0_rn(double r, double sigma, double T, double K, double S0) {
    double sqT = std::sqrt(T);
    double d1 = (std::log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqT);
    double d2 = d1 - sigma * sqT;
    return S0 * phi_cdf(d1) - K * std::exp(-r * T) * phi_cdf(d2);
}

// Форматирование числа с фиксированной точностью
static inline std::string fmt(double v, int prec = 2) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(prec) << v;
    return ss.str();
}

// Цвета распределений для заданных значений μ:
//   -0.05 → голубой,  0.05 → серый,  0.15 → оранжевый,  0.22 → красный
static inline std::string distribution_color(double mu) {
    if (std::abs(mu + 0.05) < 1e-6) return "#17becf";   // cyan
    if (std::abs(mu - 0.05) < 1e-6) return "#7f7f7f";   // grey
    if (std::abs(mu - 0.15) < 1e-6) return "#ff7f0e";   // orange
    if (std::abs(mu - 0.22) < 1e-6) return "#d62728";   // red
    return "#4D4D4D";                                    // fallback dark grey
}

// =============================================================================
// Нормализация гистограммы на плотность вероятности (площадь = 1)
// =============================================================================
static inline std::vector<double> density_normalize(const std::vector<int> &counts, double total, double bin_width) {
    std::vector<double> h(counts.size());
    for (size_t i = 0; i < counts.size(); ++i)
        h[i] = static_cast<double>(counts[i]) / (total * bin_width);
    return h;
}

// =============================================================================
// Разметка горизонтальной оси с красивым шагом и принудительными метками
// =============================================================================
static inline void apply_x_ticks(axes_handle ax, double x_lo, double x_hi, bool force_zero,
                                 double special_val, const std::string &special_label) {
    const double x_range = x_hi - x_lo;
    if (x_range < 1e-12) return;

    const double raw = x_range / 10.0;
    double step = [&]() -> double {
        double mag = std::pow(10.0, std::floor(std::log10(raw)));
        double norm = raw / mag;
        if (norm < 1.5) return 1.0 * mag;
        if (norm < 2.25) return 2.0 * mag;
        if (norm < 3.75) return 2.5 * mag;
        if (norm < 7.5) return 5.0 * mag;
        return 10.0 * mag;
    }();

    double first = std::ceil(x_lo / step) * step;
    std::vector<double> pos;
    std::vector<std::string> lbl;

    for (double t = first; t <= x_hi + step * 0.01; t += step) {
        double tr = std::round(t / step) * step;
        pos.push_back(tr);
        std::ostringstream s;
        if (step >= 1.0)      s << std::fixed << std::setprecision(0) << tr;
        else if (step >= 0.1) s << std::fixed << std::setprecision(1) << tr;
        else                  s << std::fixed << std::setprecision(2) << tr;
        lbl.push_back(s.str());
    }

    // Принудительно x=0
    if (force_zero && x_lo < 0.0 && x_hi > 0.0) {
        bool found = false;
        for (double p : pos) if (std::abs(p) < step * 0.05) { found = true; break; }
        if (!found) { pos.push_back(0.0); lbl.push_back("0"); }
    }

    // Принудительно special_val
    if (!special_label.empty() && special_val >= x_lo && special_val <= x_hi) {
        bool found = false;
        for (double p : pos) if (std::abs(p - special_val) < step * 0.05) { found = true; break; }
        if (!found) { pos.push_back(special_val); lbl.push_back(special_label); }
    }

    std::vector<std::pair<double, std::string>> pairs;
    pairs.reserve(pos.size());
    for (size_t i = 0; i < pos.size(); ++i) pairs.push_back({pos[i], lbl[i]});
    std::sort(pairs.begin(), pairs.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });

    pos.clear(); lbl.clear();
    for (const auto &[p, l] : pairs) { pos.push_back(p); lbl.push_back(l); }
    ax->x_axis().tick_values(pos);
    ax->x_axis().ticklabels(lbl);
}

// Построение гистограммы плотности
struct HistData {
    std::vector<double> centers;
    std::vector<double> heights;
    double x_lo, x_hi, y_top;
};

static inline HistData make_histogram(const std::vector<double> &data, int num_bins) {
    const double x_min = *std::min_element(data.begin(), data.end());
    const double x_max = *std::max_element(data.begin(), data.end());
    const double margin = 0.08 * (x_max - x_min + 1e-6);
    const double x_lo = x_min - margin;
    const double x_hi = x_max + margin;
    const double bw   = (x_hi - x_lo) / num_bins;

    std::vector<double> centers(num_bins);
    std::vector<int>    counts(num_bins, 0);
    for (int i = 0; i < num_bins; ++i) centers[i] = x_lo + (i + 0.5) * bw;
    for (double v : data) {
        int b = static_cast<int>((v - x_lo) / bw);
        b = std::max(0, std::min(num_bins - 1, b));
        counts[b]++;
    }
    auto heights = density_normalize(counts, data.size(), bw);
    const double y_top = *std::max_element(heights.begin(), heights.end()) * 1.25;
    return {centers, heights, x_lo, x_hi, y_top};
}

static inline std::array<float,4> hex_to_rgba(const std::string &hex, float alpha) {
    unsigned int r=0,g=0,b=0;
    std::sscanf(hex.c_str()+1, "%02x%02x%02x", &r,&g,&b);
    return {b/255.f, g/255.f, r/255.f, alpha};   // порядок BGR — так нужно Matplot++
}

static inline void filled_band(axes_handle ax_h,
                               const std::vector<double> &x,
                               const std::vector<double> &y_lo,
                               const std::vector<double> &y_hi,
                               const std::string &hex_color, float alpha = 0.75f) {
    const int m = x.size();
    if (m < 2) return;
    std::vector<double> px, py;
    px.reserve(2*m+2); py.reserve(2*m+2);
    for (int i=0; i<m; ++i) { px.push_back(x[i]); py.push_back(y_hi[i]); }
    for (int i=m-1; i>=0; --i) { px.push_back(x[i]); py.push_back(y_lo[i]); }
    px.push_back(px[0]); py.push_back(py[0]);
    auto ar = fill(ax_h, px, py);
    ar->color(hex_to_rgba(hex_color, alpha));
}

// Структура вероятностей P&L для одного значения μ
struct PnLProbPoint {
    double mu, lambda, c0_p, c0_q;
    double p_profit, p_itm_loss, p_otm;
    double p_profit_itm, p_loss_itm;
    double expected_pnl;
};

// =============================================================================
// G1a – Цена колла C0 в зависимости от физического дрейфа μ
// =============================================================================
static inline void plot_physical_vs_rn_prices(
        const std::vector<PhysicalPricingExperiment::Point> &pts,
        const std::string &filename, double r = 0.05) {
    const int n = pts.size();
    std::vector<double> mu_all(n), p_all(n);
    for (int i=0; i<n; ++i) {
        mu_all[i] = pts[i].mu;
        p_all[i]  = pts[i].price_physical_analytic;
    }
    const double price_rn = pts[0].price_rn;
    const double y_top    = *std::max_element(p_all.begin(), p_all.end()) * 1.07;

    auto f = figure(true); f->size(2200, 1250);
    auto ax = gca(); hold(ax, true);

    auto lp = plot(ax, mu_all, p_all, "-");
    lp->line_width(3.5); lp->color("#0A8C0A");

    auto lrn = plot(ax, mu_all, std::vector<double>(n, price_rn), "--");
    lrn->line_width(2.5); lrn->color("#3A3A3A");

    auto lv = plot(ax, {r, r}, {0.0, y_top}, "--");
    lv->line_width(2.0); lv->color({0.5f,0.5f,0.5f,1.f});

    ylim(ax, {0.0, y_top});
    auto leg = legend(ax, {
        "C0_P(μ)  (physical measure)",
        "C0_Q = " + fmt(price_rn) + "  (risk‑neutral BS)",
        "μ = r = " + fmt(r,3) + "  (P = Q, prices coincide)"
    });
    leg->inside(false);
    leg->location(legend::general_alignment::topright);
    leg->font_size(8);
    leg->num_columns(1);

    ax->title("Call price C0 vs physical drift μ\nC0_Q (dashed) is independent of μ — Girsanov theorem");
    ax->xlabel("μ  (physical drift of GBM under P)");
    ax->ylabel("Call price  C0");
    grid(ax, true);
    apply_x_ticks(ax, mu_all.front(), mu_all.back(), false, r, "r=" + fmt(r,3));
    save(filename); cla();
}

// =============================================================================
// G1b – Распределение дисконтированного выигрыша (только ITM-траектории)
//       X = e^{-rT}(S_T – K)
//       Вертикали: C0_P (цена покупателя), C0_Q (рыночная цена BS)
// =============================================================================
static inline std::string plot_single_payoff_distribution(
        double r, double sigma, double T, double K, double S0, double mu,
        int N_samples, const std::string &filename, int num_bins = 100) {
    const double discount = std::exp(-r * T);
    const double sqrt_T   = std::sqrt(T);
    const double c0_p     = c0_physical(mu, r, sigma, T, K, S0);
    const double bs_price = c0_rn(r, sigma, T, K, S0);
    const std::string color = distribution_color(mu);

    std::vector<double> payoffs;
    payoffs.reserve(N_samples/3);

    std::mt19937_64 rng(7777u + static_cast<unsigned>(mu*1000) + 12345u);
    std::normal_distribution<double> nd(0.0, 1.0);
    const double drift = (mu - 0.5*sigma*sigma)*T;
    const double vol   = sigma * sqrt_T;

    for (int i=0; i<N_samples; ++i) {
        double S_T = S0 * std::exp(drift + vol*nd(rng));
        if (S_T <= K) continue;
        payoffs.push_back(discount * (S_T - K));
    }
    if (payoffs.empty()) return "";

    auto hd = make_histogram(payoffs, num_bins);

    auto f = figure(true); f->size(2200, 1250);
    auto ax = gca(); hold(ax, true);

    auto bar_h = bar(ax, hd.centers, hd.heights);
    bar_h->face_color(color);
    bar_h->edge_color("none");

    auto lp = plot(ax, {c0_p, c0_p}, {0.0, hd.y_top}, "-");
    lp->line_width(2.0); lp->color({0.4f,0.4f,0.4f,1.f});

    auto lq = plot(ax, {bs_price, bs_price}, {0.0, hd.y_top}, "--");
    lq->line_width(2.0); lq->color("black");

    xlim(ax, {hd.x_lo, hd.x_hi});
    ylim(ax, {0.0, hd.y_top});

    ax->title("Payoff distribution (S_T > K)  for μ = " + fmt(mu,2) +
              "\nDiscounted payoff conditional on ITM");
    ax->xlabel("Discounted payoff  e^{-rT}(S_T – K)");
    ax->ylabel("Probability density  (area = 1)");

    // Легенда компактная, справа вверху (аналогично P&L‑графикам)
    auto leg = legend(ax, {
        "μ = " + fmt(mu,2),
        "C0_P = " + fmt(c0_p,2),
        "C0_Q = " + fmt(bs_price,2)
    });
    leg->inside(false);
    leg->location(legend::general_alignment::topright);      // <-- изменено на topright
    leg->font_size(8);
    leg->num_columns(1);

    grid(ax, true);
    apply_x_ticks(ax, hd.x_lo, hd.x_hi, false, 0.0, "");
    save(filename); cla();
    return filename;
}

// =============================================================================
// G1c – Распределение P&L покупателя (только ITM-траектории)
//       X = e^{-rT}(S_T–K) – C0_P
//       Вертикали: P&L = 0 (безубыточность), E[P&L] = C0_P – C0_Q
// =============================================================================
static inline std::string plot_single_pnl_distribution(
        double r, double sigma, double T, double K, double S0, double mu,
        int N_samples, const std::string &filename, int num_bins = 100) {
    const double discount = std::exp(-r * T);
    const double sqrt_T   = std::sqrt(T);
    const double c0_p     = c0_physical(mu, r, sigma, T, K, S0);
    const double bs_price = c0_rn(r, sigma, T, K, S0);
    const double expected_pnl = c0_p - bs_price;
    const std::string color = distribution_color(mu);

    std::vector<double> pnl;
    pnl.reserve(N_samples/3);

    std::mt19937_64 rng(8888u + static_cast<unsigned>(mu*1000) + 54321u);
    std::normal_distribution<double> nd(0.0, 1.0);
    const double drift = (mu - 0.5*sigma*sigma)*T;
    const double vol   = sigma * sqrt_T;

    for (int i=0; i<N_samples; ++i) {
        double S_T = S0 * std::exp(drift + vol*nd(rng));
        if (S_T <= K) continue;
        pnl.push_back(discount*(S_T - K) - c0_p);
    }
    if (pnl.empty()) return "";

    auto hd = make_histogram(pnl, num_bins);

    auto f = figure(true); f->size(2200, 1250);
    auto ax = gca(); hold(ax, true);

    auto bar_h = bar(ax, hd.centers, hd.heights);
    bar_h->face_color(color);
    bar_h->edge_color("none");

    auto l0 = plot(ax, {0.0, 0.0}, {0.0, hd.y_top}, "--");
    l0->line_width(2.0); l0->color("black");

    auto le = plot(ax, {expected_pnl, expected_pnl}, {0.0, hd.y_top}, "-");
    le->line_width(2.0); le->color({0.4f,0.4f,0.4f,1.f});

    xlim(ax, {hd.x_lo, hd.x_hi});
    ylim(ax, {0.0, hd.y_top});

    ax->title("P&L distribution (S_T > K)  for μ = " + fmt(mu,2) +
              "\nP&L = e^{-rT}(S_T – K) – C0_P(μ)");
    ax->xlabel("P&L = e^{-rT}(S_T – K) – C0_P(μ)");
    ax->ylabel("Probability density  (area = 1)");

    auto leg = legend(ax, {
        "μ = " + fmt(mu,2),
        "P&L = 0  (breakeven)",
        "E[P&L] = " + fmt(expected_pnl,2)
    });
    leg->inside(false);
    leg->location(legend::general_alignment::topright);
    leg->font_size(8);
    leg->num_columns(1);

    grid(ax, true);
    apply_x_ticks(ax, hd.x_lo, hd.x_hi, true, 0.0, "");
    save(filename); cla();
    return filename;
}

// =============================================================================
// G1d – Накопительная диаграмма: структура исходов P(OTM) + P(ITM убыток) + P(прибыль) = 100%
//
//   P(OTM)    — акция не выросла выше страйка K, опцион истёк worthless (серый).
//   P(ITM убыток) — акция выросла выше K, но payoff не окупил цену (красный).
//   P(прибыль) — акция выросла выше порога безубыточности K_eff (зелёный).
// =============================================================================
static inline void plot_pnl_stacked(const std::vector<PnLProbPoint> &pts,
                                    const std::string &filename, double r = 0.05) {
    const int n = pts.size();
    std::vector<double> mu_v(n), p_otm(n), bot_itm(n), bot_profit(n);
    for (int i=0; i<n; ++i) {
        mu_v[i]       = pts[i].mu;
        p_otm[i]      = pts[i].p_otm * 100.0;
        bot_itm[i]    = pts[i].p_otm * 100.0;
        bot_profit[i] = (pts[i].p_otm + pts[i].p_itm_loss) * 100.0;
    }
    std::vector<double> zeros(n,0.0), hundreds(n,100.0);

    auto f = figure(true); f->size(2200, 1250);
    auto ax = gca(); hold(ax, true);

    // Слой 1: P(OTM) — серый
    filled_band(ax, mu_v, zeros, p_otm, "#7F8C8D", 0.82f);
    // Слой 2: P(ITM убыток) — красный
    filled_band(ax, mu_v, bot_itm, bot_profit, "#E74C3C", 0.75f);
    // Слой 3: P(прибыль) — зелёный
    filled_band(ax, mu_v, bot_profit, hundreds, "#27AE60", 0.75f);

    auto lb1 = plot(ax, mu_v, p_otm, "-");
    lb1->line_width(1.2); lb1->color("#5D6D7E");
    auto lb2 = plot(ax, mu_v, bot_profit, "-");
    lb2->line_width(1.2); lb2->color("#A93226");

    auto lh = plot(ax, mu_v, std::vector<double>(n,50.0), ":");
    lh->line_width(1.6); lh->color({0.25f,0.25f,0.25f,0.55f});

    auto lv = plot(ax, {r,r}, {0.0,100.0}, "--");
    lv->line_width(2.2); lv->color("black");

    ylim(ax, {0.0, 100.0});
    // Легенда в левом нижнем углу, чуть мельче
    auto leg = legend(ax, {
        "P(OTM): S_T < K",
        "P(ITM loss): K < S_T < K_eff",
        "P(profit): S_T > K_eff",
        "50% reference",
        "μ = r = " + fmt(r,3) + "  [C0_P = C0_Q]"
    });
    leg->inside(false);
    leg->location(legend::general_alignment::bottomleft);   // левый нижний
    leg->font_size(7);                                      // меньше, 7 pt
    leg->num_columns(1);

    ax->title("Buyer outcome structure vs physical drift μ\n"
              "OTM + ITM-loss + profit = 100%");
    ax->xlabel("μ  (physical drift of GBM under P)");
    ax->ylabel("Probability  (%)");
    grid(ax, true);
    apply_x_ticks(ax, mu_v.front(), mu_v.back(), false, r, "r=" + fmt(r,3));
    save(filename); cla();
}

// =============================================================================
// G1d2 – Ожидаемый P&L покупателя = C0_P – C0_Q в зависимости от μ
// =============================================================================
static inline void plot_expected_pnl(const std::vector<PnLProbPoint> &pts,
                                     const std::string &filename, double r = 0.05) {
    const int n = pts.size();
    std::vector<double> mu_v(n), epnl(n);
    for (int i=0; i<n; ++i) {
        mu_v[i] = pts[i].mu;
        epnl[i] = pts[i].expected_pnl;
    }

    std::vector<double> mu_pos, ep_pos, mu_neg, ep_neg;
    for (int i=0; i<n; ++i) {
        if (epnl[i] >= 0.0) {
            mu_pos.push_back(mu_v[i]); ep_pos.push_back(epnl[i]);
        } else {
            mu_neg.push_back(mu_v[i]); ep_neg.push_back(epnl[i]);
        }
    }

    auto f = figure(true); f->size(2200, 1250);
    auto ax = gca(); hold(ax, true);

    filled_band(ax, mu_pos, std::vector<double>(mu_pos.size(),0.0), ep_pos,
                "#27AE60", 0.28f);
    filled_band(ax, mu_neg, ep_neg, std::vector<double>(mu_neg.size(),0.0),
                "#C0392B", 0.28f);

    auto l0 = plot(ax, mu_v, std::vector<double>(n,0.0), "--");
    l0->line_width(1.8); l0->color("black");

    auto le = plot(ax, mu_v, epnl, "o-");
    le->line_width(2.8); le->color("#D07000"); le->marker_size(7);

    double y_abs = 0.0;
    for (double v : epnl) y_abs = std::max(y_abs, std::abs(v));
    auto lv = plot(ax, {r,r}, {-y_abs*1.15, y_abs*1.15}, "--");
    lv->line_width(2.2); lv->color("#888888");

    ylim(ax, {-y_abs*1.18, y_abs*1.18});
    auto leg = legend(ax, {
        "E[P&L] = 0  (fair price)",
        "E[P&L] = C0_P – C0_Q",
        "μ = r = " + fmt(r,3)
    });
    leg->inside(false);
    leg->location(legend::general_alignment::topright);
    leg->font_size(8);
    leg->num_columns(1);

    ax->title("Expected P&L of buyer = C0_P – C0_Q  vs physical drift μ\n"
              "green: buyer gains,  red: buyer loses");
    ax->xlabel("μ  (physical drift of GBM under P)");
    ax->ylabel("E[P&L]  =  C0_P – C0_Q");
    grid(ax, true);
    apply_x_ticks(ax, mu_v.front(), mu_v.back(), false, r, "r=" + fmt(r,3));
    save(filename); cla();
}
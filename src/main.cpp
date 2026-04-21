#include "convergence_analyzer.hpp"
#include "models/cir.hpp"
#include "models/geometric_brownian_motion.hpp"
#include "models/ornstein_uhlenbeck.hpp"
#include "models/poisson_process.hpp"
#include "monte_carlo_engine.hpp"
#include "plotting.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <matplot/matplot.h>
#include <random>
#include <string>
#include <vector>


using namespace std;
using namespace matplot;

namespace A {
    const string R = "\033[0m";
    const string B = "\033[1m";
    const string D = "\033[2m";
    const string CY = "\033[36m";
    const string GR = "\033[32m";
    const string YL = "\033[33m";
} // namespace A
static size_t utf8_len(const string &s) {
    size_t len = 0;
    for (unsigned char c: s)
        if ((c & 0xC0) != 0x80)
            ++len;
    return len;
}
static string pad(const string &s, int w) {
    int sp = w - (int) utf8_len(s);
    return s + string(std::max(sp, 1), ' ');
}
void section(const string &title, const string &sub = "") {
    cout << "\n"
         << A::B << A::CY << "╔══════════════════════════════════════════════════════════════╗\n"
         << "║  " << pad(title, 60) << "║\n";
    if (!sub.empty())
        cout << "║  " << A::D << pad(sub, 60) << A::R << A::CY << "║\n";
    cout << "╚══════════════════════════════════════════════════════════════╝" << A::R << "\n";
}
void table_header(int lw = 16) {
    cout << "  " << A::B << pad("", lw + 2) << setw(12) << "MC"
         << "  " << setw(12) << "Theory"
         << "  " << setw(8) << "Err %" << A::R << "\n"
         << "  " << string(lw + 40, '-') << "\n";
}
void table_row(const string &label, double mc, double theory, int lw = 16) {
    double err = std::abs(mc - theory) / (std::abs(theory) + 1e-15) * 100.0;
    string flag = err < 2.0 ? A::GR + "✓" + A::R : A::YL + "~" + A::R;
    cout << "  " << A::D << pad(label, lw) << A::R << "  " << A::GR << fixed << setprecision(6) << setw(12) << mc
         << A::R << "  " << A::D << setw(12) << theory << A::R << "  " << A::YL << setprecision(3) << setw(7) << err
         << "%" << A::R << "  " << flag << "\n";
}

void print_saved(const std::vector<string> &files) {
    for (const auto &f: files)
        cout << "  " << A::GR << "✓" << A::R << A::D << "  " << f << A::R << "\n";
}

vector<double> simulate_poisson_path(const vector<double> &tg, double lam) {
    static thread_local mt19937 gen(random_device{}());
    exponential_distribution<double> exp_d(lam);
    vector<double> path(tg.size(), 0.0);
    double nxt = exp_d(gen);
    int cnt = 0;
    for (size_t i = 0; i < tg.size(); ++i) {
        while (nxt <= tg[i]) {
            ++cnt;
            nxt += exp_d(gen);
        }
        path[i] = cnt;
    }
    return path;
}

double cir_exact_var(double x0, double t, double theta, double mu, double sigma) {
    double e = std::exp(-theta * t);
    double s2 = sigma * sigma;
    return x0 * (s2 / theta) * e * (1.0 - e) + (mu * s2) / (2.0 * theta) * (1.0 - e) * (1.0 - e);
}


int main() {
    // freopen("/dev/null", "w", stderr);

    std::system("mkdir -p ../plots 2>/dev/null");

    const double T = 1.0;
    const int N_steps = 10000;
    const int N_paths = 1000;
    const int N_mc = 20000;

    auto tv = linspace(0.0, T, N_steps + 1);
    vector<double> time_std(tv.begin(), tv.end());

    section("GBM  +  InverseGBM  (U=1/S)", "dS=muS dt+sS dW  |  dU=(s^2-mu)U dt-sU dW");
    {
        const double mu = 0.05, sigma = 0.20, S0 = 100.0;
        const double U0 = 1.0 / S0, mu_star = sigma * sigma - mu;
        const double dt = T / N_steps, sqdt = std::sqrt(dt);

        mt19937 rng_g(42), rng_i(99);
        normal_distribution<double> nd(0.0, 1.0);

        vector<vector<double>> gbm_p(N_mc, vector<double>(N_steps + 1));
        vector<vector<double>> inv_p(N_mc, vector<double>(N_steps + 1));

        for (int j = 0; j < N_mc; ++j) {
            gbm_p[j][0] = S0;
            for (int k = 0; k < N_steps; ++k) {
                double dW = sqdt * nd(rng_g), S = gbm_p[j][k];
                gbm_p[j][k + 1] = S + mu * S * dt + sigma * S * dW + 0.5 * sigma * S * sigma * (dW * dW - dt);
            }
        }
        for (int j = 0; j < N_mc; ++j) {
            inv_p[j][0] = U0;
            for (int k = 0; k < N_steps; ++k) {
                double dW = sqdt * nd(rng_i), U = inv_p[j][k];
                inv_p[j][k + 1] = U + mu_star * U * dt + (-sigma) * U * dW + 0.5 * sigma * sigma * U * (dW * dW - dt);
            }
        }

        vector<double> gm(N_steps + 1), gv(N_steps + 1), im(N_steps + 1), iv(N_steps + 1);
        for (int i = 0; i <= N_steps; ++i) {
            double t = time_std[i];
            gm[i] = S0 * std::exp(mu * t);
            gv[i] = S0 * S0 * std::exp(2 * mu * t) * (std::exp(sigma * sigma * t) - 1.0);
            im[i] = U0 * std::exp(mu_star * t);
            iv[i] = U0 * U0 * std::exp(2 * mu_star * t) * (std::exp(sigma * sigma * t) - 1.0);
        }

        vector<vector<double>> gp(gbm_p.begin(), gbm_p.begin() + N_paths);
        vector<vector<double>> ip(inv_p.begin(), inv_p.begin() + N_paths);

        plot_gbm_and_inverse(time_std, gp, gm, gv, ip, im, iv, "../plots/C4_gbm");

        vector<double> gmc_m, gmc_v, imc_m, imc_v;
        compute_mc_moments(gbm_p, gmc_m, gmc_v);
        compute_mc_moments(inv_p, imc_m, imc_v);

        table_header(20);
        table_row("GBM  E[S_T]", gmc_m.back(), gm.back(), 20);
        table_row("GBM  Var[S_T]", gmc_v.back(), gv.back(), 20);
        table_row("Inv  E[U_T]", imc_m.back(), im.back(), 20);
        table_row("Inv  Var[U_T]", imc_v.back(), iv.back(), 20);
        print_saved({"C4_gbm_mean.pdf", "C4_gbm_var.pdf", "C4_gbm_inv_mean.pdf", "C4_gbm_inv_var.pdf"});
    }

    section("Euler  vs  Milshtein", "strong error  e(dt) = E[|X_exact - X_scheme|]");
    {
        GeometricBrownianMotion gbm(0.05, 0.20);
        ConvergenceAnalyzer analyzer(3000, {8, 16, 32, 64, 128, 256, 512});
        auto conv = analyzer.analyze_gbm(gbm.mu(), gbm.sigma(), 100.0, T);
        plot_convergence(conv, "../plots/C3_convergence.pdf");

        cout << "\n  " << A::B << "log-log slopes" << A::R << "\n";
        cout << "  " << string(44, '-') << "\n";
        auto srow = [&](const string &name, double slope, double th) {
            string flag = std::abs(slope - th) < 0.05 ? A::GR + "✓" + A::R : A::YL + "~" + A::R;
            cout << "  " << A::D << pad(name, 14) << A::R << A::GR << fixed << setprecision(3) << setw(8) << slope
                 << A::R << A::D << "   theory " << th << A::R << "  " << flag << "\n";
        };
        srow("Euler", conv.euler_slope, 0.5);
        srow("Milshtein", conv.milshtein_slope, 1.0);
        print_saved({"C3_convergence.pdf"});
    }

    section("Vasicek / Ornstein-Uhlenbeck", "dX=theta(mu-X)dt+sigma dW  |  sigma'=0 => Euler=Milshtein");
    {
        const double theta_v = 2.0, mu_v = 0.05, sigma_v = 0.02, x0_v = 0.03;
        OrnsteinUhlenbeck vasicek(theta_v, mu_v, sigma_v);

        mt19937 rng(123);
        normal_distribution<double> nd(0.0, 1.0);
        const double dt = T / N_steps, sqdt = std::sqrt(dt);

        vector<vector<double>> all(N_mc, vector<double>(N_steps + 1));
        for (int j = 0; j < N_mc; ++j) {
            all[j][0] = x0_v;
            for (int k = 0; k < N_steps; ++k) {
                double dW = sqdt * nd(rng), X = all[j][k];
                all[j][k + 1] = X + vasicek.drift(X, 0) * dt + vasicek.diffusion(X, 0) * dW;
            }
        }

        vector<double> mth(N_steps + 1), vth(N_steps + 1);
        for (int i = 0; i <= N_steps; ++i) {
            mth[i] = vasicek.exact_mean(x0_v, time_std[i]);
            vth[i] = vasicek.exact_var(time_std[i]);
        }

        vector<vector<double>> pp(all.begin(), all.begin() + N_paths);
        plot_mean_reversion(time_std, pp, mth, vth, "Vasicek", "../plots/C4_vasicek");

        vector<double> mc_m, mc_v;
        compute_mc_moments(all, mc_m, mc_v);
        table_header(16);
        table_row("E[X_T]", mc_m.back(), mth.back());
        table_row("Var[X_T]", mc_v.back(), vth.back());
        print_saved({"C4_vasicek_mean.pdf", "C4_vasicek_var.pdf"});
    }

    section("Cox-Ingersoll-Ross", "dX=theta(mu-X)dt+sigma*sqrt(X)dW  |  Feller: 2*theta*mu>=sigma^2");
    {
        const double theta_c = 2.0, mu_c = 0.05, sigma_c = 0.1, x0_c = 0.03;
        CoxIngersollRoss cir(theta_c, mu_c, sigma_c);
        bool feller = cir.feller_condition();
        cout << "  Feller: 2θμ=" << 2 * theta_c * mu_c << "  σ²=" << sigma_c * sigma_c << "  "
             << (feller ? A::GR + "✓ satisfied" + A::R : A::YL + "✗ violated" + A::R) << "\n\n";

        mt19937 rng(321);
        normal_distribution<double> nd(0.0, 1.0);
        const double dt = T / N_steps, sqdt = std::sqrt(dt);

        vector<vector<double>> all(N_mc, vector<double>(N_steps + 1));
        for (int j = 0; j < N_mc; ++j) {
            all[j][0] = x0_c;
            for (int k = 0; k < N_steps; ++k) {
                double dW = sqdt * nd(rng);
                double X = std::max(all[j][k], 0.0);
                double s = cir.diffusion(X, 0), sp = cir.diffusion_derivative(X, 0);
                all[j][k + 1] = X + cir.drift(X, 0) * dt + s * dW + 0.5 * s * sp * (dW * dW - dt);
            }
        }

        OrnsteinUhlenbeck ou_for_mean(theta_c, mu_c, sigma_c);
        vector<double> mth(N_steps + 1), vth(N_steps + 1);
        for (int i = 0; i <= N_steps; ++i) {
            mth[i] = ou_for_mean.exact_mean(x0_c, time_std[i]);
            vth[i] = cir_exact_var(x0_c, time_std[i], theta_c, mu_c, sigma_c);
        }

        vector<vector<double>> pp(all.begin(), all.begin() + N_paths);
        plot_mean_reversion(time_std, pp, mth, vth, "CIR", "../plots/C4_cir");

        vector<double> mc_m, mc_v;
        compute_mc_moments(all, mc_m, mc_v);
        table_header(16);
        table_row("E[X_T]", mc_m.back(), mth.back());
        table_row("Var[X_T]", mc_v.back(), vth.back());
        print_saved({"C4_cir_mean.pdf", "C4_cir_var.pdf"});
    }

    section("Poisson processes", "E[N_t]=Var[N_t]=lambda*t  |  inter-arrival ~ Exp(lambda)");
    {
        vector<double> lambdas = {2.0, 5.0, 15.0};
        const int Np = 300, Nmc_p = 3000;

        auto tg = linspace(0.0, T, Np + 1);
        vector<double> tgrid(tg.begin(), tg.end());

        vector<vector<vector<double>>> all_p(lambdas.size());
        for (size_t li = 0; li < lambdas.size(); ++li) {
            all_p[li].resize(Nmc_p);
            for (int j = 0; j < Nmc_p; ++j)
                all_p[li][j] = simulate_poisson_path(tgrid, lambdas[li]);
        }

        plot_poisson_full(tgrid, all_p, lambdas, "../plots/C6_poisson");

        table_header(20);
        for (size_t li = 0; li < lambdas.size(); ++li) {
            double sn = 0, sn2 = 0;
            for (int j = 0; j < Nmc_p; ++j) {
                double n = all_p[li][j].back();
                sn += n;
                sn2 += n * n;
            }
            double mn = sn / Nmc_p, vn = sn2 / Nmc_p - mn * mn, lth = lambdas[li] * T;
            string L = to_string((int) lambdas[li]);
            table_row("E[N_T]   λ=" + L, mn, lth, 20);
            table_row("Var[N_T] λ=" + L, vn, lth, 20);
        }
        print_saved({"C6_poisson_paths.pdf", "C6_poisson_mean.pdf", "C6_poisson_var.pdf"});
    }

    cout << "\n  " << A::D << "All plots saved to ../plots/" << A::R << "\n\n";

    return 0;
}

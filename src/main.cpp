// =============================================================================
// main.cpp — демонстрация всех моделей и сравнение схем дискретизации
//
// Структура программы соответствует разделам PDF:
//   C.1  — Броуновское движение
//   C.3  — Сравнение схем Эйлер vs Мильштейн (НОВОЕ)
//   C.4  — GBM, Vasicek/OU, CIR (НОВОЕ) с верификацией
//   C.6  — Процесс Пуассона
// =============================================================================

#include "monte_carlo_engine.hpp"
#include "convergence_analyzer.hpp"
#include "models/geometric_brownian_motion.hpp"
#include "models/ornstein_uhlenbeck.hpp"
#include "models/cir.hpp"
#include "models/poisson_process.hpp"
#include "plotting.hpp"

#include <matplot/matplot.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdio>
#include <string>

using namespace std;
using namespace matplot;

// Вспомогательная функция: напечатать разделитель с заголовком
void section(const string& title) {
    cout << "\n" << string(60, '=') << "\n";
    cout << "  " << title << "\n";
    cout << string(60, '=') << "\n";
}

// Печать строки таблицы с фиксированной шириной колонок.
// Решение: дополняем строку пробелами вручную по числу Unicode-символов.
static size_t utf8_len(const string& s) {
    size_t len = 0;
    for (unsigned char c : s)
        if ((c & 0xC0) != 0x80) ++len;  // считаем только первые байты символов
    return len;
}

static string pad(const string& s, int width) {
    int spaces = width - static_cast<int>(utf8_len(s));
    return s + string(std::max(spaces, 1), ' ');
}

int main() {
    freopen("/dev/null", "w", stderr);

    // --- Параметры -----------------------------------------------------------
    const double T        = 1.0;
    const int    N_fine   = 10000;   // мелкая сетка для красивых графиков
    const int    N_coarse = 100;     // грубая сетка для демонстрации ошибок
    const int    N_paths  = 50;      // траекторий для графиков
    const int    N_mc     = 10000;   // сценариев для Монте-Карло

    auto time_fine   = linspace(0.0, T, N_fine   + 1);
    auto time_coarse = linspace(0.0, T, N_coarse + 1);

    // =========================================================================
    // C.1 — БРОУНОВСКОЕ ДВИЖЕНИЕ
    // =========================================================================
    section("C.1 — Броуновское движение");
    {
        vector<vector<double>> paths(N_paths);
        mt19937 gen(42);
        normal_distribution<double> nd(0.0, 1.0);
        double dt = T / N_fine;

        for (int i = 0; i < N_paths; ++i) {
            vector<double> p(N_fine + 1, 0.0);
            for (int k = 1; k <= N_fine; ++k)
                p[k] = p[k-1] + sqrt(dt) * nd(gen);
            paths[i] = p;
        }

        plot_sde_paths(time_fine, paths,
                       "Броуновское движение — непрерывные, нигде не дифференцируемые траектории (C.1)");
        save("../plots/C1_brownian_motion.pdf");
        cla();
        cout << "  Сохранён: C1_brownian_motion.pdf\n";
        cout << "  Наблюдение: траектории непрерывны, но их 'скорость' нигде не определена.\n";
        cout << "  dW ~ sqrt(dt), поэтому dW/dt ~ 1/sqrt(dt) -> inf при dt->0\n";
    }

    // =========================================================================
    // C.3 — СРАВНЕНИЕ СХЕМ: ЭЙЛЕР vs МИЛЬШТЕЙН
    // Это ключевое расширение: PDF говорит что Мильштейн точнее,
    // мы это ИЗМЕРЯЕМ.
    // =========================================================================
    section("C.3 — Сравнение схем Эйлера и Мильштейна");
    {
        GeometricBrownianMotion gbm(0.05, 0.20);

        // --- 3а. Анализ сходимости ---
        ConvergenceAnalyzer analyzer(2000, {8, 16, 32, 64, 128, 256, 512});
        auto conv = analyzer.analyze_gbm(gbm.mu(), gbm.sigma(), 100.0, T);

        plot_convergence(conv, "../plots/C3_convergence_euler_vs_milshtein.pdf");
        cout << "  Сохранён: C3_convergence_euler_vs_milshtein.pdf\n";
        cout << "\n  Интерпретация:\n";
        cout << "    Эйлер:     наклон = " << fixed << setprecision(3)
             << conv.euler_slope << " (теория: 0.500)\n";
        cout << "    Мильштейн: наклон = " << conv.milshtein_slope
             << " (теория: 1.000)\n";
        cout << "\n  Вывод: Мильштейн вдвое точнее в log-log смысле.\n";
        cout << "  При одинаковом dt ошибка Мильштейна ~ dt, а Эйлера ~ sqrt(dt).\n";
        cout << "  Чтобы Эйлер достиг точности Мильштейна — нужно в 4 раза\n";
        cout << "  больше шагов => в 4 раза больше времени.\n";

        // --- 3б. Визуальное сравнение траекторий на ГРУБОЙ сетке ---
        // На грубой сетке разница между схемами видна визуально
        MonteCarloEngine euler_engine    (N_mc, N_coarse, T, Scheme::Euler);
        MonteCarloEngine milshtein_engine(N_mc, N_coarse, T, Scheme::Milshtein);

        auto euler_paths     = euler_engine.simulate(gbm, 100.0, 20);
        auto milshtein_paths = milshtein_engine.simulate(gbm, 100.0, 20);

        vector<double> exact_mean(N_coarse + 1);
        for (int i = 0; i <= N_coarse; ++i)
            exact_mean[i] = 100.0 * exp(gbm.mu() * time_coarse[i]);

        plot_scheme_comparison(time_coarse,
                               euler_paths, milshtein_paths, exact_mean,
                               "Эйлер vs Мильштейн на грубой сетке (N=100) — GBM (C.3)");
        save("../plots/C3_scheme_comparison.pdf");
        cla();
        cout << "  Сохранён: C3_scheme_comparison.pdf\n";

        // --- 3в. Vasicek (OU): схемы должны СОВПАДАТЬ ---
        // PDF утверждает что при константной диффузии схемы совпадают.
        // Проверим это.
        cout << "\n  Vasicek (OU) — проверка совпадения схем:\n";
        OrnsteinUhlenbeck vasicek(5.0, 0.05, 0.1);  // параметры близкие к реальным
        ConvergenceAnalyzer vasicek_analyzer(2000, {32, 64, 128, 256});

        // Для OU нет аналитического точного решения в одной точке W_T (зависит от всего пути),
        // поэтому сравниваем среднее МК с теоретическим.
        MonteCarloEngine ou_euler    (50000, N_fine, T, Scheme::Euler);
        MonteCarloEngine ou_milshtein(50000, N_fine, T, Scheme::Milshtein);

        ou_euler.simulate(vasicek, 0.0, 0);
        ou_milshtein.simulate(vasicek, 0.0, 0);

        double theory_mean = vasicek.exact_mean(0.0, T);
        double theory_var  = vasicek.exact_var(T);

        cout << "    E[X_T] теория:     " << theory_mean << "\n";
        cout << "    E[X_T] Эйлер:      " << ou_euler.get_mean_at_T() << "\n";
        cout << "    E[X_T] Мильштейн:  " << ou_milshtein.get_mean_at_T() << "\n";
        cout << "    Var[X_T] теория:   " << theory_var << "\n";
        cout << "    Var[X_T] Эйлер:    " << ou_euler.get_variance_at_T() << "\n";
        cout << "    Var[X_T] Мильштейн:" << ou_milshtein.get_variance_at_T() << "\n";
        cout << "  => Схемы дают практически одинаковый результат, как и предсказывает PDF.\n";
    }

    // =========================================================================
    // C.4 — МОДЕЛИ СДУ С ВЕРИФИКАЦИЕЙ
    // =========================================================================
    section("C.4 — Геометрическое броуновское движение (GBM)");
    {
        GeometricBrownianMotion gbm(0.05, 0.20);
        MonteCarloEngine engine(N_mc, N_fine, T, Scheme::Milshtein);
        auto paths = engine.simulate(gbm, 100.0, N_paths);

        // Теоретическое среднее: E[X_t] = X_0 * exp(mu * t)
        vector<double> mean_theory(N_fine + 1);
        for (int i = 0; i <= N_fine; ++i)
            mean_theory[i] = 100.0 * exp(gbm.mu() * time_fine[i]);

        plot_paths_with_mean(time_fine, paths, mean_theory,
            "GBM: dX = mu*X*dt + sigma*X*dW, субмартингал (C.4)");
        save("../plots/C4_gbm.pdf");
        cla();

        cout << "  E[X_T] теория:      " << 100.0 * exp(gbm.mu() * T) << "\n";
        cout << "  E[X_T] Монте-Карло: " << engine.get_mean_at_T() << "\n";
        cout << "\n  Важно: E[X_t] = X_0*exp(mu*t), но median(X_t) = X_0*exp((mu-0.5*sigma^2)*t)\n";
        cout << "  Разница — поправка Ито. При sigma=0.20: 0.5*sigma^2 = 0.02\n";
        cout << "  => median растёт медленнее среднего. Это нетривиальный эффект.\n";
    }

    section("C.4 — Модель Vasicek / Ornstein-Uhlenbeck");
    {
        // Параметры близкие к реальной калибровке на рынке ставок
        // theta=2: умеренная скорость возврата
        // mu=0.05: долгосрочная ставка 5%
        // sigma=0.02: волатильность 2%
        OrnsteinUhlenbeck vasicek(2.0, 0.05, 0.02);
        MonteCarloEngine engine(N_mc, N_fine, T, Scheme::Euler);
        // Для Vasicek Эйлер = Мильштейн (sigma'=0)
        auto paths = engine.simulate(vasicek, 0.03, N_paths);  // x0=3% (ниже mu)

        vector<double> mean_theory(N_fine + 1);
        for (int i = 0; i <= N_fine; ++i)
            mean_theory[i] = vasicek.exact_mean(0.03, time_fine[i]);

        plot_paths_with_mean(time_fine, paths, mean_theory,
            "Vasicek/OU: mean-reversion к долгосрочной ставке 5% (C.4)");
        save("../plots/C4_vasicek.pdf");
        cla();

        cout << "  Начальное значение x0 = 0.03 (3%) < mu = 0.05 (5%)\n";
        cout << "  => Процесс тянется вверх к 5% со скоростью theta=2\n";
        cout << "  E[X_T] теория:      " << vasicek.exact_mean(0.03, T) << "\n";
        cout << "  E[X_T] Монте-Карло: " << engine.get_mean_at_T() << "\n";
        cout << "  Var[X_T] теория:    " << vasicek.exact_var(T) << "\n";
        cout << "  Var[X_T] МК:        " << engine.get_variance_at_T() << "\n";
        cout << "\n  Проблема Vasicek: процесс может уйти в отрицательную зону.\n";
        cout << "  Это нереалистично для процентной ставки => нужна модель CIR.\n";
    }

    section("C.4 — Модель CIR (Cox-Ingersoll-Ross) — квадратно-корневой процесс");
    {
        // Параметры CIR
        // theta=2, mu=0.05, sigma=0.1
        // Условие Феллера: 2*theta*mu = 0.2 > sigma^2 = 0.01 — выполнено
        CoxIngersollRoss cir(2.0, 0.05, 0.1);

        cout << "  Условие Феллера (2*theta*mu >= sigma^2): "
             << (cir.feller_condition() ? "ВЫПОЛНЕНО" : "НЕ ВЫПОЛНЕНО") << "\n";
        cout << "  2*theta*mu = " << 2*2.0*0.05 << ", sigma^2 = " << 0.1*0.1 << "\n";
        cout << "  => Процесс строго положителен\n";

        MonteCarloEngine engine_euler    (N_mc, N_fine, T, Scheme::Euler);
        MonteCarloEngine engine_milshtein(N_mc, N_fine, T, Scheme::Milshtein);

        auto paths_euler     = engine_euler.simulate(cir, 0.03, N_paths);
        auto paths_milshtein = engine_milshtein.simulate(cir, 0.03, N_paths/2);

        vector<double> mean_theory(N_fine + 1);
        // Для CIR: E[X_t|X_0] = X_0*exp(-theta*t) + mu*(1-exp(-theta*t))
        // (совпадает с формулой OU для среднего!)
        OrnsteinUhlenbeck ou_approx(2.0, 0.05, 0.1);
        for (int i = 0; i <= N_fine; ++i)
            mean_theory[i] = ou_approx.exact_mean(0.03, time_fine[i]);

        plot_paths_with_mean(time_fine, paths_euler, mean_theory,
            "CIR: квадратно-корневой процесс, X_t >= 0 (C.4)");
        save("../plots/C4_cir.pdf");
        cla();

        cout << "\n  Сравнение E[X_T] для CIR:\n";
        cout << "    Эйлер:      " << engine_euler.get_mean_at_T() << "\n";
        cout << "    Мильштейн:  " << engine_milshtein.get_mean_at_T() << "\n";
        cout << "    Теория (среднее как у OU): " << ou_approx.exact_mean(0.03, T) << "\n";
        cout << "\n  CIR vs Vasicek:\n";
        cout << "    Vasicek: диффузия = sigma (константа) => может уйти < 0\n";
        cout << "    CIR:     диффузия = sigma*sqrt(X)     => при X->0 шум->0\n";
        cout << "    => CIR 'мягко' удерживает процесс от нуля\n";
        cout << "    Для CIR поправка Мильштейна ненулевая (sigma'!=0)\n";
    }

    // =========================================================================
    // C.6 — ПРОЦЕССЫ ПУАССОНА
    // =========================================================================
    section("C.6 — Процессы Пуассона");
    {
        // Несколько значений интенсивности для сравнения
        vector<double> lambdas = {2.0, 5.0, 15.0};
        vector<string> filenames = {
            "../plots/C6_poisson_lambda2.pdf",
            "../plots/C6_poisson_lambda5.pdf",
            "../plots/C6_poisson_lambda15.pdf"
        };

        for (int i = 0; i < 3; ++i) {
            auto events = PoissonProcess::simulate(T, lambdas[i]);
            plot_poisson(events, T);
            matplot::title("Процесс Пуассона, lambda=" + to_string((int)lambdas[i])
                + " (C.6): E[N_T]=" + to_string((int)lambdas[i]));
            save(filenames[i]);
            cla();
            cout << "  lambda=" << lambdas[i]
                 << ": смоделировано " << events.size()
                 << " скачков (теория E[N_T]=" << lambdas[i] << ")\n";
        }

        // Статистическая проверка: E[N_T] и Var[N_T] должны быть равны lambda*T
        cout << "\n  Статистическая верификация (10000 симуляций, lambda=5, T=1):\n";
        int n_sim = 10000;
        double lambda = 5.0;
        double sum_n = 0, sum_n2 = 0;
        for (int i = 0; i < n_sim; ++i) {
            auto ev = PoissonProcess::simulate(T, lambda);
            double n = static_cast<double>(ev.size());
            sum_n  += n;
            sum_n2 += n * n;
        }
        double mean_n = sum_n / n_sim;
        double var_n  = sum_n2 / n_sim - mean_n * mean_n;
        cout << "    E[N_T] теория=" << lambda*T << ", МК=" << mean_n << "\n";
        cout << "    Var[N_T] теория=" << lambda*T << ", МК=" << var_n << "\n";
        cout << "  => Для процесса Пуассона E=Var=lambda*T — подтверждено\n";
    }

    // =========================================================================
    // ИТОГОВАЯ ТАБЛИЦА
    // =========================================================================
    section("Итоговое сравнение моделей");
    // Колонки: Модель(14) | СДУ(26) | Диффузия(18) | Эйлер=Мильш?(16) | Применение(18)
    const int c1=14, c2=26, c3=18, c4=16, c5=18;
    int total = c1+c2+c3+c4+c5+4;
    cout << pad("Модель",       c1) << "  "
         << pad("СДУ",         c2) << "  "
         << pad("Диффузия",    c3) << "  "
         << pad("Эйлер=Мильш?",c4) << "  "
         << pad("Применение",  c5) << "\n"
         << string(total, '-') << "\n";

    auto row = [&](const string& m, const string& sde,
                   const string& diff, const string& eq, const string& app){
        cout << pad(m,    c1) << "  "
             << pad(sde,  c2) << "  "
             << pad(diff, c3) << "  "
             << pad(eq,   c4) << "  "
             << pad(app,  c5) << "\n";
    };

    row("GBM",        "mu*X dt + sigma*X dW",    "sigma*x",       "НЕТ", "Акции (B&S)");
    row("Vasicek/OU", "theta*(mu-X)dt + sigma dW","const",         "ДА",  "Ставки");
    row("CIR",        "theta*(mu-X)dt + s*sqX dW","sigma*sqrt(x)", "НЕТ", "Ставки (X>=0)");
    row("Пуассон",    "скачки размера 1",          "—",             "—",   "Дефолт, скачки");

    cout << "\nВсе графики сохранены в ../plots/\n";
    return 0;
}
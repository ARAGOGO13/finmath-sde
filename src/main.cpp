// =============================================================================
// Точка входа. Последовательно запускает все эксперименты и сохраняет
// графики в папку ../plots.
// =============================================================================

#include "console_utils.hpp"
#include "plotting.hpp"

#include "experiments/convergence_exp.hpp"
#include "experiments/gbm_inverse_exp.hpp"
#include "experiments/vasicek_exp.hpp"
#include "experiments/cir_exp.hpp"
#include "experiments/poisson_exp.hpp"
#include "experiments/ratio_gbm_exp.hpp"
#include "experiments/physical_pricing_exp.hpp"
#include "experiments/hedging_exp.hpp"

#include <cstdlib>
#include <iostream>
#include <matplot/matplot.h>

using namespace std;
using namespace matplot;

int main() {
    // Подавляем предупреждения Matplot++ в stderr
    freopen("/dev/null", "w", stderr);

    // Создаём папку для графиков, если её нет
    std::system("mkdir -p ../plots 2>/dev/null");

    const double T = 1.0;
    std::array<double,2> shared_ylim = compute_shared_mean_ylim(T);

    // C3: Сильная сходимость схем Эйлера и Мильштейна
    run_convergence_experiment();

    // C4: Модели стохастических процессов
    run_gbm_inverse_experiment(shared_ylim);
    run_vasicek_experiment(shared_ylim);
    run_cir_experiment(shared_ylim);

    // C6: Пуассоновские процессы
    run_poisson_experiment();

    // C_ratio: Отношение двух коррелированных GBM
    run_ratio_gbm_experiment();

    // G1: Ценообразование колла под физической мерой P vs риск-нейтральной Q
    run_physical_pricing_experiment();

    // G2: Дискретное дельта-хеджирование
    run_hedging_experiment();

    std::cout << "\nВсе эксперименты завершены. Графики сохранены в ../plots/\n";

    return 0;
}
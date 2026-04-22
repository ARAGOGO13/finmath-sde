#include "console_utils.hpp"
#include "plotting.hpp"

#include "experiments/convergence_exp.hpp"
#include "experiments/gbm_inverse_exp.hpp"
#include "experiments/vasicek_exp.hpp"
#include "experiments/cir_exp.hpp"
#include "experiments/poisson_exp.hpp"
#include "experiments/ratio_gbm_exp.hpp"
#include "experiments/girsanov_exp.hpp"
#include "experiments/hedging_exp.hpp"

#include <cstdlib>
#include <iostream>
#include <matplot/matplot.h>

using namespace std;
using namespace matplot;

int main() {
    freopen("/dev/null", "w", stderr);

    std::system("mkdir -p ../plots 2>/dev/null");

    const double T = 1.0;

    std::array<double,2> shared_ylim = compute_shared_mean_ylim(T);

    run_convergence_experiment();

    run_gbm_inverse_experiment(shared_ylim);
    run_vasicek_experiment(shared_ylim);
    run_cir_experiment(shared_ylim);

    run_poisson_experiment();

    run_ratio_gbm_experiment();

    run_girsanov_experiment();

    run_hedging_experiment();

    std::cout << "\nAll experiments completed. Plots saved to ../plots/\n";

    return 0;
}
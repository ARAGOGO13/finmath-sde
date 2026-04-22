#pragma once
#include <random>
#include <utility>
#include <vector>

// =============================================================================
// Модель пуассоновского процесса
//
// Пуассоновский процесс N_t с интенсивностью (параметром) λ > 0:
//   • N_0 = 0
//   • Приращения независимы и распределены по Пуассону: N_t - N_s ~ Poisson(λ(t-s))
//   • Времена между скачками независимы и имеют экспоненциальное распределение Exp(λ)
//
// Параметр λ определяет среднее число событий в единицу времени.
// Математическое ожидание и дисперсия: E[N_t] = Var[N_t] = λ t.
//
// Метод simulate возвращает последовательность событий в виде вектора пар
// (время события, накопленное число событий к этому моменту).
// =============================================================================

class PoissonProcess {
public:
    static std::vector<std::pair<double, int>> simulate(double T, double lambda) {
        std::mt19937 gen(std::random_device{}());
        std::exponential_distribution<double> exp_dist(lambda);
        std::vector<std::pair<double, int>> events;
        double t = 0.0;
        int count = 0;
        while (true) {
            t += exp_dist(gen);
            if (t >= T)
                break;
            events.emplace_back(t, ++count);
        }
        return events;
    }
};

// -----------------------------------------------------------------------------
// Генерация траектории пуассоновского процесса на заданной временной сетке
// -----------------------------------------------------------------------------
// Параметры:
//   tg     – вектор моментов времени (должен включать 0 и T, отсортирован)
//   lambda – интенсивность процесса (среднее число событий в единицу времени)
//
// Возвращает вектор значений N_t в те же моменты времени.
// Межскачковые интервалы моделируются экспоненциальным распределением.
// -----------------------------------------------------------------------------
inline std::vector<double> simulate_poisson_path(const std::vector<double>& tg, double lambda) {
    static thread_local std::mt19937 gen(std::random_device{}());
    std::exponential_distribution<double> exp_dist(lambda);
    std::vector<double> path(tg.size(), 0.0);
    double next_arrival = exp_dist(gen);
    int count = 0;
    for (size_t i = 0; i < tg.size(); ++i) {
        while (next_arrival <= tg[i]) {
            ++count;
            next_arrival += exp_dist(gen);
        }
        path[i] = static_cast<double>(count);
    }
    return path;
}


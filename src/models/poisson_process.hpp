// =============================================================================
// Пуассоновский процесс.
//
// Моделирует процесс с независимыми приращениями, распределёнными по Пуассону:
//   N_t ~ Poisson(λ t)
// Межскачковые интервалы распределены экспоненциально.
// =============================================================================

#pragma once
#include <random>
#include <vector>

// Генерация траектории пуассоновского процесса на заданной временной сетке.
// tg     – моменты времени (должен включать 0 и T, отсортирован по возрастанию)
// lambda – интенсивность процесса (среднее число событий в единицу времени)
// Возвращает вектор значений N_t в те же моменты времени.
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
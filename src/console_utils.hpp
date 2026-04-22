#pragma once

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// =============================================================================
// Утилиты для форматированного консольного вывода результатов экспериментов
// =============================================================================

// Длина строки в символах UTF-8 (для корректного выравнивания)
static size_t utf8_len(const std::string &s) {
    size_t len = 0;
    for (unsigned char c : s)
        if ((c & 0xC0) != 0x80) ++len;
    return len;
}

// Дополняет строку пробелами до заданной ширины
static std::string pad(const std::string &s, int w) {
    int sp = w - static_cast<int>(utf8_len(s));
    return s + std::string(std::max(sp, 1), ' ');
}

// Печатает разделитель из дефисов заданной длины (с отступом в 2 пробела)
static void print_separator(int length = 70) {
    std::cout << "  " << std::string(length, '-') << "\n";
}

// Печатает заголовок раздела
static void section(const std::string &title, const std::string &sub = "") {
    std::cout << "\n===== " << title << " =====\n";
    if (!sub.empty())
        std::cout << "      " << sub << "\n";
}

// Печатает заголовок таблицы с произвольными колонками
static void print_table_header(const std::vector<std::string>& headers,
                               const std::vector<int>& widths) {
    std::cout << "  ";
    for (size_t i = 0; i < headers.size(); ++i) {
        std::cout << pad(headers[i], widths[i]);
    }
    std::cout << "\n";
    int total = 0;
    for (int w : widths) total += w;
    print_separator(total);
}

// Печатает строку таблицы (все значения как строки)
static void print_table_row(const std::vector<std::string>& values,
                            const std::vector<int>& widths) {
    std::cout << "  ";
    for (size_t i = 0; i < values.size(); ++i) {
        std::cout << pad(values[i], widths[i]);
    }
    std::cout << "\n";
}

// Заголовок таблицы сравнения MC и теории (устаревший вариант, сохранён для совместимости)
static void table_header(int lw = 16) {
    std::cout << "\n  " << pad("", lw + 2) << std::setw(12) << "MC"
              << "  " << std::setw(12) << "Theory"
              << "  " << std::setw(8) << "Err %" << "\n"
              << "  " << std::string(lw + 40, '-') << "\n";
}

// Строка таблицы с расчётом относительной ошибки
static void table_row(const std::string &label, double mc, double theory, int lw = 16) {
    double err = std::abs(mc - theory) / (std::abs(theory) + 1e-15) * 100.0;
    std::cout << "  " << pad(label, lw)
              << "  " << std::fixed << std::setprecision(6) << std::setw(12) << mc
              << "  " << std::setw(12) << theory
              << "  " << std::setprecision(3) << std::setw(7) << err << "%\n";
}

// Выводит список сохранённых файлов
static void print_saved(const std::vector<std::string> &files) {
    for (const auto &f : files)
        std::cout << "  " << f << "\n";
}
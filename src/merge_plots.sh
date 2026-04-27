#!/usr/bin/env bash
set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PLOTS_DIR="$SCRIPT_DIR/../plots"

if [ ! -d "$PLOTS_DIR" ]; then
    echo "❌ Папка $PLOTS_DIR не найдена"
    exit 1
fi

cd "$PLOTS_DIR"

# Объединяет PDF-файлы в одну многостраничную PDF, удаляет исходники
merge_and_remove() {
    local output_name="$1"
    shift
    local input_files=("$@")

    local existing=()
    for f in "${input_files[@]}"; do
        [ -f "$f" ] && existing+=("$f")
    done

    if [ ${#existing[@]} -eq 0 ]; then
        echo "  ⚠️  Нет файлов для $output_name — пропускаем"
        return 0
    fi

    if pdfunite "${existing[@]}" "$output_name" 2>/dev/null; then
        for f in "${existing[@]}"; do rm -f "$f"; done
        echo "  ✓  $output_name  (${#existing[@]} страниц)"
    else
        echo "  ❌  Не удалось создать $output_name"
        return 1
    fi
}

# Склеивает до 4 PDF в таблицу 2×2 на одной странице, удаляет исходники
merge_2x2() {
    local output_name="$1"
    shift
    local input_files=("$@")

    local existing=()
    for f in "${input_files[@]}"; do
        [ -f "$f" ] && existing+=("$f")
    done

    if [ ${#existing[@]} -eq 0 ]; then
        echo "  ⚠️  Нет файлов для $output_name — пропускаем"
        return 0
    fi

    local ok=0

    local MAKE_GRID="$SCRIPT_DIR/make_grid.py"
    if [ -f "$MAKE_GRID" ]; then
        if python3 "$MAKE_GRID" "$output_name" "${existing[@]}" 2>/dev/null; then
            ok=1
        fi
    fi

    if [ $ok -eq 0 ] && command -v pdfjam &>/dev/null; then
        if pdfjam --quiet --nup '2x2' --landscape \
                  "${existing[@]}" --outfile "$output_name" 2>/dev/null; then
            ok=1
        fi
    fi

    if [ $ok -eq 1 ]; then
        for f in "${existing[@]}"; do rm -f "$f"; done
        echo "  ✓  $output_name  (2×2 grid)"
    else
        if pdfunite "${existing[@]}" "$output_name" 2>/dev/null; then
            for f in "${existing[@]}"; do rm -f "$f"; done
            echo "  ✓  $output_name  (${#existing[@]} стр. подряд — fallback)"
        fi
    fi
}

echo ""
echo "── C4: мат. ожидание — 2×2 (GBM | InvGBM / Vasicek | CIR) ─────────────────"
merge_2x2 "C4_mean_2x2.pdf" \
    "C4_gbm_mean.pdf" \
    "C4_gbm_inv_mean.pdf" \
    "C4_vasicek_mean.pdf" \
    "C4_cir_mean.pdf"

echo ""
echo "── C4: дисперсия — 2×2 (GBM | InvGBM / Vasicek | CIR) ─────────────────────"
merge_2x2 "C4_var_2x2.pdf" \
    "C4_gbm_var.pdf" \
    "C4_gbm_inv_var.pdf" \
    "C4_vasicek_var.pdf" \
    "C4_cir_var.pdf"

echo ""
echo "── C3: Euler vs Milshtein (уже один файл) ───────────────────────────────────"
[ -f "C3_convergence.pdf" ] && echo "  ✓  C3_convergence.pdf  (без изменений)"

echo ""
echo "── C6: Poisson ──────────────────────────────────────────────────────────────"
merge_and_remove "C6_poisson_merged.pdf" \
    "C6_poisson_paths.pdf" \
    "C6_poisson_mean.pdf" \
    "C6_poisson_var.pdf"

echo ""
echo "── C_ratio: траектории + ±1σ  —  2×2 + 1 ───────────────────────────────────"
MAKE_GRID="$SCRIPT_DIR/make_grid.py"
_ratio_files=(
    "C_ratio_paths_m10.pdf"
    "C_ratio_paths_m05.pdf"
    "C_ratio_paths_00.pdf"
    "C_ratio_paths_p05.pdf"
    "C_ratio_paths_p10.pdf"
)
_ratio_exist=()
for f in "${_ratio_files[@]}"; do
    [ -f "$f" ] && _ratio_exist+=("$f")
done

if [ ${#_ratio_exist[@]} -gt 0 ]; then
    _ok=0
    if [ -f "$MAKE_GRID" ]; then
        if python3 "$MAKE_GRID" "C_ratio_paths_grid.pdf" \
                  --layout 2x2+1 "${_ratio_exist[@]}" 2>/dev/null; then
            _ok=1
        fi
    fi
    if [ $_ok -eq 1 ]; then
        for f in "${_ratio_exist[@]}"; do rm -f "$f"; done
        echo "  ✓  C_ratio_paths_grid.pdf  (2×2 + 1)"
    else
        pdfunite "${_ratio_exist[@]}" "C_ratio_paths_grid.pdf" 2>/dev/null
        for f in "${_ratio_exist[@]}"; do rm -f "$f"; done
        echo "  ✓  C_ratio_paths_grid.pdf  (${#_ratio_exist[@]} стр. подряд — fallback)"
    fi
fi

echo ""
echo "── C_ratio: Var[V_t] все rho (уже один файл) ────────────────────────────────"
[ -f "C_ratio_var_all.pdf" ] && echo "  ✓  C_ratio_var_all.pdf  (без изменений)"

echo ""
echo "── G1: Physical pricing vs RN ───────────────────────────────────────────────"

# G1b: 4 payoff-гистограммы → 2×2
merge_2x2 "G1_payoff_2x2.pdf" \
    "G1_payoff_mu_0.pdf" \
    "G1_payoff_mu_1.pdf" \
    "G1_payoff_mu_2.pdf" \
    "G1_payoff_mu_3.pdf"

# G1c: 4 P&L-гистограммы → 2×2
merge_2x2 "G1_pnl_2x2.pdf" \
    "G1_pnl_mu_0.pdf" \
    "G1_pnl_mu_1.pdf" \
    "G1_pnl_mu_2.pdf" \
    "G1_pnl_mu_3.pdf"

# Финальный merge — каждый аналитический график на отдельной странице:
#   стр. 1 — G1a: цена C0_P(mu) vs mu
#   стр. 2 — G1d-stacked: OTM + ITM-loss + profit = 100%
#   стр. 3 — G1d-epnl: E[P&L] = C0_Q - C0_P vs mu
#   стр. 4 — G1e: P(profit|ITM) и P(loss|ITM) vs mu
#   стр. 5 — G1b: payoff-гистограммы 2×2
#   стр. 6 — G1c: P&L-гистограммы 2×2
merge_and_remove "G1_physical_pricing_merged.pdf" \
    "G1_physical_vs_rn_prices.pdf" \
    "G1_pnl_stacked.pdf" \
    "G1_expected_pnl.pdf" \
    "G1_conditional_probs.pdf" \
    "G1_payoff_2x2.pdf" \
    "G1_pnl_2x2.pdf"

echo ""
echo "── G2: Delta hedging (4 графика → 2×2) ─────────────────────────────────────"
merge_2x2 "G2_hedging_merged.pdf" \
    "G2_hedging_pnl_hist.pdf" \
    "G2_hedging_std_vs_dt.pdf" \
    "G2_hedging_vol_mismatch_mean.pdf" \
    "G2_hedging_vol_mismatch_std.pdf"

# ── Итог ──────────────────────────────────────────────────────────────────────
echo ""
echo "📁 Итоговые файлы в $PLOTS_DIR:"
ls -1 *.pdf 2>/dev/null \
    | sed 's/^/    /' \
    || echo "    (нет PDF)"
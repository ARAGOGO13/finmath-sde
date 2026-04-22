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
# Использует make_grid.py (Python3 + Pillow + pdftoppm)
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

    # Попытка 1: Python-скрипт make_grid.py (Pillow + pdftoppm)
    local MAKE_GRID="$SCRIPT_DIR/make_grid.py"
    if [ -f "$MAKE_GRID" ]; then
        if python3 "$MAKE_GRID" "$output_name" "${existing[@]}" 2>/dev/null; then
            ok=1
        fi
    fi

    # Попытка 2: pdfjam (MacTeX / TeX Live)
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
        # Fallback: страницы подряд
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
echo "── C_ratio: траектории + ±1σ  —  2×2 + 1 (rho: -1 -½ / 0 +½ / +1 centre) ──"
# 5 plots: 2x2 grid + 1 centered below (via make_grid.py --layout 2x2+1)
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
        # fallback: sequential pages
        pdfunite "${_ratio_exist[@]}" "C_ratio_paths_grid.pdf" 2>/dev/null
        for f in "${_ratio_exist[@]}"; do rm -f "$f"; done
        echo "  ✓  C_ratio_paths_grid.pdf  (${#_ratio_exist[@]} стр. подряд — fallback)"
    fi
fi

echo ""
echo "── C_ratio: Var[V_t] все rho (уже один файл) ────────────────────────────────"
[ -f "C_ratio_var_all.pdf" ] && echo "  ✓  C_ratio_var_all.pdf  (без изменений)"

echo ""
echo "── G1: Girsanov verification (2 графика → склейка) ───────────────────────────"
merge_and_remove "G1_girsanov_merged.pdf" \
    "G1_girsanov_prices.pdf" \
    "G1_girsanov_variance.pdf"

echo ""
echo "── G2: Delta hedging (3 графика → склейка) ──────────────────────────────────"
merge_and_remove "G2_hedging_merged.pdf" \
    "G2_hedging_pnl_hist.pdf" \
    "G2_hedging_std_vs_dt.pdf" \
    "G2_hedging_vol_mismatch.pdf"

# ── Итог ──────────────────────────────────────────────────────────────────────
echo ""
echo "📁 Итоговые файлы в $PLOTS_DIR:"
ls -1 *.pdf 2>/dev/null \
    | sed 's/^/    /' \
    || echo "    (нет PDF)"
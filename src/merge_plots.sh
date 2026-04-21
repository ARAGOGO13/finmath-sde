#!/bin/bash
set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PLOTS_DIR="$SCRIPT_DIR/../plots"

if [ ! -d "$PLOTS_DIR" ]; then
    echo "❌ Папка $PLOTS_DIR не найдена"
    exit 1
fi

cd "$PLOTS_DIR"

merge_and_remove() {
    local output_name="$1"
    shift
    local input_files=("$@")


    if pdfunite "${input_files[@]}" "$output_name" 2>/dev/null; then
        for f in "${input_files[@]}"; do
            if [ -f "$f" ]; then
                rm "$f"
            fi
        done
    else
        return 1
    fi
}

merge_and_remove "C4_gbm_merged.pdf" \
    "C4_gbm_mean.pdf" \
    "C4_gbm_var.pdf" \
    "C4_gbm_inv_mean.pdf" \
    "C4_gbm_inv_var.pdf"

merge_and_remove "C4_cir_merged.pdf" \
    "C4_cir_mean.pdf" \
    "C4_cir_var.pdf"

merge_and_remove "C4_vasicek_merged.pdf" \
    "C4_vasicek_mean.pdf" \
    "C4_vasicek_var.pdf"

# 4. Poisson (3 файла: paths, mean, var)
merge_and_remove "C6_poisson_merged.pdf" \
    "C6_poisson_paths.pdf" \
    "C6_poisson_mean.pdf" \
    "C6_poisson_var.pdf"



ls -1 *.pdf | grep -E "C.*_merged\.pdf|C3_convergence\.pdf" || echo "  (нет)"
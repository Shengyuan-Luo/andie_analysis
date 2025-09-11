#!/bin/bash
# submit_summarize_lcc_trend_dual.sh
set -euo pipefail

# 自动上溯 Sample-* 目录
ascend_sample_dir() {
  local cur="$(pwd -P)"
  while [[ "$cur" != "/" ]]; do
    local base="$(basename "$cur")"
    if [[ "$base" == Sample-* ]]; then
      echo "$cur"; return 0
    fi
    cur="$(dirname "$cur")"
  done
  return 1
}

BASE_DIR="$(ascend_sample_dir || true)"
if [[ -z "${BASE_DIR:-}" ]]; then
  HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
  BASE_DIR="$(cd "$HERE" && ascend_sample_dir || true)"
fi
if [[ -z "${BASE_DIR:-}" ]]; then
  echo "[ERROR] 无法从当前路径上溯找到 Sample-* 目录。" >&2; exit 2
fi

PY="$BASE_DIR/summarize_lcc_trend_dual.py"
LOG="$BASE_DIR/viz_results_ver2_dual_summary/logs"
mkdir -p "$LOG"

cat > "$LOG/sbatch_summarize_lcc.sh" <<'EOF'
#!/bin/bash
#SBATCH --job-name=lcc_trend_dual
#SBATCH --output=__LOG__/lcc_trend_dual_%j.out
#SBATCH --error=__LOG__/lcc_trend_dual_%j.err
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

PYTHON=$(command -v python3 || true)
if [[ -z "$PYTHON" ]]; then echo "python3 not found in PATH"; exit 2; fi

"$PYTHON" "__PY__" --thresholds "1.5,1.75,2.0,2.25,2.5"
EOF

sed -i -e "s|__LOG__|$LOG|g" -e "s|__PY__|$PY|g" "$LOG/sbatch_summarize_lcc.sh"
chmod +x "$LOG/sbatch_summarize_lcc.sh"
sbatch "$LOG/sbatch_summarize_lcc.sh"
echo "[SUBMITTED] summarize_lcc_trend_dual → $LOG"


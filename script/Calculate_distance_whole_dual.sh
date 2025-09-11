#!/bin/bash
# Calculate_distance_whole_dual.sh
# Usage: bash Calculate_distance_whole_dual.sh [NUM_TASKS]
# NUM_TASKS≤1: 单任务，用KD-tree；>1：切分为NUM_TASKS个SLURM任务
set -euo pipefail

TASKS=${1:-1}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

# 上溯找到 Sample- 目录（也许 BASE_DIR 就是）
ascend_to_sample_dir() {
  local cur="$1"
  while true; do
    local base="$(basename "$cur")"
    if [[ "$base" == Sample-* ]]; then echo "$cur"; return 0; fi
    local parent="$(dirname "$cur")"
    if [[ "$parent" == "$cur" ]]; then
      echo "ERROR: 未能从 $1 上溯找到 Sample-* 目录" >&2; exit 1
    fi
    cur="$parent"
  done
}
SAMPLE_DIR="$(ascend_to_sample_dir "$SCRIPT_DIR")"
SAMPLE_NAME="$(basename "$SAMPLE_DIR")"

# 定位两类输入（优先精确命名，其次glob）
EU_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}.euchromatin_cluster.txt"
[[ -f "$EU_FILE" ]] || EU_FILE="$(ls -1 "${SAMPLE_DIR}"/*euchromatin_cluster.txt 2>/dev/null | head -n1 || true)"
H3_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}.h3k4me3_cluster.txt"
[[ -f "$H3_FILE" ]] || H3_FILE="$(ls -1 "${SAMPLE_DIR}"/*h3k4me3_cluster.txt 2>/dev/null | head -n1 || true)"

# 双份输出根目录
OUT_EU="${BASE_DIR}/Whole_genome_distance_dual_euchr"
OUT_H3="${BASE_DIR}/Whole_genome_distance_dual_h3k4"
mkdir -p "$OUT_EU" "$OUT_H3" "${SCRIPT_DIR}/tmp"

# 统计有效行数（剔除可能的表头：homolog/chrom）
count_bins() {
  local f="$1"
  awk 'NR==1{if($1=="homolog"||$1=="chrom") next}{c++}END{print c+0}' "$f"
}

echo "========== PLAN（输入读取结构） =========="
echo "Sample name : ${SAMPLE_NAME}"
echo "Sample dir  : ${SAMPLE_DIR}"
echo "Inputs:"
printf "  euchr     : %s\n" "${EU_FILE:-未找到}"
printf "  h3k4me3   : %s\n" "${H3_FILE:-未找到}"
echo "Outputs:"
echo "  euchr → ${OUT_EU}"
echo "  h3k4 → ${OUT_H3}"
echo "Tasks per file: ${TASKS}"
echo "========================================="

submit_one() {
  local label="$1"               # euchr | h3k4
  local infile="$2"
  local outroot="$3"
  [[ -f "$infile" ]] || { echo "[Skip] ${label}: 输入不存在"; return; }

  local fname base
  fname="$(basename "$infile")"
  base="${fname%.euchromatin_cluster.txt}"
  if [[ "$base" == "$fname" ]]; then
    base="${fname%.h3k4me3_cluster.txt}"
  fi

  echo "[WholeGenome-${label}] Computing distances for $fname (tasks=$TASKS)…"
  if [[ ${TASKS} -le 1 ]]; then
    local out_file="${outroot}/${base}_distance_filtered.txt"
    python3 "${SCRIPT_DIR}/Calculate_distance_whole_dual.py" \
      "$infile" "$out_file" --source_label "$label"
    [[ $? -eq 0 ]] || { echo "[Error] compute failed for $fname"; exit 1; }
  else
    local total_bins bins_per_task part_dir
    total_bins=$(count_bins "$infile")
    bins_per_task=$(( (total_bins + TASKS - 1) / TASKS ))
    part_dir="${outroot}/partial_${base}"
    mkdir -p "$part_dir"

    for (( t=1; t<=TASKS; t++ )); do
      start_index=$(( (t-1) * bins_per_task ))
      [[ $start_index -lt $total_bins ]] || break
      end_index=$(( t * bins_per_task - 1 ))
      [[ $end_index -lt $total_bins ]] || end_index=$(( total_bins - 1 ))

      part_out="${part_dir}/${base}_dist_part${t}.txt"
      job_script=$(mktemp /tmp/dist_whole_${label}_${t}_XXXX.sh)
      cat > "$job_script" << EOF
#!/bin/bash
#SBATCH --job-name=dist_whole_${label}_${t}
#SBATCH --output=${SCRIPT_DIR}/tmp/whole_${label}_${t}.out
#SBATCH --error=${SCRIPT_DIR}/tmp/whole_${label}_${t}.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G

python3 "${SCRIPT_DIR}/Calculate_distance_whole_dual.py" \
  "$infile" "$part_out" --source_label "$label" --range $start_index $end_index
EOF
      sbatch "$job_script" >/dev/null
      rm -f "$job_script"
      echo "  Submitted ${label} task $t for bins $start_index–$end_index"
    done

    echo "[WholeGenome-${label}] All ${TASKS} tasks submitted for ${base}."
    echo "  Merge after completion:"
    echo "    cat ${part_dir}/${base}_dist_part*.txt > ${outroot}/${base}_distance_filtered.txt"
  fi
}

# 分别提交 euchr 与 h3k4
submit_one "euchr" "${EU_FILE:-}" "$OUT_EU"
submit_one "h3k4"  "${H3_FILE:-}" "$OUT_H3"


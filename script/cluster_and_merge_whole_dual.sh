#!/bin/bash
# cluster_and_merge_whole_dual.sh
# Usage: bash cluster_and_merge_whole_dual.sh [THRESHOLDS...]
# - 严格：只读带后缀 Whole_genome_distance_dual_{euchr|h3k4}/
# - 双输出：graph_matrix_dual_{euchr|h3k4}/
set -euo pipefail

if [ $# -gt 0 ]; then
  THRESHOLDS=("$@")
else
  THRESHOLDS=(1 1.25 1.5 1.75 2 2.25 2.5 2.75 3)
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# 上溯找到 Sample- 目录与样本名
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
BASE_NAME="$(basename "$SAMPLE_DIR")"

LOG_DIR="${SCRIPT_DIR}/logs_dual"
mkdir -p "$LOG_DIR"

# 输入（严格只读 dual 目录）
WHOLE_DIR_EU="${PROJECT_DIR}/Whole_genome_distance_dual_euchr"
WHOLE_DIR_H3="${PROJECT_DIR}/Whole_genome_distance_dual_h3k4"
MERGED_EU="${WHOLE_DIR_EU}/${BASE_NAME}_distance_filtered.txt"
MERGED_H3="${WHOLE_DIR_H3}/${BASE_NAME}_distance_filtered.txt"

# 对应 cluster 文件（all_bins 模式需要）
CLUSTER_EU="$(ls -1 "${SAMPLE_DIR}"/*euchromatin_cluster.txt 2>/dev/null | head -n1 || true)"
CLUSTER_H3="$(ls -1 "${SAMPLE_DIR}"/*h3k4me3_cluster.txt 2>/dev/null | head -n1 || true)"

NODE_MODE="all_bins"  # 你原来的默认

echo "========== PLAN（输入读取结构） =========="
echo "Sample name : ${BASE_NAME}"
echo "Sample dir  : ${SAMPLE_DIR}"
echo "Inputs:"
printf "  euchr  distance : %s\n" "${MERGED_EU:-未找到}"
printf "  euchr  cluster  : %s\n" "${CLUSTER_EU:-未找到}"
printf "  h3k4   distance : %s\n" "${MERGED_H3:-未找到}"
printf "  h3k4   cluster  : %s\n" "${CLUSTER_H3:-未找到}"
echo "Outputs:"
echo "  euchr → ${PROJECT_DIR}/graph_matrix_dual_euchr/"
echo "  h3k4  → ${PROJECT_DIR}/graph_matrix_dual_h3k4/"
echo "Thresholds: ${THRESHOLDS[*]}"
echo "Node mode : ${NODE_MODE}"
echo "========================================="

# 严格校验输入
[[ -f "$MERGED_EU" ]] || { echo "[ERROR] 缺少 euchr 距离文件：$MERGED_EU"; exit 2; }
[[ -f "$MERGED_H3" ]] || { echo "[ERROR] 缺少 h3k4 距离文件：$MERGED_H3"; exit 2; }
if [[ "$NODE_MODE" == "all_bins" ]]; then
  [[ -f "$CLUSTER_EU" ]] || { echo "[ERROR] 缺少 euchr cluster 文件"; exit 2; }
  [[ -f "$CLUSTER_H3" ]] || { echo "[ERROR] 缺少 h3k4 cluster 文件"; exit 2; }
fi

submit_group() {
  local label="$1"    # euchr | h3k4
  local dist_file="$2"
  local cluster_file="$3"
  local deps=()
  for d in "${THRESHOLDS[@]}"; do
    tag="${d//./p}"
    job_script="$(mktemp /tmp/cluster_${label}_${tag}_XXXX.sh)"
    cat > "$job_script" << EOF
#!/bin/bash
#SBATCH --job-name=cluster_${label}_${tag}
#SBATCH --output=${LOG_DIR}/cluster_${label}_${tag}.out
#SBATCH --error=${LOG_DIR}/cluster_${label}_${tag}.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --chdir=${SCRIPT_DIR}

python3 "${SCRIPT_DIR}/analyze_graph_parallel_dual.py" \
  --distance_file "${dist_file}" \
  --threshold ${d} \
  --output_prefix whole${d} \
  --node_mode "${NODE_MODE}" \
  --source_label "${label}" \
  --cluster_file "${cluster_file}" \
  --out_root "${PROJECT_DIR}/graph_matrix_dual_${label}"
EOF
    jid=$(sbatch "$job_script" | awk '{print $NF}')
    rm -f "$job_script"
    deps+=("$jid")
    echo "Submitted ${label} clustering thr=${d} (jobid=${jid})"
  done

  local dep_str; dep_str=$(IFS=:; echo "${deps[*]}")
  merge_job="$(mktemp /tmp/merge_${label}_XXXX.sh)"
  cat > "$merge_job" << EOF
#!/bin/bash
#SBATCH --job-name=merge_${label}
#SBATCH --output=${LOG_DIR}/merge_${label}.out
#SBATCH --error=${LOG_DIR}/merge_${label}.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --chdir=${SCRIPT_DIR}

python3 "${SCRIPT_DIR}/merge_components_dual.py" --source_label "${label}"
EOF
  sbatch --dependency=afterok:${dep_str} "$merge_job" >/dev/null
  rm -f "$merge_job"
  echo "Scheduled merge for ${label} after all thresholds complete."
}

submit_group "euchr" "$MERGED_EU" "$CLUSTER_EU"
submit_group "h3k4"  "$MERGED_H3" "$CLUSTER_H3"

echo "All cluster jobs submitted; each label has its own merge job."
echo "Logs under ${LOG_DIR}/"


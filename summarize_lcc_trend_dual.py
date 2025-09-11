#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, argparse, re, glob
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- 公共：自动识别 Sample-* 目录 ----------
def ascend_sample_dir(p: str):
    cur = os.path.abspath(p)
    while True:
        base = os.path.basename(cur)
        if base.startswith("Sample-"):
            return cur, base
        parent = os.path.dirname(cur)
        if parent == cur:
            raise RuntimeError(f"无法从 {p} 上溯找到 Sample-* 目录")
        cur = parent

# ---------- 路径构造 & 容错 ----------
def thr_variants(thr_str: str):
    """为阈值提供若干候选文本形式，用于路径容错：'2.0' → ['2.0','2']；'1.5' → ['1.5']"""
    try:
        v = float(thr_str)
        if v.is_integer():
            return [thr_str, str(int(v))]
        else:
            return [thr_str]
    except Exception:
        return [thr_str]

def find_metrics_file(base_dir: str, base_name: str, label: str, thr_str: str):
    """
    优先寻找：
      graph_matrix_dual_<label>/whole{tag}/<base>_whole{tag}_metrics.txt
    对 2.0 这类整数阈值，既尝试 '2.0' 也尝试 '2' 版本；
    若以上都没找到，则退化为在 whole{tag} 目录里 glob *_metrics.txt 并比对文件里 'threshold' 值。
    返回 (metrics_file_path or None, threshold_value or None, num_nodes or None, largest_cc_size or None)
    """
    root = os.path.join(base_dir, f"graph_matrix_dual_{label}")
    for tag in thr_variants(thr_str):
        wdir = os.path.join(root, f"whole{tag}")
        if not os.path.isdir(wdir):
            continue
        # 尝试标准文件名
        candidate = os.path.join(wdir, f"{base_name}_whole{tag}_metrics.txt")
        if os.path.isfile(candidate):
            vals = parse_metrics(candidate)
            if vals["threshold"] is None or float(vals["threshold"]) == float(thr_str):
                return candidate, vals["threshold"], vals["num_nodes"], vals["largest_cc_size"]
        # 退化：glob
        for mf in glob.glob(os.path.join(wdir, "*_metrics.txt")):
            vals = parse_metrics(mf)
            if vals["threshold"] is not None and abs(float(vals["threshold"]) - float(thr_str)) < 1e-9:
                return mf, vals["threshold"], vals["num_nodes"], vals["largest_cc_size"]
    return None, None, None, None

def parse_metrics(path: str):
    """
    解析 metrics 文本，容错读取两列 'key  value'；忽略多余内容。
    关心：threshold, num_nodes, largest_cc_size
    """
    out = {"threshold": None, "num_nodes": None, "largest_cc_size": None}
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"): 
                    continue
                parts = re.split(r"\s+", line)
                if len(parts) < 2:
                    continue
                key = parts[0].strip().lower()
                val = parts[1].strip()
                if key == "threshold":
                    try: out["threshold"] = float(val)
                    except: pass
                elif key == "num_nodes":
                    try: out["num_nodes"] = int(float(val))
                    except: pass
                elif key == "largest_cc_size":
                    try: out["largest_cc_size"] = int(float(val))
                    except: pass
    except Exception as e:
        print(f"[WARN] 解析失败：{path} -> {e}", file=sys.stderr)
    return out

# ---------- 作图 ----------
def plot_grouped_bars(thrs, ratios_eu, ratios_h3, out_png, color_eu="#b0d9a5", color_h3="#fdd379"):
    """
    画一张“分组柱状图”：x 轴为阈值，euchr/h3k4 两组对比各自最大 CC 占比
    """
    import numpy as np
    x = np.arange(len(thrs))
    width = 0.36

    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    ax.bar(x - width/2, ratios_eu, width, label="euchr", color=color_eu, alpha=0.9)
    ax.bar(x + width/2, ratios_h3, width, label="h3k4me3", color=color_h3, alpha=0.9)

    ax.set_title("Largest connected component ratio vs. distance threshold")
    ax.set_xlabel("distance threshold")
    ax.set_ylabel("largest_cc_size / num_nodes")
    ax.set_xticks(x)
    ax.set_xticklabels(thrs)
    ax.set_ylim(0, 1.05)
    ax.legend(loc="upper left", frameon=True)
    ax.grid(axis="y", linestyle="--", alpha=0.3)

    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"[OK] 保存：{out_png}")

def plot_single_bars(thrs, ratios, out_png, title, color="#b0d9a5"):
    import numpy as np
    x = np.arange(len(thrs))

    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    ax.bar(x, ratios, width=0.55, color=color, alpha=0.9)

    ax.set_title(title)
    ax.set_xlabel("distance threshold")
    ax.set_ylabel("largest_cc_size / num_nodes")
    ax.set_xticks(x)
    ax.set_xticklabels(thrs)
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", linestyle="--", alpha=0.3)

    plt.tight_layout()
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"[OK] 保存：{out_png}")

# ---------- 主程序 ----------
def main():
    ap = argparse.ArgumentParser(
        description="统计 graph_matrix_dual_euchr/h3k4 下 whole{thr} 的 largest_cc 占比并画柱状图（自动识别样本名与路径）"
    )
    ap.add_argument("--thresholds", default="1.5,1.75,2.0,2.25,2.5",
                    help="逗号分隔的阈值列表（字符串），默认：1.5,1.75,2.0,2.25,2.5")
    ap.add_argument("--color_euchr", default="#b0d9a5", help="euchr 柱状颜色")
    ap.add_argument("--color_h3k4", default="#fdd379", help="h3k4 柱状颜色")
    ap.add_argument("--out_suffix", default="ver2_dual", help="输出目录后缀，避免覆盖旧结果")
    args = ap.parse_args()

    sample_dir, sample_name = ascend_sample_dir(os.getcwd())
    thrs = [t.strip() for t in args.thresholds.split(",") if t.strip()]
    print(f"[INFO] Sample dir: {sample_dir}")
    print(f"[INFO] Sample name: {sample_name}")
    print(f"[INFO] Thresholds: {thrs}")

    # 逐 label 收集比例
    ratios = {"euchr": [], "h3k4": []}
    used_thr = []  # 真正找到文件的阈值顺序（两侧都尽量齐全；缺失用 None）

    for t in thrs:
        # euchr
        me_path_eu, th_eu, n_eu, lcc_eu = find_metrics_file(sample_dir, sample_name, "euchr", t)
        # h3k4
        me_path_h3, th_h3, n_h3, lcc_h3 = find_metrics_file(sample_dir, sample_name, "h3k4", t)

        # 打印来源
        print(f"\n[SCAN] thr={t}")
        print(f"  euchr metrics: {me_path_eu or 'NOT FOUND'}")
        print(f"  h3k4  metrics: {me_path_h3 or 'NOT FOUND'}")

        # 计算比例（缺失填 None）
        if n_eu and lcc_eu is not None and n_eu > 0:
            ratios["euchr"].append(lcc_eu / n_eu)
        else:
            ratios["euchr"].append(None)
        if n_h3 and lcc_h3 is not None and n_h3 > 0:
            ratios["h3k4"].append(lcc_h3 / n_h3)
        else:
            ratios["h3k4"].append(None)

        used_thr.append(t)

    # 把 None 转为 0（保持柱状图不报错），同时记录缺失信息
    def fill_none_to_zero(arr, label):
        out = []
        for i, v in enumerate(arr):
            if v is None:
                print(f"[WARN] {label} 在 thr={used_thr[i]} 缺少数据，按 0 绘制。", file=sys.stderr)
                out.append(0.0)
            else:
                out.append(float(v))
        return out

    eu_vals = fill_none_to_zero(ratios["euchr"], "euchr")
    h3_vals = fill_none_to_zero(ratios["h3k4"], "h3k4")

    # 输出目录结构（全部带 ver2_dual* 后缀，避免覆盖旧结果）
    out_group = os.path.join(sample_dir, f"viz_results_{args.out_suffix}_summary", "lcc_trend_grouped.png")
    out_eu    = os.path.join(sample_dir, f"viz_results_{args.out_suffix}_euchr", "summary", "lcc_trend_euchr.png")
    out_h3    = os.path.join(sample_dir, f"viz_results_{args.out_suffix}_h3k4",  "summary", "lcc_trend_h3k4.png")

    # 画图：一张分组对比 + 各自单张
    plot_grouped_bars(used_thr, eu_vals, h3_vals, out_group, color_eu=args.color_euchr, color_h3=args.color_h3k4)
    plot_single_bars(used_thr, eu_vals, out_eu, "Largest CC ratio (euchr)", color=args.color_euchr)
    plot_single_bars(used_thr, h3_vals, out_h3, "Largest CC ratio (h3k4me3)", color=args.color_h3k4)

    # 存一份 CSV 方便你核查
    df = pd.DataFrame({
        "threshold": used_thr,
        "ratio_euchr": eu_vals,
        "ratio_h3k4": h3_vals
    })
    csv_path = os.path.join(sample_dir, f"viz_results_{args.out_suffix}_summary", "lcc_trend_values.tsv")
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    df.to_csv(csv_path, sep="\t", index=False)

    # 汇报存储结构
    print("\n=== SUMMARY (存储结构) ===")
    print(f"Sample dir   : {sample_dir}")
    print(f"Metrics in   :")
    print(f"  euchr      -> {os.path.join(sample_dir, 'graph_matrix_dual_euchr', 'whole{thr}', f'{sample_name}_whole{{thr}}_metrics.txt')}")
    print(f"  h3k4       -> {os.path.join(sample_dir, 'graph_matrix_dual_h3k4',  'whole{thr}', f'{sample_name}_whole{{thr}}_metrics.txt')}")
    print(f"Outputs to   :")
    print(f"  grouped    -> {out_group}")
    print(f"  euchr only -> {out_eu}")
    print(f"  h3k4  only -> {out_h3}")
    print(f"  values.tsv -> {csv_path}")
    print("=========================")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[FATAL] {e}", file=sys.stderr)
        sys.exit(1)


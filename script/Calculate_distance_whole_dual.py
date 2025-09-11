#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, math, argparse, re
from pathlib import Path

try:
    from scipy.spatial import cKDTree  # 用于全局KD-tree
except Exception:
    cKDTree = None

# ---------- 基础：样本定位与格式判定 ----------
def ascend_to_sample_dir(start: Path) -> Path:
    p = start.resolve()
    if p.is_file(): p = p.parent
    while True:
        if p.name.startswith("Sample-"): return p
        if p.parent == p: raise FileNotFoundError(f"上溯失败：{start}")
        p = p.parent

def detect_format_from_file(first_line: str):
    """
    返回 ('euchr'|'h3k4', has_header: bool, locus_idx: int)
    - euchr: 第一列 homolog，如 chr1(mat)，第二列 locus → locus_idx=1
    - h3k4 : 第一列 chrom, 第二列 allele, 第三列 locus → locus_idx=2
    """
    s = first_line.strip()
    parts = s.split()
    if not parts: return "euchr", False, 1
    h0 = parts[0].lower()
    if h0 == "homolog": return "euchr", True, 1
    if h0 in ("chrom", "allele"): return "h3k4", True, 2
    if re.match(r"^chr[^()]+\((mat|pat)\)$", parts[0]): return "euchr", False, 1
    if len(parts) >= 2 and parts[1] in ("mat","pat"): return "h3k4", False, 2
    return "euchr", False, 1

def parse_line(parts, fmt):
    """
    统一返回 (chrom_label, locus, x, y, z)
    - euchr: chrom_label = homolog (chr1(mat))
    - h3k4 : chrom_label = f"{chrom}({allele})"
    """
    try:
        if fmt == "euchr":
            if len(parts) < 5: return None
            chrom_label = parts[0]
            locus = parts[1]
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            return chrom_label, locus, x, y, z
        else:
            if len(parts) < 6: return None
            chrom_label = f"{parts[0]}({parts[1]})"
            locus = parts[2]
            x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
            return chrom_label, locus, x, y, z
    except ValueError:
        return None

# ---------- KD-tree ----------
def compute_all_pairs_kdtree(coords, threshold):
    points = [c[2:5] for c in coords]
    tree = cKDTree(points)
    pairs = tree.query_pairs(r=threshold, output_type='set')
    dists, idxs = tree.query(points, k=2)
    nearest_info = [(dists[i][1], idxs[i][1]) for i in range(len(points))]
    return pairs, nearest_info

# ---------- 主流程 ----------
def main():
    ap = argparse.ArgumentParser(description="Whole-genome pairwise distances (dual formats)")
    ap.add_argument("input_file", help="*.euchromatin_cluster.txt 或 *.h3k4me3_cluster.txt")
    ap.add_argument("output_file", help="输出文件路径（由批脚本放到 dual_euchr/dual_h3k4）")
    ap.add_argument("--source_label", type=str, default=None, help="euchr/h3k4；若不传则自动判定")
    ap.add_argument("--range", type=int, nargs=2, metavar=('START','END'),
                    help="可选：仅处理 i∈[START,END] 的bin（分块计算用）")
    ap.add_argument("--threshold", type=float, default=5.0,
                    help="距离阈值（默认5.0）")
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    sample_dir = ascend_to_sample_dir(script_dir)
    sample_name = sample_dir.name

    # 读取与格式判定
    coords = []  # [(chrom_label, locus, x, y, z)]
    with open(args.input_file, 'r', encoding='utf-8') as fin:
        first = fin.readline()
        fmt, has_header, locus_idx = detect_format_from_file(first)
        if args.source_label:
            if args.source_label not in ("euchr","h3k4"):
                print("[ERROR] --source_label 只能是 euchr/h3k4", file=sys.stderr); sys.exit(2)
            fmt = args.source_label  # 以显式为准
        if not has_header:
            fin.seek(0)
        for line in fin:
            if not line.strip(): continue
            parts = line.split()
            parsed = parse_line(parts, fmt)
            if parsed is None: continue
            coords.append(parsed)

    n = len(coords)
    if n == 0:
        print(f"[Info] No data points in {args.input_file}")
        sys.exit(0)

    threshold = args.threshold
    use_range = args.range is not None
    start_idx, end_idx = (args.range if use_range else (0, n-1))
    if start_idx < 0: start_idx = 0
    if end_idx >= n: end_idx = n-1

    pairs_written = 0
    added_nn = 0

    with open(args.output_file, 'w', encoding='utf-8') as fout:
        if (not use_range) and cKDTree:
            print(f"[Info] KD-tree all-pairs ≤ {threshold}")
            pairs, nearest_info = compute_all_pairs_kdtree(coords, threshold)
            for i, j in pairs:
                c1, l1, x1, y1, z1 = coords[i]
                c2, l2, x2, y2, z2 = coords[j]
                dist = math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
                fout.write(f"{c1}\t{l1}\t{c2}\t{l2}\t{dist}\n")
                pairs_written += 1
            has_edge = [False]*n
            for i,j in pairs:
                has_edge[i]=True; has_edge[j]=True
            added_pairs=set()
            for i,(d,nb) in enumerate(nearest_info):
                if not has_edge[i] and nb is not None:
                    key = tuple(sorted((i,nb)))
                    if key in added_pairs: continue
                    added_pairs.add(key)
                    c1, l1, x1, y1, z1 = coords[i]
                    c2, l2, x2, y2, z2 = coords[nb]
                    dist = math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
                    fout.write(f"{c1}\t{l1}\t{c2}\t{l2}\t{dist}\n")
                    pairs_written += 1; added_nn += 1
        else:
            if use_range:
                print(f"[Info] Brute-force range i∈[{start_idx},{end_idx}] (thr={threshold})")
            else:
                print(f"[Info] cKDTree unavailable, brute-force all (thr={threshold})")
            nearest_dist = {}  # i -> (min_dist, j)
            for i in range(start_idx, end_idx+1):
                c1, l1, x1, y1, z1 = coords[i]
                min_dist = float('inf'); min_j = None
                for j in range(i+1, n):
                    _, l2, x2, y2, z2 = coords[j]
                    dist = math.hypot(math.hypot(x1-x2, y1-y2), z1-z2)
                    if dist <= threshold:
                        fout.write(f"{c1}\t{l1}\t{coords[j][0]}\t{l2}\t{dist}\n")
                        pairs_written += 1
                        nearest_dist[i] = (0, None)
                        nearest_dist[j] = (0, None)
                    if dist < min_dist:
                        min_dist = dist; min_j = j
                if i not in nearest_dist:
                    nearest_dist[i] = (min_dist, min_j)
            for j in range(end_idx+1, n):
                if j not in nearest_dist or nearest_dist[j][0] != 0:
                    min_dist = float('inf'); min_i = None
                    for i in range(start_idx, end_idx+1):
                        dx = coords[i][2] - coords[j][2]
                        dy = coords[i][3] - coords[j][3]
                        dz = coords[i][4] - coords[j][4]
                        dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                        if dist < min_dist:
                            min_dist = dist; min_i = i
                    if j not in nearest_dist or min_dist < nearest_dist[j][0]:
                        nearest_dist[j] = (min_dist, min_i)
            added_pairs=set()
            for i,(dist,nb) in nearest_dist.items():
                if dist == 0 or nb is None: continue
                key = tuple(sorted((i,nb)))
                if key in added_pairs: continue
                added_pairs.add(key)
                c1, l1 = coords[i][0], coords[i][1]
                c2, l2 = coords[nb][0], coords[nb][1]
                fout.write(f"{c1}\t{l1}\t{c2}\t{l2}\t{dist}\n")
                pairs_written += 1; added_nn += 1

    # -------- REPORT：输出存储结构 --------
    print("\n=== Calculate_distance_whole_dual REPORT ===")
    print(f"Sample name : {sample_name}")
    print(f"Script dir  : {script_dir}")
    print(f"Input file  : {Path(args.input_file).resolve()}")
    print(f"Format      : {fmt} (locus_col={locus_idx})")
    print(f"Mode        : {'KD-tree' if (not use_range and cKDTree) else ('range' if use_range else 'bruteforce')}")
    print(f"Threshold   : {threshold}")
    print(f"Bins        : {n}  | Range: [{start_idx},{end_idx}]")
    print(f"Output file : {Path(args.output_file).resolve()}")
    print(f"Pairs<=thr  : {pairs_written}  | Added NN: {added_nn}")
    print("===========================================\n")

if __name__ == "__main__":
    main()


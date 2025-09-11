#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, os, re, sys
from pathlib import Path
import networkx as nx

def ascend_to_sample_dir(start: Path) -> Path:
    p = start.resolve()
    if p.is_file(): p = p.parent
    while True:
        if p.name.startswith("Sample-"): return p
        if p.parent == p: raise FileNotFoundError(f"上溯失败：{start}")
        p = p.parent

def iter_edges(dist_file, thr):
    with open(dist_file, encoding='utf-8') as fin:
        for line in fin:
            p = line.strip().split()
            if len(p) < 5: continue
            try:
                d = float(p[4])
            except ValueError:
                continue
            if d <= thr:
                yield f"{p[0]}:{p[1]}", f"{p[2]}:{p[3]}"

def collect_all_endpoints(dist_file):
    nodes = set()
    with open(dist_file, encoding='utf-8') as fin:
        for line in fin:
            p = line.strip().split()
            if len(p) >= 4:
                nodes.add(f"{p[0]}:{p[1]}")
                nodes.add(f"{p[2]}:{p[3]}")
    return nodes

def collect_all_bins_from_cluster(cluster_file):
    """
    支持两种表头：
    - euchr:  homolog locus x y z ...
    - h3k4 :  chrom   allele locus x y z ...
    """
    nodes = set()
    with open(cluster_file, encoding='utf-8') as f:
        first = f.readline().strip()
        has_header = first and first.split()[0].lower() in ("homolog","chrom")
        if not has_header:
            f.seek(0)
        for line in f:
            if not line.strip(): continue
            p = line.strip().split()
            if len(p) < 2: continue
            head = first.split()[0].lower() if has_header else None
            if head == "homolog" or (not has_header and re.match(r"^chr[^()]+\((mat|pat)\)$", p[0])):
                # euchr
                nodes.add(f"{p[0]}:{p[1]}")
            else:
                # h3k4
                if len(p) < 3: continue
                nodes.add(f"{p[0]}({p[1]}):{p[2]}")
    return nodes

def build_and_analyze(dist_file, thr, node_mode, cluster_file=None):
    G = nx.Graph()

    if node_mode == "all_distance_endpoints":
        nodes = collect_all_endpoints(dist_file); G.add_nodes_from(nodes)
    elif node_mode == "all_bins":
        if not cluster_file:
            raise ValueError("--node_mode all_bins 需要 --cluster_file")
        nodes = collect_all_bins_from_cluster(cluster_file); G.add_nodes_from(nodes)
    # else: leq_thr_endpoints → 不预置节点

    for u, v in iter_edges(dist_file, thr):
        G.add_edge(u, v)

    comps = list(nx.connected_components(G))
    mapping = {node: cid+1 for cid, comp in enumerate(comps) for node in comp}
    stats = {
        'threshold': thr,
        'node_mode': node_mode,
        'num_nodes': G.number_of_nodes(),
        'num_edges': G.number_of_edges(),
        'num_components': len(comps),
        'largest_cc_size': max((len(c) for c in comps), default=0),
        'avg_clustering': nx.average_clustering(G) if G.number_of_nodes() else 0.0,
        'density': nx.density(G) if G.number_of_nodes() > 1 else 0.0,
    }
    return mapping, stats

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--distance_file', required=True)
    ap.add_argument('--threshold', type=float, required=True)
    ap.add_argument('--output_prefix', required=True, help='如 whole1 或 whole1.25')
    ap.add_argument('--node_mode', default='leq_thr_endpoints',
                   choices=['leq_thr_endpoints','all_distance_endpoints','all_bins'])
    ap.add_argument('--cluster_file', default=None, help='all_bins 模式需要')
    ap.add_argument('--source_label', required=True, choices=['euchr','h3k4'],
                   help='来源标签（决定输出落地目录）')
    ap.add_argument('--out_root', default=None,
                   help='覆盖输出根目录；默认在工程根下 graph_matrix_dual_{label}')
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent
    sample_dir = ascend_to_sample_dir(script_dir)
    sample_name = sample_dir.name

    # 输出根目录（按来源分开）
    out_root = Path(args.out_root) if args.out_root else (project_dir / f"graph_matrix_dual_{args.source_label}")
    out_root.mkdir(parents=True, exist_ok=True)

    base = os.path.basename(args.distance_file).replace('_distance.txt','').replace('_distance_filtered.txt','')
    mdir = out_root / args.output_prefix
    mdir.mkdir(parents=True, exist_ok=True)
    cdir = out_root / 'components_single'
    cdir.mkdir(parents=True, exist_ok=True)

    mapping, stats = build_and_analyze(args.distance_file, args.threshold, args.node_mode, args.cluster_file)

    # metrics
    mf = mdir / f"{base}_{args.output_prefix}_metrics.txt"
    with open(mf, 'w', encoding='utf-8') as f:
        for k, v in stats.items():
            f.write(f"{k}\t{v}\n")

    # components（两列：locus_id, component_<prefix>）
    compf = cdir / f"{base}_comp_{args.output_prefix}.txt"
    colname = f"component_{args.output_prefix}"
    with open(compf,'w', encoding='utf-8') as f:
        f.write(f"locus_id\t{colname}\n")
        for locus_id, cid in mapping.items():
            f.write(f"{locus_id}\t{cid}\n")

    # REPORT（输出存储结构）
    print("\n=== analyze_graph_parallel_dual REPORT ===")
    print(f"Sample name : {sample_name}")
    print(f"Distance in : {Path(args.distance_file).resolve()}")
    print(f"Source label: {args.source_label}")
    if args.node_mode == "all_bins":
        print(f"Cluster in  : {Path(args.cluster_file).resolve() if args.cluster_file else 'N/A'}")
    print(f"Out root    : {out_root.resolve()}")
    print(f"Metrics dir : {mdir.resolve()}")
    print(f"Comp single : {cdir.resolve()}")
    print(f"Files out   : {mf.name}, {compf.name}")
    print("==========================================\n")

    print(f"Finished threshold={args.threshold} [{args.node_mode}] -> {mf}, {compf}")


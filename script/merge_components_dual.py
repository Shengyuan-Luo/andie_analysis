#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, glob, sys, re
import pandas as pd
from pathlib import Path

def ascend_to_sample_dir(start: Path) -> Path:
    p = start.resolve()
    if p.is_file(): p = p.parent
    while True:
        if p.name.startswith("Sample-"): return p
        if p.parent == p: raise FileNotFoundError(f"上溯失败：{start}")
        p = p.parent

def read_and_concat_split(label_dir: Path, label: str) -> pd.DataFrame:
    """
    读取 Split_based_on_chr_dual_{label}/*.txt 并拼接成一个 DataFrame
    - euchr: 需要列 homolog, locus
    - h3k4 : 需要列 chrom, allele, locus  -> 制作 homolog_like = f"{chrom}({allele})"
    """
    files = sorted(label_dir.glob("*.txt"))
    if not files:
        raise FileNotFoundError(f"未发现拆分文件：{label_dir}")
    dfs = []
    for fp in files:
        with open(fp, 'r', encoding='utf-8') as f:
            first = f.readline().strip()
        header = first.split()
        fdf = pd.read_csv(fp, sep=r"\s+", engine="python")
        if header and header[0].lower() == "homolog":    # euchr
            need = {"homolog","locus"}
            if not need.issubset(fdf.columns):
                raise ValueError(f"{fp} 缺少列 {need}")
            fdf["_homolog_like"] = fdf["homolog"].astype(str)
            fdf["_locus"] = fdf["locus"].astype(str)
        else:                                            # h3k4
            need = {"chrom","allele","locus"}
            if not need.issubset(fdf.columns):
                raise ValueError(f"{fp} 缺少列 {need}")
            fdf["_homolog_like"] = fdf["chrom"].astype(str) + "(" + fdf["allele"].astype(str) + ")"
            fdf["_locus"] = fdf["locus"].astype(str)
        dfs.append(fdf[["_homolog_like","_locus"]])
    merged = pd.concat(dfs, ignore_index=True)
    merged.rename(columns={"_homolog_like":"homolog_like","_locus":"locus"}, inplace=True)
    merged["locus_id"] = merged["homolog_like"] + ":" + merged["locus"]
    return merged

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="合并多阈值组件（dual来源）")
    ap.add_argument("--source_label", required=True, choices=["euchr","h3k4"])
    args = ap.parse_args()

    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent
    sample_dir = ascend_to_sample_dir(script_dir)
    base_name = sample_dir.name

    # 输入/输出目录（按来源分开）
    split_dir = project_dir / f"Split_based_on_chr_dual_{args.source_label}"
    comp_single_dir = project_dir / f"graph_matrix_dual_{args.source_label}" / "components_single"
    out_dir = project_dir / f"graph_matrix_dual_{args.source_label}" / "components"
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) 原始 bins（由多个拆分文件拼接）
    orig_df = read_and_concat_split(split_dir, args.source_label)

    # 2) 组件文件（多个阈值）
    pattern = str(comp_single_dir / f"{base_name}_comp_whole*.txt")
    comp_files = sorted(glob.glob(pattern))
    if not comp_files:
        raise FileNotFoundError(f"在 {comp_single_dir} 下未找到 {pattern} 匹配文件")

    merged = orig_df.copy()
    for comp_file in comp_files:
        df = pd.read_csv(comp_file, sep=r"\s+", engine="python")
        if "locus_id" not in df.columns or df.shape[1] != 2:
            raise ValueError(f"{comp_file} 应为两列：locus_id 与 component_wholeX")
        merged = merged.merge(df, on="locus_id", how="left")

    out_path = out_dir / f"{base_name}_components.txt"
    # 输出列顺序：homolog_like, locus, locus_id, component_whole*
    # （可按需调整；保留 homolog_like 便于核查）
    merged.to_csv(out_path, sep="\t", index=False)

    # REPORT
    print("\n=== merge_components_dual REPORT ===")
    print(f"Sample name    : {base_name}")
    print(f"Source label   : {args.source_label}")
    print(f"Split dir in   : {split_dir.resolve()}")
    print(f"Comp single in : {comp_single_dir.resolve()}")
    print(f"Output dir     : {out_dir.resolve()}")
    print(f"Output file    : {out_path.resolve()}")
    print("====================================\n")

    print(f"✅ 合并完成 → {out_path}")


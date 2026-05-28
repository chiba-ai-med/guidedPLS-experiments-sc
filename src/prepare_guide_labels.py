#!/usr/bin/env python3
# =============================================================================
# guided-PLS 入力データ準備: 発現行列 + ガイドラベル (one-hot)
# =============================================================================
# Usage: python3 src/prepare_guide_labels.py \
#   --rna-rds <rna_combined.rds> --rna-meta <rna_metadata.csv> \
#   --atac-rds <atac_processed.rds> --atac-meta <atac_metadata.csv> \
#   --condition <none|stage|germlayer|stage_germlayer> \
#   --outdir <output_dir> --seed <42>
# =============================================================================

import argparse
import numpy as np
import pandas as pd
import subprocess
import os
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna-rds", required=True)
    parser.add_argument("--rna-meta", required=True)
    parser.add_argument("--atac-rds", required=True)
    parser.add_argument("--atac-meta", required=True)
    parser.add_argument("--condition", required=True,
                        choices=["none", "stage", "germlayer", "stage_germlayer",
                                 "stage_x_germlayer", "stage_binned"])
    parser.add_argument("--dataset", default="Organogenesis",
                        help="Dataset name (Organogenesis, Testis, ...)")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    np.random.seed(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    # --- 1. メタデータ読み込み ---
    rna_meta = pd.read_csv(args.rna_meta)
    atac_meta = pd.read_csv(args.atac_meta)

    print(f"RNA cells: {len(rna_meta)}, ATAC cells: {len(atac_meta)}")

    # --- 2. 特徴量行列をRから抽出 ---
    # guided-PLS はクロスモーダル: X1(RNA genes) と X2(ATAC peaks) は異なる特徴量
    r_script = f"""
    suppressPackageStartupMessages({{library(Seurat); library(Signac)}})

    rna <- readRDS("{args.rna_rds}")
    atac <- readRDS("{args.atac_rds}")

    # Seurat v5: JoinLayers if needed
    if (inherits(rna[["RNA"]], "Assay5")) rna[["RNA"]] <- JoinLayers(rna[["RNA"]])

    # --- X1: RNA 遺伝子発現 (HVG) ---
    hvg <- VariableFeatures(rna)
    cat(sprintf("RNA HVGs: %d\\n", length(hvg)))
    x1 <- as.matrix(LayerData(rna, layer = "data")[hvg, ])

    # --- X2: ATAC ピーク行列 (TF-IDF正規化) ---
    DefaultAssay(atac) <- "ATAC"
    if (inherits(atac[["ATAC"]], "Assay5")) atac[["ATAC"]] <- JoinLayers(atac[["ATAC"]])

    # TF-IDF正規化
    atac <- RunTFIDF(atac)

    # 可変ピークを選択 (top features by variance)
    atac <- FindTopFeatures(atac, min.cutoff = "q75")  # 上位25%のピーク
    top_peaks <- VariableFeatures(atac)
    cat(sprintf("ATAC peaks (top features): %d\\n", length(top_peaks)))

    # TF-IDF正規化済みデータを抽出
    x2 <- as.matrix(LayerData(atac, layer = "data")[top_peaks, ])

    # cells × features に転置して保存
    write.csv(t(x1), "{args.outdir}/X1.csv", row.names = FALSE)
    write.csv(t(x2), "{args.outdir}/X2.csv", row.names = FALSE)
    cat(sprintf("X1: %d cells x %d genes\\n", ncol(x1), nrow(x1)))
    cat(sprintf("X2: %d cells x %d peaks\\n", ncol(x2), nrow(x2)))
    cat("Matrices saved.\\n")
    """

    print("Extracting expression matrices via R ...")
    result = subprocess.run(
        ["Rscript", "-e", r_script],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print("R stderr:", result.stderr, file=sys.stderr)
        sys.exit(1)
    print(result.stdout)

    # --- 3. ガイドラベル生成 ---
    print(f"Generating guide labels: condition = {args.condition}")

    if args.condition == "none":
        # 全細胞に同一ラベル (定数列 = unguided)
        y1 = pd.DataFrame({"const": np.ones(len(rna_meta), dtype=int)})
        y2 = pd.DataFrame({"const": np.ones(len(atac_meta), dtype=int)})

    elif args.condition == "stage":
        y1, y2 = build_onehot(rna_meta, atac_meta, ["Stage"])

    elif args.condition == "germlayer":
        # RNA: germ_layer列, ATAC: Germ.layer列
        rna_meta = rna_meta.copy()
        atac_meta = atac_meta.copy()
        rna_meta["germlayer_unified"] = rna_meta["germ_layer"]
        atac_meta["germlayer_unified"] = atac_meta["Germ.layer"]
        y1, y2 = build_onehot(rna_meta, atac_meta, ["germlayer_unified"])

    elif args.condition == "stage_germlayer":
        rna_meta = rna_meta.copy()
        atac_meta = atac_meta.copy()
        rna_meta["germlayer_unified"] = rna_meta["germ_layer"]
        atac_meta["germlayer_unified"] = atac_meta["Germ.layer"]
        y1, y2 = build_onehot(rna_meta, atac_meta,
                               ["Stage", "germlayer_unified"])

    elif args.condition == "stage_x_germlayer":
        # Interaction encoding: each (Stage, germ_layer) combination becomes a
        # single categorical level, raising the Y rank from 3+6=9 (additive
        # "stage_germlayer") up to |Stage| × |germ_layer| ≤ 18.
        rna_meta = rna_meta.copy()
        atac_meta = atac_meta.copy()
        rna_meta["stage_x_gl"] = (rna_meta["Stage"].astype(str)
                                  + "::" + rna_meta["germ_layer"].astype(str))
        atac_meta["stage_x_gl"] = (atac_meta["Stage"].astype(str)
                                   + "::" + atac_meta["Germ.layer"].astype(str))
        y1, y2 = build_onehot(rna_meta, atac_meta, ["stage_x_gl"])

    elif args.condition == "stage_binned":
        # Testis用: Stage列はビニング済み (E18, P2-3, P6-7)
        y1, y2 = build_onehot(rna_meta, atac_meta, ["Stage"])

    print(f"Y1 shape: {y1.shape}, Y2 shape: {y2.shape}")

    # --- 4. 保存 ---
    y1.to_csv(os.path.join(args.outdir, "Y1.csv"), index=False)
    y2.to_csv(os.path.join(args.outdir, "Y2.csv"), index=False)

    # メタデータもコピー
    rna_meta.to_csv(os.path.join(args.outdir, "rna_metadata.csv"), index=False)
    atac_meta.to_csv(os.path.join(args.outdir, "atac_metadata.csv"), index=False)

    print("Done.")


def build_onehot(rna_meta, atac_meta, columns):
    """
    指定列のone-hotエンコーディングを生成。
    ソースとターゲットで同じカラム集合を使用。
    """
    if len(columns) == 1:
        rna_labels = rna_meta[columns[0]].astype(str)
        atac_labels = atac_meta[columns[0]].astype(str)
    else:
        # 複数列: 各列ごとに独立にone-hotを作成し連結
        y1_parts = []
        y2_parts = []
        for col in columns:
            rna_s = rna_meta[col].astype(str)
            atac_s = atac_meta[col].astype(str)
            all_labels = sorted(
                pd.Index(rna_s.unique()).union(pd.Index(atac_s.unique()))
            )
            prefix = col + "_"
            ohe_rna = pd.get_dummies(rna_s).reindex(
                columns=all_labels, fill_value=0
            ).astype(int)
            ohe_rna.columns = [prefix + c for c in ohe_rna.columns]

            ohe_atac = pd.get_dummies(atac_s).reindex(
                columns=all_labels, fill_value=0
            ).astype(int)
            ohe_atac.columns = [prefix + c for c in ohe_atac.columns]

            y1_parts.append(ohe_rna.reset_index(drop=True))
            y2_parts.append(ohe_atac.reset_index(drop=True))

        return pd.concat(y1_parts, axis=1), pd.concat(y2_parts, axis=1)

    # 単一列のone-hot
    all_labels = sorted(
        pd.Index(rna_labels.unique()).union(pd.Index(atac_labels.unique()))
    )
    y1 = pd.get_dummies(rna_labels).reindex(
        columns=all_labels, fill_value=0
    ).astype(int)
    y2 = pd.get_dummies(atac_labels).reindex(
        columns=all_labels, fill_value=0
    ).astype(int)

    return y1.reset_index(drop=True), y2.reset_index(drop=True)


if __name__ == "__main__":
    main()

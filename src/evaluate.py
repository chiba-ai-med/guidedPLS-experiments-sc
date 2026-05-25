#!/usr/bin/env python3
# =============================================================================
# 評価: Accuracy, Balanced Accuracy, ARI, NMI, per-class F1
# =============================================================================
# Usage: python3 src/evaluate.py \
#   --atac-meta <atac_metadata.csv> \
#   --guidedpls-dir <output/guidedpls> \
#   --comparison-dir <output/comparison> \
#   --conditions none,stage,germlayer,stage_germlayer \
#   --methods seurat,harmony,liger \
#   --out-metrics <metrics.csv> \
#   --out-per-class <per_class_f1.csv>
# =============================================================================

import argparse
import os
import pandas as pd
import numpy as np
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    adjusted_rand_score,
    normalized_mutual_info_score,
    f1_score,
    classification_report,
)


def load_predictions(path):
    """predicted_labels.csv を読み込む"""
    df = pd.read_csv(path)
    return df["predicted_celltype"].values


def compute_metrics(y_true, y_pred):
    """全体指標を計算"""
    # 共通ラベルのみで評価（予測にないラベルも考慮）
    mask = pd.notna(y_true) & pd.notna(y_pred)
    y_true = y_true[mask]
    y_pred = y_pred[mask]

    return {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "ARI": adjusted_rand_score(y_true, y_pred),
        "NMI": normalized_mutual_info_score(y_true, y_pred),
        "macro_f1": f1_score(y_true, y_pred, average="macro", zero_division=0),
        "weighted_f1": f1_score(y_true, y_pred, average="weighted", zero_division=0),
        "n_cells": len(y_true),
    }


def compute_per_class_f1(y_true, y_pred):
    """クラスごとのF1スコア"""
    mask = pd.notna(y_true) & pd.notna(y_pred)
    y_true = y_true[mask]
    y_pred = y_pred[mask]

    labels = sorted(set(y_true))
    f1s = f1_score(y_true, y_pred, labels=labels, average=None, zero_division=0)
    return dict(zip(labels, f1s))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--atac-meta", required=True)
    parser.add_argument("--guidedpls-dir", required=True)
    parser.add_argument("--comparison-dir", required=True)
    parser.add_argument("--conditions", required=True)
    parser.add_argument("--methods", required=True)
    parser.add_argument("--out-metrics", required=True)
    parser.add_argument("--out-per-class", required=True)
    args = parser.parse_args()

    conditions = args.conditions.split(",")
    methods = [m for m in args.methods.split(",") if m and m != "NONE"]

    # --- Ground truth ---
    atac_meta = pd.read_csv(args.atac_meta)
    y_true = atac_meta["celltype"].values
    print(f"Ground truth: {len(y_true)} cells, {len(set(y_true))} cell types")

    # --- 各手法の評価 ---
    all_metrics = []
    all_per_class = []

    # guided-PLS 条件
    for cond in conditions:
        label = f"guidedpls_{cond}"
        pred_file = os.path.join(args.guidedpls_dir, cond, "predicted_labels.csv")
        y_pred = load_predictions(pred_file)

        metrics = compute_metrics(y_true, y_pred)
        metrics["method"] = label
        all_metrics.append(metrics)

        per_class = compute_per_class_f1(y_true, y_pred)
        for ct, f1 in per_class.items():
            all_per_class.append({
                "method": label, "celltype": ct, "f1": f1
            })

        print(f"{label}: acc={metrics['accuracy']:.3f}, "
              f"bal_acc={metrics['balanced_accuracy']:.3f}, "
              f"ARI={metrics['ARI']:.3f}, NMI={metrics['NMI']:.3f}")

    # 比較手法
    for method in methods:
        pred_file = os.path.join(args.comparison_dir, method, "predicted_labels.csv")
        y_pred = load_predictions(pred_file)

        metrics = compute_metrics(y_true, y_pred)
        metrics["method"] = method
        all_metrics.append(metrics)

        per_class = compute_per_class_f1(y_true, y_pred)
        for ct, f1 in per_class.items():
            all_per_class.append({
                "method": method, "celltype": ct, "f1": f1
            })

        print(f"{method}: acc={metrics['accuracy']:.3f}, "
              f"bal_acc={metrics['balanced_accuracy']:.3f}, "
              f"ARI={metrics['ARI']:.3f}, NMI={metrics['NMI']:.3f}")

    # --- 保存 ---
    metrics_df = pd.DataFrame(all_metrics)
    # method列を先頭に
    cols = ["method"] + [c for c in metrics_df.columns if c != "method"]
    metrics_df = metrics_df[cols]
    metrics_df.to_csv(args.out_metrics, index=False)
    print(f"\nMetrics saved to {args.out_metrics}")

    per_class_df = pd.DataFrame(all_per_class)
    per_class_df.to_csv(args.out_per_class, index=False)
    print(f"Per-class F1 saved to {args.out_per_class}")


if __name__ == "__main__":
    main()

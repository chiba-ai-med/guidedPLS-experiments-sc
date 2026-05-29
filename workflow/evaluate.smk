# =============================================================================
# 評価 + 可視化 — 全データセット共通
# =============================================================================

def eval_inputs(wc):
    """データセットごとのguided-PLS条件 + 比較手法の予測ラベルを列挙"""
    ds = wc.dataset
    inputs = {
        "atac_meta": f"output/{ds}/preprocess/atac_metadata.csv",
    }
    conds = guide_conditions_for(ds)
    for i, c in enumerate(conds):
        inputs[f"gp_{i}"] = f"output/{ds}/guidedpls/{c}/predicted_labels.csv"
    for i, m in enumerate(comparison_methods_for(ds)):
        inputs[f"cm_{i}"] = f"output/{ds}/comparison/{m}/predicted_labels.csv"
    return inputs

rule evaluate:
    input:
        unpack(eval_inputs),
    output:
        metrics="output/{dataset}/evaluation/metrics.csv",
        per_class="output/{dataset}/evaluation/per_class_f1.csv",
    params:
        guidedpls_dir=lambda wc: f"output/{wc.dataset}/guidedpls",
        comparison_dir=lambda wc: f"output/{wc.dataset}/comparison",
        conditions=lambda wc: ",".join(guide_conditions_for(wc.dataset)),
        methods=lambda wc: ",".join(comparison_methods_for(wc.dataset)) or "NONE",
    shell:
        """
        mkdir -p output/{wildcards.dataset}/evaluation
        python3 src/evaluate.py \
            --atac-meta {input.atac_meta} \
            --guidedpls-dir {params.guidedpls_dir} \
            --comparison-dir {params.comparison_dir} \
            --conditions {params.conditions} \
            --methods {params.methods} \
            --out-metrics {output.metrics} \
            --out-per-class {output.per_class}
        """

rule visualize:
    input:
        metrics="output/{dataset}/evaluation/metrics.csv",
        per_class="output/{dataset}/evaluation/per_class_f1.csv",
    output:
        barplot_acc="output/{dataset}/figures/barplot_accuracy.png",
    params:
        guidedpls_dir=lambda wc: f"output/{wc.dataset}/guidedpls",
        comparison_dir=lambda wc: f"output/{wc.dataset}/comparison",
        atac_meta=lambda wc: f"output/{wc.dataset}/preprocess/atac_metadata.csv",
        conditions=lambda wc: ",".join(guide_conditions_for(wc.dataset)),
        methods=lambda wc: ",".join(comparison_methods_for(wc.dataset)) or "NONE",
    shell:
        """
        mkdir -p output/{wildcards.dataset}/figures
        Rscript src/plot_results.R \
            {input.metrics} \
            {input.per_class} \
            {params.guidedpls_dir} \
            {params.comparison_dir} \
            {params.atac_meta} \
            {params.conditions} \
            {params.methods} \
            output/{wildcards.dataset}/figures
        """

# --- 計算時間 + メモリ使用量 ---
def resources_inputs(wc):
    ds = wc.dataset
    inputs = {}
    for c in guide_conditions_for(ds):
        inputs[f"gp_{c}"] = f"output/{ds}/benchmark/guidedpls_{c}.tsv"
    for m in comparison_methods_for(ds):
        inputs[f"cm_{m}"] = f"output/{ds}/benchmark/{m}.tsv"
    return inputs

rule resource_barplot:
    input:
        unpack(resources_inputs),
    output:
        png="output/{dataset}/figures/method_resources.png",
    params:
        bench_dir=lambda wc: f"output/{wc.dataset}/benchmark",
        conditions=lambda wc: ",".join(guide_conditions_for(wc.dataset)),
        methods=lambda wc: ",".join(comparison_methods_for(wc.dataset)) or "NONE",
    shell:
        """
        Rscript src/plot_method_resources.R \
            {params.bench_dir} {params.conditions} {params.methods} \
            {wildcards.dataset} {output.png}
        """

# --- kNN accuracy + Time + Memory の統合 3-panel (Dark2, 透明背景) ---
rule method_summary:
    input:
        unpack(resources_inputs),
        metrics="output/{dataset}/evaluation/metrics.csv",
    output:
        acc="output/{dataset}/figures/method_accuracy.png",
        time="output/{dataset}/figures/method_time.png",
        mem="output/{dataset}/figures/method_memory.png",
    params:
        bench_dir=lambda wc: f"output/{wc.dataset}/benchmark",
        outdir=lambda wc: f"output/{wc.dataset}/figures",
        conditions=lambda wc: ",".join(guide_conditions_for(wc.dataset)),
        methods=lambda wc: ",".join(comparison_methods_for(wc.dataset)) or "NONE",
    shell:
        """
        Rscript src/plot_method_summary.R \
            {input.metrics} {params.bench_dir} \
            {params.conditions} {params.methods} \
            {wildcards.dataset} {params.outdir}
        """

# --- 手法-色対応の横向き legend (別ファイル、透明背景) ---
rule method_legend:
    output:
        png="output/{dataset}/figures/method_legend.png",
    params:
        conditions=lambda wc: ",".join(guide_conditions_for(wc.dataset)),
        methods=lambda wc: ",".join(comparison_methods_for(wc.dataset)) or "NONE",
    shell:
        """
        Rscript src/plot_method_legend.R \
            {params.conditions} {params.methods} \
            {wildcards.dataset} {output.png}
        """

# --- Per-class F1 ヒートマップ ---
rule per_class_f1_heatmap:
    input:
        per_class="output/{dataset}/evaluation/per_class_f1.csv",
    output:
        png="output/{dataset}/figures/per_class_f1_heatmap.png",
    params:
        conditions=lambda wc: ",".join(guide_conditions_for(wc.dataset)),
        methods=lambda wc: ",".join(comparison_methods_for(wc.dataset)) or "NONE",
    shell:
        """
        Rscript src/plot_per_class_f1_heatmap.R \
            {input.per_class} {params.conditions} {params.methods} {output.png}
        """

# --- gPLS UMAP (条件ごと) ---
rule umap_gpls:
    input:
        rdata="output/{dataset}/guidedpls/{condition}/guidedpls.RData",
        rna_meta="output/{dataset}/preprocess/rna_metadata.csv",
        atac_meta="output/{dataset}/preprocess/atac_metadata.csv",
        preds="output/{dataset}/guidedpls/{condition}/predicted_labels.csv",
    output:
        png="output/{dataset}/figures/umap_gpls/{condition}/umap_combined.png",
    params:
        outdir=lambda wc: f"output/{wc.dataset}/figures/umap_gpls/{wc.condition}",
    shell:
        """
        Rscript src/plot_gpls_umap.R \
            {input.rdata} {input.rna_meta} {input.atac_meta} {input.preds} \
            {wildcards.dataset} {wildcards.condition} {params.outdir}
        """

# --- ラベル流れ (alluvial) ---
rule label_flow_gpls:
    input:
        preds="output/{dataset}/guidedpls/{condition}/predicted_labels.csv",
        atac_meta="output/{dataset}/preprocess/atac_metadata.csv",
    output:
        png="output/{dataset}/figures/label_flow/guidedpls_{condition}.png",
    shell:
        """
        Rscript src/plot_label_flow.R \
            {input.preds} {input.atac_meta} \
            "gPLS ({wildcards.condition})" {output.png}
        """

rule label_flow_comparison:
    input:
        preds="output/{dataset}/comparison/{method}/predicted_labels.csv",
        atac_meta="output/{dataset}/preprocess/atac_metadata.csv",
    output:
        png="output/{dataset}/figures/label_flow/{method}.png",
    wildcard_constraints:
        method="seurat|harmony|scanorama"
    shell:
        """
        Rscript src/plot_label_flow.R \
            {input.preds} {input.atac_meta} \
            "{wildcards.method}" {output.png}
        """

# --- 手法間 UMAP 比較 ---
def umap_comparison_inputs(wc):
    ds = wc.dataset
    inputs = {}
    for c in guide_conditions_for(ds):
        inputs[f"gp_{c}"] = f"output/{ds}/guidedpls/{c}/embeddings.csv"
    for m in comparison_methods_for(ds):
        inputs[f"cm_{m}"] = f"output/{ds}/comparison/{m}/embeddings.csv"
    return inputs

def umap_comparison_pairs(wc):
    ds = wc.dataset
    pairs = []
    for c in guide_conditions_for(ds):
        pairs.append(f"gPLS_{c}=output/{ds}/guidedpls/{c}/embeddings.csv")
    for m in comparison_methods_for(ds):
        pairs.append(f"{m}=output/{ds}/comparison/{m}/embeddings.csv")
    return " ".join(pairs)

rule umap_method_comparison:
    input:
        unpack(umap_comparison_inputs),
        rna_meta="output/{dataset}/preprocess/rna_metadata.csv",
        atac_meta="output/{dataset}/preprocess/atac_metadata.csv",
    output:
        png="output/{dataset}/figures/umap_methods_combined.png",
        leg_mod="output/{dataset}/figures/umap_legend_modality.png",
        leg_ct ="output/{dataset}/figures/umap_legend_celltype.png",
    params:
        pairs=umap_comparison_pairs,
        outdir=lambda wc: f"output/{wc.dataset}/figures",
    shell:
        """
        Rscript src/plot_method_umap_comparison.R \
            --meta={input.rna_meta},{input.atac_meta} \
            {params.pairs} {wildcards.dataset} {params.outdir}
        """

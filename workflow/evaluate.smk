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

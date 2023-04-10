
rule feature_overlaps:
    input:
        genes = refs("results/flybase-anno/transcriptome.gtf.gz"),
        ins = rules.collect_all_tidal_insertions.output
    output:
        tsv="results/overlaps/overlaps.tsv.gz"
    params:
        maxgap =  config.get("MAXGAP_FOR_TE_OVERLAP")
    resources:
        time=20,
        mem=12000,
        cpus=1
    script:
        "../scripts/per_strain_te_gene_overlaps.R"

rule excluded_gene_te_pairs:
    """
    list of gene/TE pairs that are excluded from the analysis
    because they are fixed in tidal, making it impossible to estimate coefficients
    """
    input:
        tsv = rules.feature_overlaps.output.tsv
    output:
        tsv="results/overlaps/excluded_gene_te_pairs.tsv"
    resources:
        time=20,
        mem=12000,
        cpus=1
    script:
        "../scripts/excluded_gene_te_pairs.R"

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

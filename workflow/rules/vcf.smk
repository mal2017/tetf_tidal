rule bed_to_var:
    """
    processes tidal output into a format that vep can be easily turned into a vcf
    with VariantAnnotation pkg.
    """
    input:
        poly_dir = rules.extract_polymorphic_insertions.output
    params:
        poly_ins_tsv = rules.extract_polymorphic_insertions.output.dir + "/{tidal_strain}_result/{tidal_strain}_Inserts_Annotated.txt",
        poly_deps_tsv = rules.extract_polymorphic_insertions.output.dir + "/{tidal_strain}_result/{tidal_strain}_Depletion_Annotated_TEonly.txt",
        sample_name = "{tidal_strain}",
    output:
        rds = "results/vcf/{tidal_group}/{tidal_strain}_polymorphic.rds",
    script:
        "../scripts/per_strain_vars.R"

def aggregate_variants_within_group(wildcards):
    checkpoint_output = checkpoints.extract_polymorphic_insertions.get(**wildcards).output.dir
    #print(checkpoint_output)
    wc_path = os.path.join(checkpoint_output, "{d}_result/")
    #print(wc_path)
    x = glob_wildcards(wc_path)
    #print(x)
    return expand("results/vcf/{tidal_group}/{i}_polymorphic.rds",
           tidal_group=wildcards.tidal_group,
           i=x.d)


rule collect_tidal_vcf:
    input:
        aggregate_variants_within_group
    params:
        vcf = "results/vcf/{tidal_group}_polymorphic_tes.vcf"
    output:
        vcf = "results/vcf/{tidal_group}_polymorphic_tes.vcf.bgz"
    script:
        "../scripts/collect_as_vcf.R"

"""
    vep accepts a specific tab-delimited format for input
    see: https://useast.ensembl.org/info/docs/tools/vep/vep_formats.html#default

    specifically: chr start end genotype strand identifier

    identifier is a unique identifier for the variant

    genotype should be INS for insertions, DEL for deletions
"""
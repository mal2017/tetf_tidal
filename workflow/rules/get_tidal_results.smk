rule get_tidal_zip:
    params:
        tidal = config.get("TIDAL_ZIP")
    output:
        temp("results/TIDAL.zip")
    shell:
        """
        wget {params.tidal} -O {output}
        """

rule extract_polymorphic_insertions_0:
    input:
        rules.get_tidal_zip.output
    output:
        dgn = temp("results/polymorphic_insertions/DGN_flies.zip"),
        labstrain = temp("results/polymorphic_insertions/LabStrain_flies.zip"),
        dgrp = temp("results/polymorphic_insertions/DGRP_flies.zip"),
        lines = temp("results/polymorphic_insertions/CellLines.zip"),
        pools = temp("results/polymorphic_insertions/Pool_Flies.zip"),
    shell:
        """
        unzip {input} -d results/polymorphic_insertions/
        """

TIDAL_STANDARD_OUTPUTS = ["_Depletion_Annotated.bed","_Depletion_Annotated_TEonly.bed",
"_Depletion_Annotated_TEonly.txt","_Depletion_Annotated_TEonly.txt","_fixed_bin.txt",
"_Inserts_Annotated.bed","_Inserts_Annotated.txt","_map_insertion_depletion.txt",
"_ReadDepletion.txt","_ReadInsertion.txt","_summary.txt","_TE_Indel_genome_plot.pdf"]


checkpoint extract_polymorphic_insertions:
    input:
        "results/polymorphic_insertions/{tidal_group}.zip",
    output:
        dir = directory("results/polymorphic_insertions/dl/{tidal_group}/")
    shell:
        """
        mkdir -p {output.dir} &&
        unzip {input} -d {output}
        """

rule get_per_strain_insertions:
    input:
        ref_copies = refs("results/repeatmasker/reference_insertions.bed"),
        poly_dir = rules.extract_polymorphic_insertions.output
    params:
        poly_ins_tsv = rules.extract_polymorphic_insertions.output.dir + "/{tidal_strain}_result/{tidal_strain}_Inserts_Annotated.txt",
        poly_deps_tsv = rules.extract_polymorphic_insertions.output.dir + "/{tidal_strain}_result/{tidal_strain}_Depletion_Annotated_TEonly.txt",
    output:
        tsv="results/polymorphic_insertions/{tidal_group}/{tidal_strain}.tsv"
    script:
        '../scripts/per_strain_ins.R'

def aggregate_known_insertions_within_group(wildcards):
    checkpoint_output = checkpoints.extract_polymorphic_insertions.get(**wildcards).output.dir
    #print(checkpoint_output)
    wc_path = os.path.join(checkpoint_output, "{d}_result/")
    #print(wc_path)
    x = glob_wildcards(wc_path)
    #print(x)
    return expand("results/polymorphic_insertions/{tidal_group}/{i}.tsv",
           tidal_group=wildcards.tidal_group,
           i=x.d)

localrules: collect_all_tidal_insertions

rule collect_known_insertions_within_group:
    input:
        aggregate_known_insertions_within_group
    output:
        tsv="results/collected_insertions/{tidal_group}.tsv"
    resources:
        time=60,
        mem=24000,
        cpus=1
    script:
        "../scripts/collect_known_insertions_within_group.R"

rule collect_all_tidal_insertions:
    input:
        expand("results/collected_insertions/{g}.tsv",g=["DGRP_flies"])
    output:
        "results/collected_insertions/insertions_by_strain.tsv.gz"
    resources:
        time=60,
        mem=24000,
        cpus=1
    shell:
        "xsv cat rows -d '\t' {input} | tr ',' '\t' | gzip -c > {output}"

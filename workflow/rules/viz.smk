rule annotate_insertion_penetrance:
    """
    gets all of the reconstituted insertion data into a bigbed for easy sanity checking
    and access without boilerplate importing stuff in downstream scripts
    """
    input:
        genome_fa = refs("results/flybase-anno/genome.fasta.gz"),
        ins = rules.collect_all_tidal_insertions.output,
    output:
        all_ins = "results/beds/dgrp_tidal_insertions.bb",
        penetrance = "results/beds/dgrp_tidal_insertions.unique.bb"
    resources:
        time=20,
        mem=32000,
        cpus=1
    script:
        "../scripts/annotate_insertion_penetrance.R"
rule arriba:
    input:
        aln = rules.star_aln.output,
        genome = config.get("GENOME"),
        gtf = config.get("GTF")
    output:
        fus = "results/arriba/{sample}/fusions.tsv",
        discarded = "results/arriba/{sample}/fusions-discarded.tsv"
    singularity:
        "docker://quay.io/biocontainers/arriba:2.1.0--h3198e80_1"
    resources:
        time=480,
        mem=64000,
        cpus=24
    shell:
        """
        arriba -x {input.aln}/Aligned.sortedByCoord.out.bam \
            -o {output.fus} -O {output.discarded} \
            -a {input.genome} -g {input.gtf} -f blacklist
        """

rule arriba:
    input:
        aln = rules.star_aln.output,
        genome = custom_genome('results/custom-genome/combined.fasta'),
        gtf = custom_genome('results/custom-genome/combined.fixed.gtf')
    output:
        fus = "results/arriba/{sample}/fusions.tsv",
        discarded = "results/arriba/{sample}/fusions-discarded.tsv"
    singularity:
        "docker://quay.io/biocontainers/arriba:2.1.0--h3198e80_1"
    resources:
        time=60,
        mem=20000,
        cpus=2
    shell:
        """
	VIRAL=$(grep ">" {input.genome} | grep -v "type" | sed 's/>\+//' | tr "\n" ",")
        arriba -x {input.aln}/Aligned.sortedByCoord.out.bam \
            -v $VIRAL \
            -o {output.fus} -O {output.discarded} \
            -a {input.genome} -g {input.gtf} -U 32766 -f blacklist,uninteresting_contigs,read_through,intronic,long_gap,intragenic_exonic,no_coverage,mismatches,duplicates,top_expressed_viral_contigs,low_coverage_viral_contigs
        """

rule arriba_y_fusions_filt:
    input:
        rules.arriba.output.fus
    output:
        "results/arriba/{sample}/fusions.y.tsv"
    shell:
        """
        head -n 1 {input} > {output} &&
        grep "Y:" {input} >> {output}
        """


rule arriba_plot:
    input:
        bam = rules.star_aln.output,
        bai = rules.samtools_idx.output,
        fusions = rules.arriba_y_fusions_filt.output,
        genome = custom_genome('results/custom-genome/combined.fasta'),
        gtf = custom_genome('results/custom-genome/combined.fixed.gtf')
    output:
        "results/arriba-viz/{sample}/fusions.pdf"
    singularity:
        "docker://quay.io/biocontainers/arriba:2.1.0--h3198e80_1"
    shell:
        """
	draw_fusions.R --fusions={input.fusions} --alignments={input.bam}/Aligned.sortedByCoord.out.bam --annotation={input.gtf} --output={output}
        """

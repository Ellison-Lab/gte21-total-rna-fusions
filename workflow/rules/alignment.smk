def get_fqs(s,end="r1"):
    tmp = SUBSAMPLE_TABLE[SUBSAMPLE_TABLE['sample_name'] == s]
    fqs = tmp.get("fastq_"+ end)
    return fqs

rule get_fastqs:
    output:
        r1 = "results/fastq/{sample}/{subsample}_r1.fastq.gz",
        r2 = "results/fastq/{sample}/{subsample}_r2.fastq.gz"
    params:
        r1 = lambda wc: get_fastqs(wc, "r1"),
        r2 = lambda wc: get_fastqs(wc, "r2"),
    shell:
        """
        wget -O {output.r1} {params.r1} &&
        wget -O {output.r2} {params.r2}
        """

rule concat_fqs:
    input:
        lambda wc: expand("results/fastq/{s}/{ss}_{end}.fastq.gz",s=wc.sample,end=wc.end, ss = SUBSAMPLE_TABLE[SUBSAMPLE_TABLE['sample_name'] == wc.sample]["sample_name"].tolist())
    output:
        temp("results/fastq-concat/{sample}_{end}.fq.gz")
    shell:
        "cat {input} > {output}"

rule trim_pe:
    input:
        #fq1 = lambda wc: get_from_pep(wc,"fastq_r1"),
        #fq2 = lambda wc: get_from_pep(wc,"fastq_r2"),
        fq1 = "results/fastq-concat/{sample}_r1.fq.gz",
        fq2 = "results/fastq-concat/{sample}_r2.fq.gz"
    output:
        r1 = "results/fastq/{sample}_r1.trimmed.fq.gz",
        r2 = "results/fastq/{sample}_r2.trimmed.fq.gz",
        html = "results/fastq/{sample}_fastp.html",
        json = "results/fastq/{sample}_fastp.json"
    threads:
        12
    resources:
        time=60,
        mem=20000,
        cpus=12
    conda:
        "../envs/fastp.yaml"
    singularity:
        "docker://quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
    shell:
        "fastp --in1 {input.fq1} --in2 {input.fq2} "
        "--out1 {output.r1} --out2 {output.r2} "
        "-j {output.json} -h {output.html} "
        "-w {threads} -L -R {wildcards.sample}_fastp"


rule star_idx:
    input:
        fa = custom_genome('results/custom-genome/combined.fasta'),
        gtf = custom_genome('results/custom-genome/combined.fixed.gtf')
    output:
        directory("results/idx/star")
    singularity:
        "docker://quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
    threads:
        24
    resources:
        time=60,
        mem=20000,
        cpus=24
    shell:
        """
        mkdir -p {output} &&
        STAR --runMode genomeGenerate --genomeFastaFiles {input.fa} \
            --genomeDir {output} --sjdbGTFfile {input.gtf} \
            --runThreadN {threads} --sjdbOverhang 149 --genomeSAsparseD 1 \
            --genomeSAindexNbases 12
        """

rule star_aln:
    input:
        idx=rules.star_idx.output,
        r1=rules.trim_pe.output.r1,
        r2=rules.trim_pe.output.r2
    output:
        directory("results/star/{sample}/")
    threads:
        12
    singularity:
        "docker://quay.io/biocontainers/star:2.7.9a--h9ee0642_0"
    resources:
        time=960,
        mem=64000,
        cpus=12
    shell:
        """
        mkdir -p {output}
        STAR \
        --genomeDir {input.idx} --runThreadN {threads} \
        --outFileNamePrefix {output}/ \
        --genomeLoad NoSharedMemory --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
        --readFilesIn {input.r1} {input.r2} \
        --outFilterMultimapNmax 50 \
        --peOverlapNbasesMin 10 \
        --alignTranscriptsPerReadNmax 50000 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 \
        --chimOutType WithinBAM HardClip \
        --chimJunctionOverhangMin 10 \
        --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3 \
        --chimMultimapNmax 50
        """

rule samtools_idx:
    input:
        rules.star_aln.output
    output:
        touch("results/{sample}.bam.indexed")
    #singularity:
        #"docker://quay.io/biocontainers/samtools:0.1.19--h270b39a_9"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools index -@ {threads} {input}/Aligned.sortedByCoord.out.bam
        """

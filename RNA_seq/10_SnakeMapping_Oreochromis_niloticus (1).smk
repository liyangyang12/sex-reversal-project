#!/usr/bin/env runsnakemake
configfile: "config2.yaml"
samples = config["samples"]
threads = config["threads"]
indir = "data2"
outdir = "results/mapping"

rule all:
    input:
        # Mapping
        outdir + "/star_index",
        expand(outdir + "/mapped.1st/{sample}", sample=samples),
        expand(outdir + "/mapped.2nd/{sample}", sample=samples),
        # Filter
        expand(outdir + "/filtered/{sample}.bam", sample=samples),
        expand(outdir + "/filtered/{sample}.bam.bai", sample=samples),
       

#trim_galore version 0.6.6
rule removing_adapter:
    input:
        fq1 = indir + "/{sample}_1.fastq.gz",
        fq2 = indir + "/{sample}_2.fastq.gz"
    output:
        fq1 = "results/trimmed/{sample}_1_val_1.fq.gz",
        fq2 = "results/trimmed/{sample}_2_val_2.fq.gz",
    shell:
        """
        trim_galore -q 20 --stringency 3 --length 20 -e 0.1 --paired {input.fq1} {input.fq2} --gzip -o results/trimmed
        """

# Mapping to genome
# https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
# STAR 2.7.6a
rule star_index:
    input:
        fas = "genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna",
        gtf = "genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gff"
    output:
        directory(outdir + "/star_index")
    log:
        outdir + "/star_index.log"
    threads:
        threads
    shell:
        """
        mkdir {output}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --sjdbGTFfile {input.gtf} \
            --sjdbGTFtagExonParentTranscript Parent \
            --genomeSAindexNbases 13 \
            --genomeFastaFiles {input.fas} &> {log}
        """

rule star_mapping_1st:
    input:
        fq1 = "results/trimmed/{sample}_1_val_1.fq.gz",
        fq2 = "results/trimmed/{sample}_2_val_2.fq.gz",
        idx = rules.star_index.output
    output:
        directory(outdir + "/mapped.1st/{sample}")
    log:
        outdir + "/mapped.1st/{sample}.log"
    threads:
        threads
    params:
        prefix = outdir + "/mapped.1st/{sample}",
    shell:
        """
        mkdir {output}
        STAR --runMode alignReads \
            --outSAMtype None \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {output}/ \
            --genomeDir {input.idx} \
            --readFilesIn {input.fq1} {input.fq2} &> {log}
        """

rule star_index2pass:
    input:
        fas = "genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna",
        gtf = "genome/GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.gff",
        tag = directory(expand(outdir + "/mapped.1st/{sample}",sample=samples))
    output:
        directory(outdir + "/star_index2pass")
    log:
        outdir + "/star_index2pass.log"
    threads:
        threads
    params:
        tabs = " ".join([outdir + "/mapped.1st/%s/SJ.out.tab" % sample for sample in samples])
    shell:
        """
        mkdir {output}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --sjdbGTFfile {input.gtf} \
            --sjdbGTFtagExonParentTranscript Parent \
            --genomeSAindexNbases 13 \
            --limitSjdbInsertNsj 1131226 \
            --sjdbFileChrStartEnd {params.tabs} \
            --genomeFastaFiles {input.fas} &> {log}
        """

rule star_mapping_2nd:
    input:
        fq1 = "results/trimmed/{sample}_1_val_1.fq.gz",
        fq2 = "results/trimmed/{sample}_2_val_2.fq.gz", 
        idx = outdir + "/star_index2pass",
    output:
        directory(outdir + "/mapped.2nd/{sample}")
    log:
        outdir + "/mapped.2nd/{sample}.log"
    threads:
        threads
    params:
        prefix = outdir + "/mapped.2nd/{sample}",
        tabs = " ".join([outdir + "/mapped.1st/%s/SJ.out.tab" % sample for sample in samples])
    shell:
        """
        mkdir {output}
        STAR --runMode alignReads \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS nM NM MD jM jI MC ch XS \
            --limitBAMsortRAM 10000000000 \
            --limitSjdbInsertNsj 2000000 \
            --readFilesCommand zcat \
            --runThreadN {threads} \
            --outFileNamePrefix {output}/ \
            --genomeDir {input.idx} \
            --sjdbFileChrStartEnd {params.tabs} \
            --readFilesIn {input.fq1} {input.fq2} &> {log}
        """


# Filtered
# bamtools 2.5.1

rule filter_alignment:
    input:
        rules.star_mapping_2nd.output,
    output:
        bam = outdir + "/filtered/{sample}.bam"
    shell:
        """
        bamtools filter -in {input}/Aligned.sortedByCoord.out.bam -out {output} \
            -tag NH:1 -isProperPair true -isPrimaryAlignment true
        """

# Common rules

rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    shell:
        """
        bamtools index -in {input}
        """


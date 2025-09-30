SAMPLE=["PVO_rep1","VO_rep1","VO_rep2","VO_rep3","MO_rep1","MO_rep2","MO_rep3","ISE_rep1","ISE_rep2","ISE_rep3","ISM_rep1","ISM_rep2","ISM_rep3","ISL_rep1","ISL_rep2","ISL_rep3","ES_rep1","ES_rep2","ES_rep3","MS_rep1","MS_rep2","MS_rep3","LS_rep1","LS_rep2","LS_rep3"]
reads=["1","2"]
SAMPLE2=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]
rule all:
    input:
        "genome_methylpy/malbus_r.4.bt2",
        expand("results/methylpy/merge/allc_{name}.tsv.gz",name = SAMPLE2),
        "results/methylpy/DMR/dmr_rms_results.tsv.gz",
        "results/methylpy/DMR_reidentify/dmr_rms_results_recollapsed.tsv",
        "results/methylpy/methylation-level/all_dmr_methylationlevel.tsv",
        

rule methylpy_build_reference:
    input:
        genome="genome_methylpy/malbus.fa"
    output:
        "genome_methylpy/malbus_r.4.bt2"
    params:
        prefix="genome_methylpy/malbus"
    shell:
        """
        methylpy build-reference --input-files {input.genome} --output-prefix {params.prefix} --num-procs 40
        """

rule methylpy_processing:
    input:
        genome="genome_methylpy/malbus.fa",
        r1 = "data/{sample}_1.fq.gz",
        r2 = "data/{sample}_2.fq.gz"
    output:
        directory("results/methylpy/processing/{sample}")
    log:
        "results/log_methylpy/{sample}.log"
    shell:
        """
        methylpy paired-end-pipeline --read1-files {input.r1} --read2-files {input.r2} --sample {wildcards.sample} --forward-ref genome_methylpy/malbus_f --reverse-ref genome_methylpy/malbus_r --ref-fasta {input.genome} --num-procs 20 --remove-clonal True --path-to-picard="/home/lyy/miniconda3/pkgs/picard-2.25.2-hdfd78af_0/share/picard-2.25.2-0/" --path-to-output {output} &> {log}
        """

rule merge:
    input:
        rep1 = "results/methylpy/processing/{name}_rep1/allc_{name}_rep1.tsv.gz",
        rep2 = "results/methylpy/processing/{name}_rep2/allc_{name}_rep2.tsv.gz",
        rep3 = "results/methylpy/processing/{name}_rep3/allc_{name}_rep3.tsv.gz",
    output:
        "results/methylpy/merge/allc_{name}.tsv.gz",
    log:
        "results/methylpy/merge/{name}.log"
    shell:
        """
        methylpy merge-allc \
	--allc-files {input.rep1} {input.rep2} {input.rep3} \
	--output-file {output} \
	--num-procs 20 \
	--compress-output True &> {log}
        """

rule DMR:
    input:
        expand("results/methylpy/merge/allc_{name}.tsv.gz",name=SAMPLE2),
    output:
        "results/methylpy/DMR/dmr_rms_results.tsv.gz"
    log:
        "results/methylpy/dmr.log"
    shell:
        """
        methylpy DMRfind --allc-files {input} --num-procs 60 --output-prefix results/methylpy/DMR/dmr &> {log}
        """


rule bw:
    input:
        allc="results/methylpy/merge/allc_{name}.tsv.gz",
        genome="genome_methylpy/malbus.fa",
    output:
        "results/methylpy/bw/{name}_meth.bw",
    log:
        "results/methylpy/bw/{name}.log",
    shell:
        """
        methylpy allc-to-bigwig \
	--allc-file {input.allc} \
	--output-file {output} \
	--ref-fasta {input.genome} \
	--mc-type CGN \
	--bin-size 1 &> {log}
        """


rule bismark_reference:
    input:
        genome="bismarkgenome/"
    output:
        "bismarkgenome/Bisulfite_Genome/CT_conversion/BS_CT.2.bt2"
    log:
        "bismarkgenome/bismark_reference.log"
    shell:
        """
        bismark_genome_preparation {input.genome} &> {log}
        """
        
rule bismark_run:
    input:
        r1 = "data/{sample}_1.fq.gz",
        r2 = "data/{sample}_2.fq.gz"
    output:
        directory("results/bismark/processing/{sample}")
    shell:
        """
        bismark --parallel 10 --output_dir {output} --gzip /beegfs/zhoulab/lyy/dna_methylation/bismarkgenome/ -1 {input.r1} -2 {input.r2} &> "results/log/{wildcards.sample}.log"
        """


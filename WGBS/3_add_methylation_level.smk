SAMPLE=["PVO_rep1","VO_rep1","VO_rep2","VO_rep3","MO_rep1","MO_rep2","MO_rep3","ISE_rep1","ISE_rep2","ISE_rep3","ISM_rep1","ISM_rep2","ISM_rep3","ISL_rep1","ISL_rep2","ISL_rep3","ES_rep1","ES_rep2","ES_rep3","MS_rep1","MS_rep2","MS_rep3","LS_rep1","LS_rep2","LS_rep3"]
reads=["1","2"]
SAMPLE2=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]
rule all:
    input:
        "results/methylpy/methylation-level_every_stage/all_genome_methylationlevel.tsv",
        

rule add_methylation_level_dmr:
    input:
        allc = expand("../results/methylpy/processing/{sample}/allc_{sample}.tsv.gz",sample = SAMPLE),
        tsv = "../genome_methylpy/genome_region.tsv",   
    output:
        "results/methylpy/methylation-level_every_stage/all_genome_methylationlevel.tsv",
    shell:
        """
        methylpy add-methylation-level \
        --input-tsv-file {input.tsv} \
        --output-file {output} \
        --allc-files {input.allc} \
        --mc-type CGN \
        --num-procs 60
        """

SAMPLE2=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]

rule all:
    input:
        "results/diffbind_DAR_bed/result_filterbound_sort.bed",
        "results/diffbind_DAR_bed/result_filterbound_sort_methylationlevel.tsv"

rule sort:
    input:
        "results/diffbind_DAR_bed/result_filterbound.bed"
    output:
        "results/diffbind_DAR_bed/result_filterbound_sort.bed"
    shell:
        """
        sort -k1,1 -k2,2n {input} > {output}
        """

rule add_methylation_level_promotor:
    input:
        allc = expand("../results/methylpy/merge/allc_{sample}.tsv.gz",sample = SAMPLE2),
        bed = "results/diffbind_DAR_bed/result_filterbound_sort.bed"
    output:
        "results/diffbind_DAR_bed/result_filterbound_sort_methylationlevel.tsv"
    shell:
        """
        methylpy add-methylation-level \
        --input-tsv-file {input.bed} \
        --output-file {output} \
        --allc-files {input.allc} \
        --mc-type CGN \
        --num-procs 60 \
        --input-no-header True
        """


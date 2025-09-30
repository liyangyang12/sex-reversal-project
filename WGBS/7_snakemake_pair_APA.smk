SAMPLE=["PVO_rep1","VO_rep1","VO_rep2","VO_rep3","MO_rep1","MO_rep2","MO_rep3","ISE_rep1","ISE_rep2","ISE_rep3","ISM_rep1","ISM_rep2","ISM_rep3","ISL_rep1","ISL_rep2","ISL_rep3","ES_rep1","ES_rep2","ES_rep3","MS_rep1","MS_rep2","MS_rep3","LS_rep1","LS_rep2","LS_rep3"]
reads=["1","2"]
SAMPLE2=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]

diff = []

for i in range(0,len(SAMPLE2)-1):
    for j in range(i+1,len(SAMPLE2)):
        diff.append(SAMPLE2[i] + '_vs_' + SAMPLE2[j])
print(diff)
diff2 = []

for i in range(0,len(SAMPLE2)-1):
    diff2.append(SAMPLE2[i] + '_vs_' + SAMPLE2[i+1])


rule all:
    input:
        #expand("results/methylpy/DMRpair/{diffname}/dmr_rms_results.tsv.gz",diffname=diff),
        expand("results/methylpy/DMRpair/DMRfilter/{diffname}_dmr_filter.txt",diffname=diff), 
        expand("results/methylpy/DMR/DMRfilter/{diffname}_dmr_filter.txt",diffname=diff),
        "results/methylpy/DMR/DMRfilter_allstage/allstage_dmr_filter.txt",
        expand("results/methylpy/DMR/DMRfilter_length/{diffname}_dmr_filter.txt",diffname=diff)



rule DMR_pair:
    input:
        allc1 = "results/methylpy/merge/allc_{sample1}.tsv.gz",
        allc2 = "results/methylpy/merge/allc_{sample2}.tsv.gz"
    output:
        "results/methylpy/DMRpair/{sample1}_vs_{sample2}/dmr_rms_results.tsv.gz"
    log:
        "results/methylpy/DMRpair/{sample1}_vs_{sample2}/dmr.log"
    shell:
        """
        methylpy DMRfind --allc-files {input.allc1} {input.allc2} --num-procs 8 --output-prefix results/methylpy/DMRpair/{wildcards.sample1}_vs_{wildcards.sample2}/dmr &> {log}
        """



rule filterDMR:
    input:
        "../results/methylpy/DMRpair/{sample1}_vs_{sample2}/dmr_rms_results_collapsed.tsv"
    output:
        "results/methylpy/DMRpair/DMRfilter/{sample1}_vs_{sample2}_dmr_filter.txt"
    shell:
        """
        python scripts/filterDMR_v2.py --dmr {input} --result {output}
        """

rule filterDMR_allstage:
    input:
        "../results/methylpy/DMR/dmr_rms_results_collapsed.tsv"
    output:
        "results/methylpy/DMR/DMRfilter/{sample1}_vs_{sample2}_dmr_filter.txt"
    shell:
        """
        python scripts/filterDMR_v3.py --dmr {input} --result {output} --sample1 {wildcards.sample1} --sample2 {wildcards.sample2}
        """

rule filterDMR_allstage_length:
    input:
        "../results/methylpy/DMR/dmr_rms_results_collapsed.tsv"
    output:
        "results/methylpy/DMR/DMRfilter_length/{sample1}_vs_{sample2}_dmr_filter.txt"
    shell:
        """
        python scripts/filterDMR_v5.py --dmr {input} --result {output} --sample1 {wildcards.sample1} --sample2 {wildcards.sample2}
        """

rule filterDMR_allstage_dapars2:
    input:
        "../results/methylpy/DMR/dmr_rms_results_collapsed.tsv"
    output:
        "results/methylpy/DMR/DMRfilter_allstage/allstage_dmr_filter.txt"
    shell:
        """
        python scripts/filterDMR_v4.py --dmr {input} --result {output} 
        """

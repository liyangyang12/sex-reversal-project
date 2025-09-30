

SAMPLE2=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]
diff2 = []
for i in range(0,len(SAMPLE2)-1):
    for j in range(i+1,len(SAMPLE2)):
        diff2.append(SAMPLE2[i] + '_vs_' + SAMPLE2[j])
print(len(diff2))

rule all:
    input:
        "results/peak_merge/result_filterbound.bed",


rule mergebed:
    input:
        expand("results/DARfilter/{sample}_dmr_filter.txt",sample=diff2)
    output:
        merge1 = "results/peak_merge/merge.bed",
        merge2 = "results/peak_merge/result.bed",
    shell:
        """
        cat {input} > {output.merge1}
        sort -k1,1 -k2,2n {output.merge1} |cut -f 1,2,3|sort| uniq > {output.merge2}
        """


rule filterbed:
    input:
        txt = "results/peak_merge/result.bed",
        chrsixe = "/data/liyangyang/ATAC-seq/genome/chrsize.txt"
    output:
        "results/peak_merge/result_filterbound.bed",
    shell:
        """
        python ../scripts/filter_bound.py --txt {input.txt} --chrsize {input.chrsixe} --result {output}
        """
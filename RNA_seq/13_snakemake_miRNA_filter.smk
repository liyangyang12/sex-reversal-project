rule all:
    input:
        "results/mfuzz_apa/cluster6_filter_dpas.fa",
        "results/miRNA/malbus_miR_Family_Info.fa",
        "results/miranda/filter_predict_malbus.txt",
        "results/miranda/stats.txt",
        "results/miranda/top5_TPM.txt"


#cluster6_pdui_filter.txt删除了title
rule filter_cluster:
    input:
        "results/mfuzz_apa/cluster6_pdui_filter.txt"
    output:
        "results/mfuzz_apa/cluster6_filter.txt"
    shell:
        """
        cat {input} |cut -f 1 | cut -d "|" -f 1 > {output}
        """

rule dpas_bed:
    input:
        cluster = "results/mfuzz_apa/cluster6_filter.txt",
        pdui = "../results/dapar2_merge/all_PDUI_NA.txt",
    output:
        result = "results/mfuzz_apa/cluster6_filter_dpas.bed"
    shell:
        """
        python ../scripts/dpAS_bed.py --cluster {input.cluster} --pdui {input.pdui} --result {output.result}
        """

#有些基因是定义在forward链上的，意思就是基因对应的转录本序列刚好和forward链上5‘到3’的碱基序列一致，而另一些基因定义在reverse链上，就是说，这些基因的转录本序列（以及对应的氨基酸序列）和reverse链上5‘到3’方向的序列一致

rule get_fasta:
    input:
        bed = "results/mfuzz_apa/cluster6_filter_dpas.bed",
        fa = "../../1_expression_dynamics/genome/GCF_001952655.1_M_albus_1.0_genomic.fna"
    output:
        result = "results/mfuzz_apa/cluster6_filter_dpas.fa"
    shell:
        """
        bedtools getfasta -fi {input.fa} -bed {input.bed} -name -s  > {output.result}
        """


#cat predict.txt | grep '^>>'|sed  's/^M//g'| sed 's/>>//g' >miRanda_final_res.txt
#https://zhuanlan.zhihu.com/p/406919586
rule get_malbus_mirfamily:
    input:
        "../miRNA/Dapeng_Li_2014/malmiRNA.txt"
    output:
        result = "results/miRNA/malbus_miR_Family_Info.fa"
    shell:
        """
        python ../scripts/get_malbus_mir.py --mir1 {input} --mir2 {output.result}
        """

#conda activate py34
#miranda v3.3a 
rule mirandamalbus:
    input:
        mirna = "results/miRNA/malbus_miR_Family_Info.fa",
        utr = "results/mfuzz_apa/cluster6_filter_dpas.fa"
    output:
        "results/miranda/filter_predict_malbus.txt"
    shell:
        """
        miranda {input.mirna} {input.utr} -out {output} -en -1 
        """

#>和>>区别见filter_predict_malbus.txt mal-miR-101a  一个microRNA可能在一段序列上有多个位点，>>将这些位点结合起来，最后统计的是每个mircoRNA调控多少APA转录本
rule filter_miranda_result:
    input:
        "results/miranda/filter_predict_malbus.txt"
    output:
        result = "results/miranda/miRanda_final_filter.txt",
        stat = "results/miranda/stats.txt"
    shell:
        """
        cat {input} | grep '>>'|sed 's/>>//g' > {output.result}
        cat {output.result} |cut -f 1 | sort |uniq -c | sort -rn > {output.stat}
        """

#筛选靶基因top5的microRNA,获取表达量TPM
rule filter_top5:
    input:
        "../miRNA/Dapeng_Li_2014/malmiRNA.txt"
    output:
        "results/miranda/top5_TPM.txt"
    shell:
        """
        cat {input} | grep -E "mal-mir-737|mal-miR-203a|mal-miR-144|mal-miR-203b|mal-miR-2184" | cut -f 1,7,8,9 > {output}
        """

#2023/6/10  没跑后面的
#用DNA序列算和用RNA序列算结果一样，只有几个热力学稳定性分分数不同，没影响最后统计结果
rule get_RNA_sequence:
    input:
        "results/mfuzz_apa/cluster6_filter_dpas.fa"
    output:
        "results/mfuzz_apa/cluster6_filter_dpas.rna.fa"
    shell:
        """
        python scripts/get_rna_v2.py --fa {input} --result {output} 
        """

rule mirandamalbusrna:
    input:
        mirna = "results/miRNA/malbus_miR_Family_Info.fa",
        utr = "results/mfuzz_apa/cluster6_filter_dpas.rna.fa"
    output:
        "results/miranda/filter_predict_malbus_rna.txt"
    shell:
        """
        miranda {input.mirna} {input.utr} -out {output} -en -1 
        """
#############################################################


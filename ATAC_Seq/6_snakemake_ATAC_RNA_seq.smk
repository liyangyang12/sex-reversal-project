SAMPLE3 = ["PVO_rep1","VO_rep1","VO_rep2","MO_rep1","MO_rep3","ISE_rep1","ISE_rep2","ISM_rep1","ISM_rep3","ISL_rep1","ISL_rep2","ES_rep3","ES_rep4","MS_rep4","MS_rep6","LS_rep1","LS_rep3"]

rule all:
    input:
        "results/closet/closet.bed",
        "results/closet/closet_count_cpm.txt"




rule get_promotor:
    input:
        gff = "../genome/GCF_001952655.1_M_albus_1.0_genomic.gff"
    output:
        result = "results/promoter/malbus_promoter.bed",
        result2 = "results/promoter/malbus_promoter_sort.bed",
        result3 = "results/promoter/malbus_promoter_sort_filter.bed"
    shell:
        """
        awk 'BEGIN{FS="\t|=|;"}($3=="gene"){OFS="\t"; if ($7~/+/){print $1,$4-2000,$4+500,$14}; if ($7~/-/){print $1,$5-500,$5+2000,$14}}' {input.gff}| sed 's/[";]//g;' > {output.result}
        sort -k1,1 -k2,2n {output.result} > {output.result2}
        awk '{OFS="\t"; if ($2>0){print $1,$2,$3,$4}}' {output.result2} > {output.result3}
        """

rule closet:
    input:
        promotor = "results/promoter/malbus_promoter_sort_filter.bed",
        DAR = "results/peak_merge/result_filterbound_cpm.txt"
    output:
        sortDAR = "results/peak_merge/result_filterbound_sort_cpm.txt",
        result = "results/closet/closet.txt",
        filte = "results/closet/closet.bed"
    shell:
        """
        sort -k1,1 -k2,2n {input.DAR} > {output.sortDAR}
        bedtools closest -d -a {output.sortDAR} -b {input.promotor} > {output.result}
        cat {output.result}  | awk  '$26<=50000{{print $0}}'|cut -f 1,2,3,4,25 > {output.filte}
        """


rule multicov:
    input:
        bam = expand("../results/filter/{sample}.filter.bam",sample=SAMPLE3),
        bed = "results/closet/closet.bed",
    output:
        txt = "results/closet/closet_count.txt"
    shell:
        """
        bedtools multicov -bams {input.bam} -bed {input.bed} > {output}
        """

rule bam_stat:
    input:
        "../results/filter/{sample}.filter.bam"
    output:
        "results/stats/{sample}.stats"
    shell:
        """
        bamtools stats -in {input} > {output}
        """

rule get_library_size:
    input:
        txt = "results/stats/{sample}.stats"
    output:
        txt = "results/stats/{sample}.library_size.txt"
    shell:
        """
        echo {wildcards.sample} "`cat {input.txt} | grep 'Total reads:' | awk '{{print $3}}'`" > {output.txt}
        """

rule merge_library:
    input:
        expand("results/stats/{sample}.library_size.txt",sample=SAMPLE3)
    output:
        "results/stats/library_size.txt"
    shell:
        """
        cat {input} > {output}
        """

rule normalize:
    input:
        count = "results/closet/closet_count.txt",
        library = "results/stats/library_size.txt"
    output:
        "results/closet/closet_count_cpm.txt"
    shell:
        """
        python ../scripts/calculate_cpm2.py --count {input.count} --library {input.library} --result {output}
        """

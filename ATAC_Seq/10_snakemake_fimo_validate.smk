SAMPLE=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]

rule all:
    input:
        expand("results/fimo_prepare/{sample}_fimo_prepare.fa",sample=SAMPLE),
        directory(expand("results/fimo/{sample}",sample=SAMPLE))



rule fimoprepare:
    input:
        promotor = "results/promoter/malbus_promoter_sort_filter.bed",
        file2 = "results/ananse_binding/{sample}/regions_combined.bed",
        fa = "../genome/GCF_001952655.1_M_albus_1.0_genomic.fna",
    output:
        sort = "results/fimo_prepare/{sample}_sort.bed",
        file3 = "results/fimo_prepare/{sample}_closet.bed",
        prepare = "results/fimo_prepare/{sample}_prepare.bed",
        file4 = "results/fimo_prepare/{sample}_fimo_prepare.fa"
    shell:
        """
        sort -k1,1 -k2,2n {input.file2} > {output.sort}
        bedtools closest -d -a {output.sort} -b {input.promotor} > {output.file3}
        cat {output.file3} | awk  '$8<=50000{{print $0}}'|cut -f 1,2,3,7 > {output.prepare}
        bedtools getfasta -name -fi {input.fa} -bed {output.prepare} > {output.file4}
        """

rule fimo:
    input:
        fa = "results/fimo_prepare/{sample}_fimo_prepare.fa",
        meme = "JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme"
    output:
        result = directory("results/fimo/{sample}")
    shell:
       """
       fimo --o results/fimo/{wildcards.sample} --max-stored-scores 10000000  --thresh 1e-5 --verbosity 4 {input.meme} {input.fa} 
       """
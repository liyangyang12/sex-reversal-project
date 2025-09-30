SAMPLE2=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]
loc = ["promoter","genebody"]

rule all:
    input:
        #"results/promoter/malbus_promotor_sort_filter.bed",
        #"results/gene_body/malbus_gene_body.bed",
        #"results/methylpy/promoter_methylation-level/promoter_methylationlevel.tsv",
        #"results/methylpy/genebody_methylation-level/genebody_methylationlevel.tsv",
        #expand("results/cor_methlevel_exp/{loc}_positive.txt",loc=loc),
        "results/methylpy/promoter_methylation-level/promoter_DMR_methylationlevel.tsv",
        "results/methylpy/genebody_methylation-level/genebody_DMR_methylationlevel.tsv",
        "results/methylpy/promoter_methylation-level_filter/promoter_DMR_methylationlevel.tsv",
        expand("results/cor_DMR_methlevel_exp/{loc}_positive.txt",loc=loc),
        expand("results/cor_DMR_methlevel_exp_length/{loc}_positive.txt",loc=loc),




rule get_promotor:
    input:
        gff = "../../perl/1_expression_dynamics/genome/GCF_001952655.1_M_albus_1.0_genomic.gff"
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

rule get_genebody:
    input:
        gff = "../../perl/1_expression_dynamics/genome/GCF_001952655.1_M_albus_1.0_genomic.gff"
    output:
        result = "results/gene_body/malbus_gene_body.bed",
        result2 = "results/promoter/malbus_gene_body_sort.bed",
        result3 = "results/promoter/malbus_gene_body_sort_filter.bed"
    shell:
        """
        awk 'BEGIN{FS="\t|=|;"}($3=="gene"){OFS="\t"; if ($7~/+/){print $1,$4+500,$5,$14}; if ($7~/-/){print $1,$4,$5-500,$14}}' {input.gff}| sed 's/[";]//g;' > {output.result}
        sort -k1,1 -k2,2n {output.result} > {output.result2}
        awk '{OFS="\t"; if ($2>0){print $1,$2,$3,$4}}' {output.result2} > {output.result3}
        """

rule add_methylation_level_promotor:
    input:
        allc = expand("../results/methylpy/merge/allc_{sample}.tsv.gz",sample = SAMPLE2),
        promoter = "results/promoter/malbus_promoter_sort_filter.bed"
    output:
        "results/methylpy/promoter_methylation-level/promoter_methylationlevel.tsv"
    shell:
        """
        methylpy add-methylation-level \
        --input-tsv-file {input.promoter} \
        --output-file {output} \
        --allc-files {input.allc} \
        --mc-type CGN \
        --num-procs 60
        """

rule add_methylation_level_genebody:
    input:
        allc = expand("../results/methylpy/merge/allc_{sample}.tsv.gz",sample = SAMPLE2),
        gene_body = "results/gene_body/malbus_gene_body.bed"
    output:
       "results/methylpy/genebody_methylation-level/genebody_methylationlevel.tsv"
    shell:
        """
        methylpy add-methylation-level \
        --input-tsv-file {input.gene_body} \
        --output-file {output} \
        --allc-files {input.allc} \
        --mc-type CGN \
        --num-procs 60
        """

rule gb_pr_exp_cor:
    input:
        counts = "../../perl/2_alternative_splicing_Whippet/exp_psi_file/counts_local.txt",
        methlevel = "results/methylpy/{loc}_methylation-level/{loc}_methylationlevel.tsv"
    output:
        positive = "results/cor_methlevel_exp/{loc}_positive.txt",
        negative = "results/cor_methlevel_exp/{loc}_negative.txt",
        allresult = "results/cor_methlevel_exp/{loc}_all.txt",
    shell:
        """
        Rscript scripts/methylation_level_exp_cor.R {input.counts} {input.methlevel} {output.positive} {output.negative} {output.allresult}
        """


rule get_promoter_DMR_intersect:
    input:
        promoter = "results/promoter/malbus_promoter_sort_filter.bed",
        DMR = "results/methylpy/promoter_methylation-level/dmr_rms_results_collapsed.tsv"
    output:
        "results/methylpy/promoter_methylation-level/promoter_DMR_methylationlevel.tsv"
    shell:
        """
        bedtools intersect -a {input.promoter} -b {input.DMR} -wa -wb -F 0.10 > {output}
        """

rule get_genebody_DMR_intersect:
    input:
        genebody = "results/gene_body/malbus_gene_body_sort_filter.bed",
        DMR = "results/methylpy/promoter_methylation-level/dmr_rms_results_collapsed.tsv"
    output:
        "results/methylpy/genebody_methylation-level/genebody_DMR_methylationlevel.tsv"
    shell:
        """
        bedtools intersect -a {input.genebody} -b {input.DMR} -wa -wb  -F 0.10 > {output}
        """

rule filter_DMR_intersect:
    input:
        promoter = "results/methylpy/promoter_methylation-level/promoter_DMR_methylationlevel.tsv",
        genebody = "results/methylpy/genebody_methylation-level/genebody_DMR_methylationlevel.tsv"
    output:
        filter_promoter = "results/methylpy/promoter_methylation-level_filter/promoter_DMR_methylationlevel.tsv",
        filter_genebody = "results/methylpy/genebody_methylation-level_filter/genebody_DMR_methylationlevel.tsv"
    shell:
        """
        Rscript scripts/DMR_pro_gb_classify.R {input.promoter} {input.genebody} {output.filter_promoter} {output.filter_genebody}
        """

rule gb_pr_DMR_exp_cor:
    input:
        counts = "../../perl/2_alternative_splicing_Whippet/exp_psi_file/counts_local.txt",
        methlevel = "results/methylpy/{loc}_methylation-level_filter/{loc}_DMR_methylationlevel.tsv"
    output:
        positive = "results/cor_DMR_methlevel_exp/{loc}_positive.txt",
        negative = "results/cor_DMR_methlevel_exp/{loc}_negative.txt",
        allresult = "results/cor_DMR_methlevel_exp/{loc}_all.txt",
    shell:
        """
        Rscript scripts/DMR_methylation_level_exp_cor.R {input.counts} {input.methlevel} {output.positive} {output.negative} {output.allresult}
        """

rule gb_pr_DMR_exp_cor_length:
    input:
        counts = "../../perl/2_alternative_splicing_Whippet/exp_psi_file/counts_local.txt",
        methlevel = "results/methylpy/{loc}_methylation-level_filter/{loc}_DMR_methylationlevel.tsv"
    output:
        positive = "results/cor_DMR_methlevel_exp_length/{loc}_positive.txt",
        negative = "results/cor_DMR_methlevel_exp_length/{loc}_negative.txt",
        allresult = "results/cor_DMR_methlevel_exp_length/{loc}_all.txt",
    shell:
        """
        Rscript scripts/DMR_methylation_level_exp_cor_length.R {input.counts} {input.methlevel} {output.positive} {output.negative} {output.allresult}
        """


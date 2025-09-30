SAMPLE1 = ["O_I","O_IV","O_V","OT_I","OT_II","OT_III","T_I","T_III","T_IV"]
diff1 = []
for i in range(0,len(SAMPLE1)-1):
    for j in range(i+1,len(SAMPLE1)):
        diff1.append(SAMPLE1[i] + '_vs_' + SAMPLE1[j])
#print(diff1)

SAMPLE2=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]
diff2 = []
for i in range(0,len(SAMPLE2)-1):
    for j in range(i+1,len(SAMPLE2)):
        diff2.append(SAMPLE2[i] + '_vs_' + SAMPLE2[j])

diff3=[]
for i in range(0,len(diff1)):
    diff3.append(diff1[i] + '_' + diff2[i])



rule all:
    input:
        expand("results/dapars1_APA_dpas_bed/{diffname}.bed",diffname = diff1),
        expand("results/intersect_dapars1_APA_dpas_DMR/{diffname2}.txt",diffname2=diff3),
        #"results/intersect_dapars1_APA_dpas_DMR/homer/test.flag",
        "results/intersect_dapars1_APA_dpas_DMR/homer_DMR/test.flag",
        "results/intersect_dapars1_APA_dpas_DMR/fimo/results.txt",
        "results/dapars2_APA_dpas_cluster6/cluster6.bed",
        "results/intersect_dapars2_APA_dpas_DMR/intersect.txt",
        "results/intersect_dapars2_APA_dpas_DMR/homer_DMR/test.flag",
        "results/AME_veterbrate_nobg/test.flag",
        "results/intersect_dapars1_APA_dpas_DMR/background.fa",
        "results/AME_veterbrate_bg/test.flag",
        expand("results/intersect_dapars1_APA_dpas_DMR_length/{diffname2}.txt",diffname2=diff3),
        #expand("results/intersect_dapars1_APA_dpas_DMR/filter/{diffname2}_filter.txt",diffname2=diff3),
        "results/intersect_dapars1_APA_dpas_DMR/filter/homer_DMR/test.flag",
        "results/intersect_dapars1_APA_dpas_DMR/filter/fimo/results.txt",
        "results/computeMatrix/results.pdf",
        expand("results/computeMatrix3/{diffname2}.pdf",diffname2 = diff3)



rule get_dapars1_APA_dpas_bed:
    input:
        "results/dapars1/output/{sample1}_vs_{sample2}/DaPars_data_All_Prediction_Results.txt"
    output:
        bed = "results/dapars1_APA_dpas_bed/{sample1}_vs_{sample2}.bed",
    shell:
        """
        python scripts/dpAS_bed_dapars1.py --dapars1 {input} --result {output.bed}
        """


rule intersect_dapars1_APA_dpas_DMR:
    input:
        APA_dpas_bed ="results/dapars1_APA_dpas_bed/{sample1}_vs_{sample2}.bed",
        dmr = "../../dna_methylation/core_scripts/results/methylpy/DMR/DMRfilter/{sample3}_vs_{sample4}_dmr_filter.txt"
    output:
        "results/intersect_dapars1_APA_dpas_DMR/{sample1}_vs_{sample2}_{sample3}_vs_{sample4}.txt"
    shell:
        """
        bedtools intersect -a {input.APA_dpas_bed} -b {input.dmr} -wa -wb -F 0.70 > {output}
        """

rule filter_intersect_dapars1_APA_dpas_DMR:
    input:
        "results/intersect_dapars1_APA_dpas_DMR/{sample1}_vs_{sample2}_{sample3}_vs_{sample4}.txt"
    output:
        "results/intersect_dapars1_APA_dpas_DMR/filter/{sample1}_vs_{sample2}_{sample3}_vs_{sample4}_filter.txt"
    shell:
        """
        Rscript scripts/filter.R {input} {output}
        """

#cd /DATA/lyy/perl/3_alternative_polyadenylation/results/intersect_dapars1_APA_dpas_DMR/filter
#cat *.txt | awk -v OFS="\t" '{print $9,$10,$11,$4}' |sort |uniq > DMR_name.bed  至此还剩25个基因
rule homer4:
    input:
        fa = "../1_expression_dynamics/genome/GCF_001952655.1_M_albus_1.0_genomic.fna",
        bed = "results/intersect_dapars1_APA_dpas_DMR/filter/DMR_name.bed",
    output:
        tag = "results/intersect_dapars1_APA_dpas_DMR/filter/homer_DMR/test.flag"
    shell:
        """
        touch {output.tag}
        findMotifsGenome.pl {input.bed} {input.fa} results/intersect_dapars1_APA_dpas_DMR/filter/homer_DMR -size given  -p 20 
        """



rule get_fasta2:
    input:
        fa = "../1_expression_dynamics/genome/GCF_001952655.1_M_albus_1.0_genomic.fna",
        bed = "results/intersect_dapars1_APA_dpas_DMR/filter/DMR_name.bed",
    output:
        fa2 = "results/intersect_dapars1_APA_dpas_DMR/filter/DMR_name.fa"
    shell:
        """
        bedtools getfasta -fi {input.fa} -bed {input.bed} -name > {output}
        """

rule fimo2:
    input:
        fa = "results/intersect_dapars1_APA_dpas_DMR/filter/DMR_name.fa",
        meme = "../7_DMAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms.meme"
    output:
        result = "results/intersect_dapars1_APA_dpas_DMR/filter/fimo/results.txt"
    shell:
       """
       fimo --o results/intersect_dapars1_APA_dpas_DMR/filter/fimo --max-stored-scores 10000000 --text --thresh 1e-5 --verbosity 4 {input.meme} {input.fa} > {output.result}
       """

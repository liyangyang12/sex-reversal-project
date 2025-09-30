SAMPLE=["PVO_rep1","VO_rep1","VO_rep2","VO_rep3","MO_rep1","MO_rep2","MO_rep3","ISE_rep1","ISE_rep2","ISE_rep3","ISM_rep1","ISM_rep2","ISM_rep3","ISL_rep1","ISL_rep2","ISL_rep3","ES_rep1","ES_rep2","ES_rep3","MS_rep1","MS_rep2","MS_rep3","LS_rep1","LS_rep2","LS_rep3"]
reads=["1","2"]
SAMPLE2=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]

SAMPLE_RNA_seq = ["O_I","O_IV","O_V","OT_I","OT_II","OT_III","T_I","T_III","T_IV"]

configfile: "config.yaml"
samplesRNASEQ = config["samples"]

mergename = []
for i in range(0,len(SAMPLE2)):
    mergename.append(SAMPLE2[i] + '_vs_' + SAMPLE_RNA_seq[i])
print(mergename)


rule all:
    input:
        expand("results/stringtie_jarid2/{name1}.txt",name1 = samplesRNASEQ),
        "results/stringtie_jarid2/merge.gtf",
        expand("results/stringtie2/{name1}.txt",name1 = samplesRNASEQ),
        expand("results/stringtie3/{name1}.txt",name1 = samplesRNASEQ),
        expand("results/stringtie3/modified/{name1}.gtf",name1 = samplesRNASEQ),
        expand("results/stringtie2/modified/{name1}.gtf",name1 = samplesRNASEQ),

#https://blog.csdn.net/weixin_43569478/article/details/108079257
rule stringtie_RNA_seq:
    input:
        bam="../../1_expression_dynamics/core_scripts/results/mapping/filtered/{name1}.bam",
        gtf = "../../1_expression_dynamics/genome/GCF_001952655.1_M_albus_1.0_genomic_lgr4_modified.gtf"
    output:
        quan = "results/stringtie_jarid2/{name1}.txt",
        ballgown = directory("results/stringtie_jarid2/{name1}_ballgown"),
        gtf = "results/stringtie_jarid2/{name1}.gtf",
    shell:
        """
        stringtie {input.bam} -p 10 -G {input.gtf} -A {output.quan} -o {output.gtf} -b {output.ballgown}
        """

rule stringtie_merge:
    input:
        gtf = expand("results/stringtie_jarid2/{name1}.gtf",name1 = samplesRNASEQ),
        refgtf = "../../1_expression_dynamics/genome/GCF_001952655.1_M_albus_1.0_genomic_lgr4_modified.gtf"
    output:
        "results/stringtie_jarid2/merge.gtf"
    shell:
        """
        stringtie --merge  -o {output} -p 20 -G {input.refgtf} {input.gtf}
        """

rule stringtie_RNA_seq2:
    input:
        bam="../../1_expression_dynamics/core_scripts/results/mapping/filtered/{name1}.bam",
        gtf = "results/stringtie_jarid2/merge.gtf"
    output:
        quan = "results/stringtie2/{name1}.txt",
        ballgown = directory("results/stringtie2/{name1}_ballgown"),
        gtf = "results/stringtie2/{name1}.gtf",
    shell:
        """
        stringtie {input.bam} -p 10 -e -G {input.gtf} -A {output.quan} -o {output.gtf} -b {output.ballgown}
        """

#最后没有用组装的，用的自己添加的retain intron的jarid2的gtf，然后用stringtie定量的
rule stringtie_RNA_seq3:
    input:
        bam="../../1_expression_dynamics/core_scripts/results/mapping/filtered/{name1}.bam",
        gtf = "../../1_expression_dynamics/genome/GCF_001952655.1_M_albus_1.0_genomic_lgr4_jarid2modified.gtf"
    output:
        quan = "results/stringtie3/{name1}.txt",
        ballgown = directory("results/stringtie3/{name1}_ballgown"),
        gtf = "results/stringtie3/{name1}.gtf",
    shell:
        """
        stringtie {input.bam} -p 10 -e -G {input.gtf} -A {output.quan} -o {output.gtf} -b {output.ballgown}
        """

rule process:
    input:
        gtf = "results/stringtie3/{name1}.gtf"
    output:
        gtf = "results/stringtie3/modified/{name1}.gtf",
    shell:
        """
        python results/stringtie3/stringtie.py --gtf {input.gtf} --result {output.gtf}
        """



#cd ~/project/2_alternative_splicing_Whippet/core_scripts3_filter/results/stringtie3
#Rscript merge_stringtie2.R 
#cd ~/project/2_alternative_splicing_Whippet/core_scripts3_filter/results/stringtie3/modified
#cat mergerep.txt |grep "rna-XM_020597488" > jarid2_isoform.txt
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
        #expand("results/methylpy/bw/{name}_meth.bw",name = SAMPLE2),
        expand("results/stringtie/{name1}.txt",name1 = samplesRNASEQ),
        expand("results/gene_classify/{name2}/{name2}_high.bed12",name2=SAMPLE_RNA_seq),
        expand("results/methylpy/computematrix/{mergename}_results.gz",mergename=mergename),
        expand("results/methylpy/computematrix/{mergename}_results.pdf",mergename=mergename),
        expand("results/bismark/cpglevel/{sample}.txt",sample=SAMPLE),
        "results/bismark/cpglevel/merge.txt"
       


rule stats:
    input:
        "../results/bismark/processing/{sample}/{sample}_1_bismark_bt2_PE_report.txt"
    output:
        "results/bismark/cpglevel/{sample}.txt"
    shell:
        """
        cat {input} | grep "C methylated in CpG context:"|awk '{{print $6,"{wildcards.sample}"}}' > {output}
        """

rule merge_stats:
    input:
        expand("results/bismark/cpglevel/{sample}.txt",sample=SAMPLE)
    output:
        merge = "results/bismark/cpglevel/merge.txt",
    shell:
        """
        cat {input} > {output.merge}
        sed -i 's/%//g' {output.merge}
        """

rule stringtie_RNA_seq:
    input:
        bam="../../perl/1_expression_dynamics/core_scripts/results/mapping/filtered/{name1}.bam",
        gtf = "../../perl/1_expression_dynamics/genome/GCF_001952655.1_M_albus_1.0_genomic.gtf"
    output:
        quan = "results/stringtie/{name1}.txt",
        ballgown = directory("results/stringtie/{name1}_ballgown"),
        gtf = "results/stringtie/{name1}.gtf",
    shell:
        """
        stringtie {input.bam} -p 10 -e -G {input.gtf} -A {output.quan} -o {output.gtf} -b {output.ballgown}
        """


rule gene_classify:
    input:
        fpkm = "results/stringtie/mergerep.txt",
        bed = "../../perl/1_expression_dynamics/genome/Malbus.bed12",
        refid = "../../perl/3_alternative_polyadenylation/results/refid/ref_id.txt"
    output:
        high = "results/gene_classify/{name2}/{name2}_high.bed12",
        medium = "results/gene_classify/{name2}/{name2}_medium.bed12",
        low = "results/gene_classify/{name2}/{name2}_low.bed12"
    shell:
        """
        Rscript scripts/gene_classify.R {input.fpkm} {input.bed} {input.refid} {output.high} {output.medium} {output.low} {wildcards.name2}
        """

rule computematrix:
    input:
        bw ="../results/methylpy/bw/{name}_meth.bw",
        high = "results/gene_classify/{name2}/{name2}_high.bed12",
        medium = "results/gene_classify/{name2}/{name2}_medium.bed12",
        low = "results/gene_classify/{name2}/{name2}_low.bed12"
    output:
        gz = "results/methylpy/computematrix/{name}_vs_{name2}_results.gz",
        tab = "results/methylpy/computematrix/{name}__vs_{name2}_results.tab"
    shell:
        """
        computeMatrix scale-regions -R {input.high} {input.medium} {input.low} -S {input.bw} -b 5000 -a 5000 --regionBodyLength 8000 --binSize 200 --skipZeros -o {output.gz} --outFileNameMatrix {output.tab}
        """

rule plot:
    input:
        "results/methylpy/computematrix/{name}_vs_{name2}_results.gz",
    output:
        "results/methylpy/computematrix/{name}_vs_{name2}_results.pdf"
    shell:
        """
        plotProfile -m {input} -out {output} --plotHeight 10 --plotWidth 12 --plotFileFormat pdf --colors "#9c3123" "#e4bd64" "#436ca2" 
        """

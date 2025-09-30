#对于蛋白domain预测，首先根据bam文件传到IGV上的read看reads特征，然后制作JARID2的gff文件，并计算各个物种的IR ratio,然后gffread jarid2_IR.gff -g GCF_001858045.2_O_niloticus_UMD_NMBU_genomic.fna -S -y jarid2_ir.fa，在interproscan或者https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi看domain

#!/usr/bin/env runsnakemake
configfile: "config.yaml"
samples = config["samples"]
threads = config["threads"]

sampleORY= ["female","male"]
sampleorenil = ["female","male","intersex"]

rule all:
    input:
        expand("results/IRratio/malbus_yy/{sample}.txt",sample = samples),
        "results/IRratio/malbus_yy/IRratio.txt",
        expand("results/IRratio/Oryzias_latipes/{sample1}.txt",sample1 = sampleORY),
        expand("results/IRratio/Oreochromis_niloticus/{sample2}.txt",sample2 = sampleorenil)

rule calculate_IR_ratio:
    input:
        bam = "../../1_expression_dynamics/core_scripts/results/mapping/filtered/{sample}.bam"
    output:
        result = "results/IRratio/malbus_yy/{sample}.txt"
    shell:
        """
        python scripts/IR_ratio_yy.py --bam {input.bam} --result {output.result}
        """

rule merge:
    input:
        expand("results/IRratio/malbus_yy/{sample}.txt",sample = samples)
    output:
        result = "results/IRratio/malbus_yy/IRratio.txt"
    shell:
        """
        tail -n 1 {input} >> {output.result}
        #cat IRratio.txt | grep -v "txt"|grep -v '^\s*$' > IRratio_plot.txt
        """

rule calculate_IR_ratio_oryzia:
    input:
        bam = "Oryzias_latipes/results/mapping/filtered/{sample1}.bam"
    output:
        result = "results/IRratio/Oryzias_latipes/{sample1}.txt"
    shell:
        """
        python scripts/IR_ratio_oryzias.py --bam {input.bam} --result {output.result}
        """

rule calculate_IR_ratio_orenil:
    input:
        bam = "Oreochromis_niloticus/results/mapping/filtered/{sample2}.bam"
    output:
        result = "results/IRratio/Oreochromis_niloticus/{sample2}.txt"
    shell:
        """
        python scripts/IR_ratio_Orenil.py --bam {input.bam} --result {output.result}
        """

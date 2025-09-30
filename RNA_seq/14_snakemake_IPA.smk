#RNA5 conda activate IPAFinder
REPS = {"O_I":["../../1_expression_dynamics/core_scripts/results/mapping/filtered/O_I1.bam"],"O_IV":["../../1_expression_dynamics/core_scripts/results/mapping/filtered/O_IV1.bam","../../1_expression_dynamics/core_scripts/results/mapping/filtered/O_IV2.bam","../../1_expression_dynamics/core_scripts/results/mapping/filtered/O_IV3.bam"],"O_V":["../../1_expression_dynamics/core_scripts/results/mapping/filtered/O_V1.bam","../../1_expression_dynamics/core_scripts/results/mapping/filtered/O_V2.bam","../../1_expression_dynamics/core_scripts/results/mapping/filtered/O_V3.bam"],"OT_I":["../../1_expression_dynamics/core_scripts/results/mapping/filtered/OT_I1.bam","../../1_expression_dynamics/core_scripts/results/mapping/filtered/OT_I2.bam"],"OT_II":["../../1_expression_dynamics/core_scripts/results/mapping/filtered/OT_II1.bam","../../1_expression_dynamics/core_scripts/results/mapping/filtered/OT_II2.bam"],"OT_III":["../../1_expression_dynamics/core_scripts/results/mapping/filtered/OT_III1.bam","../../1_expression_dynamics/core_scripts/results/mapping/filtered/OT_III2.bam","../../1_expression_dynamics/core_scripts/results/mapping/filtered/OT_III3.bam"],"T_I":["../../1_expression_dynamics/core_scripts/results/mapping/filtered/T_I1.bam","../../1_expression_dynamics/core_scripts/results/mapping/filtered/T_I2.bam","../../1_expression_dynamics/core_scripts/results/mapping/filtered/T_I3.bam"],"T_III":["../../1_expression_dynamics/core_scripts/results/mapping/filtered/T_III1.bam","../../1_expression_dynamics/core_scripts/results/mapping/filtered/T_III2.bam"],"T_IV":["../../1_expression_dynamics/core_scripts/results/mapping/filtered/T_IV1.bam","../../1_expression_dynamics/core_scripts/results/mapping/filtered/T_IV2.bam","../../1_expression_dynamics/core_scripts/results/mapping/filtered/T_IV3.bam"]}

SAMPLE = ["O_I","O_IV","O_V","OT_I","OT_II","OT_III","T_I","T_III","T_IV"]
#SAMPLE = ["O_I","OT_II"]
diff = []

for i in range(0,len(SAMPLE)-1):
    for j in range(i+1,len(SAMPLE)):
        diff.append(SAMPLE[i] + '_vs_' + SAMPLE[j])

rule all:
    input:
        "results/IPA/annotation/IPAFinder_anno_malbus.txt",
        "results/IPA/annotation/IPAFinder_anno_malbus_filter.txt",
        expand("results/IPA/DetectIPA/{diffname}_IPAFinder_IPUI.txt",diffname=diff),
        expand("results/IPA/configure/{diffname}.txt",diffname=diff),
        "results/IPA/DetectIPA/all_IPAFinder_IPUI.txt",
        expand("results/IPA/DUIPA/{diffname}_IPAFinder_DUIPA.txt",diffname=diff),
        "results/IPA/filter/results.txt",
        "results/IPA/filter/interate.txt",
        expand("results/IPA/filtergene/{diffname}/passfiltergene.txt",diffname=diff),
        expand("results/IPA/skipped_filter/{diffname}_passfiltergene.txt",diffname=diff),
        expand("results/IPA/expression_filter/{diffname}_passfiltergene.txt",diffname=diff),
        expand("results/IPA/filterfinal/{diffname}/passfilter.txt",diffname=diff),
        "results/IPA/filterfinal/interate.txt",
        "results/IPA/expression_filter/IPA_zebrafish_blat.txt"

#touch results/IPA/annotation/IPAFinder_anno_malbus.txt
rule annotation:
    input:
        gtf = "../../1_expression_dynamics/genome/GCF_001952655.1_M_albus_1.0_genomic_noNC.gtf"
    output:
        result = "results/IPA/annotation/IPAFinder_anno_malbus.txt"
    shell:
        """
        python ../IPAFinder/IPAFinder/IPAFinder_GetAnno.py -gtf {input.gtf} -output {output.result}
        """ 

rule filter_annotation:
    input:
        "results/IPA/annotation/IPAFinder_anno_malbus.txt"
    output:
        "results/IPA/annotation/IPAFinder_anno_malbus_filter.txt"
    shell:
        """
        python ../scripts/filter_annotation.py --anno {input} --result {output}
        """
#/home/lyy/miniconda3/envs/IPAFinder/bin/python ../IPAFinder/IPAFinder/IPAFinder_DetectIPA_modified2.py  -b {output.configure} -anno {input.anno} -p 28 -o {output.results} -output_directory {output.dir1}
# /home/lyy/miniconda3/envs/IPAFinder/bin/python ../IPAFinder/IPAFinder/IPAFinder_DetectIPA_modified2.py  -b {output.configure} -anno {input.anno} -p 28 -o {output.results} -output_directory {output.dir1}
rule DetectIPA:
    input:
        #anno = "results/IPA/annotation/test200000.txt",
        anno = "results/IPA/annotation/IPAFinder_anno_malbus_filter.txt",
        bam1 =lambda wildcards: REPS[wildcards.sample1],
        bam2 =lambda wildcards: REPS[wildcards.sample2]
    output:
        configure = "results/IPA/configure/{sample1}_vs_{sample2}.txt",
        results = "results/IPA/DetectIPA/{sample1}_vs_{sample2}_IPAFinder_IPUI.txt",
        #dir1 = directory("{sample1}_vs_{sample2}/dir")
    shell:
        """
        echo -e 'condition1={input.bam1}\ncondition2={input.bam2}' >{output.configure}
        /home/lyy/miniconda3/envs/IPAFinder/bin/python ../IPAFinder/IPAFinder/IPAFinder_DetectIPA_modified.py  -b {output.configure} -anno {input.anno} -p 20 -o {output.results} -output_directory {wildcards.sample1}_vs_{wildcards.sample2}/
        """

rule DetectallIPA:
    input:
        #anno = "results/IPA/annotation/test200000.txt",
        anno = "results/IPA/annotation/IPAFinder_anno_malbus_filter.txt",
        configure = "results/IPA/configure/all.txt",
    output:
        results = "results/IPA/DetectIPA/all_IPAFinder_IPUI.txt",
        #dir1 = directory("{sample1}_vs_{sample2}/dir")
    shell:
        """
        /home/lyy/miniconda3/envs/IPAFinder/bin/python ../IPAFinder/IPAFinder/IPAFinder_DetectIPA_modified.py  -b {input.configure} -anno {input.anno} -p 50 -o {output.results} -output_directory all/
        """

rule DUIPA:
    input:
        IPUI = "results/IPA/DetectIPA/{sample1}_vs_{sample2}_IPAFinder_IPUI.txt",
        configure = "results/IPA/configure/{sample1}_vs_{sample2}.txt",
    output:
        results = "results/IPA/DUIPA/{sample1}_vs_{sample2}_IPAFinder_DUIPA.txt"
    shell:
        """
        Rscript ../IPAFinder/IPAFinder/Infer_DUIPA_modified2.R -b {input.configure} -I {input.IPUI} -d {wildcards.sample1}_vs_{wildcards.sample2}/ -o {output.results}
        """

rule filterpass:
    input:
        alldata="results/IPA/DUIPA/{sample1}_vs_{sample2}_IPAFinder_DUIPA.txt"
    output:
        passfilter="results/IPA/filter/{sample1}_vs_{sample2}/passfilter.txt"
    shell:
        """
        Rscript ../scripts/valcano_IPA.R {input.alldata} {wildcards.sample1} {wildcards.sample2} {output.passfilter}
        """

rule mergetxt:
    input:
        passfilter=expand("results/IPA/filter/{diffname}/passfilter.txt",diffname=diff)
    output:
        results="results/IPA/filter/results.txt"
    shell:
        """
        cat {input.passfilter} > {output.results}
        """

rule integrate:
    input:
        "results/IPA/filter/results.txt"
    output:
        "results/IPA/filter/interate.txt"
    shell:
        """
        python ../scripts/intergrate_IPA.py --merge {input} --result {output}
        """


rule filterpassgene:
    input:
        alldata="results/IPA/DUIPA/{sample1}_vs_{sample2}_IPAFinder_DUIPA.txt",
    output:
        passfiltergene="results/IPA/filtergene/{sample1}_vs_{sample2}/passfiltergene.txt"
    shell:
        """
        Rscript ../scripts/allIPA.R {input.alldata} {wildcards.sample1} {wildcards.sample2} {output.passfiltergene}
        """

rule skipped_filter:
    input:
        passfiltergene="results/IPA/filtergene/{sample1}_vs_{sample2}/passfiltergene.txt",
        anno = "results/IPA/annotation/IPAFinder_anno_malbus_filter.txt",
        bam1 = lambda wildcards: REPS[wildcards.sample1],
        bam2 =lambda wildcards: REPS[wildcards.sample2],
    output:
        "results/IPA/skipped_filter/{sample1}_vs_{sample2}_passfiltergene.txt"
    shell:
        """
        python ../scripts/IPA_skipped_filter.py --bam1 {input.bam1} --bam2 {input.bam2} --IPA {input.passfiltergene} --anno {input.anno} --result {output}
        """

rule expression_filter:
    input:
        passfiltergene="results/IPA/skipped_filter/{sample1}_vs_{sample2}_passfiltergene.txt",
        anno = "results/IPA/annotation/IPAFinder_anno_malbus_filter.txt",
        gtf = "../../1_expression_dynamics/genome/GCF_001952655.1_M_albus_1.0_genomic_noNC.gtf",
        bam1 = lambda wildcards: REPS[wildcards.sample1],
        bam2 =lambda wildcards: REPS[wildcards.sample2],
    output:
        result =  "results/IPA/expression_filter/{sample1}_vs_{sample2}_passfiltergene.txt",
        result2 = "results/IPA/expression_foldchange/{sample1}_vs_{sample2}_passfiltergene.txt"
    shell:
        """
        python ../scripts/exon_expression_filter.py  --bam1 {input.bam1} --bam2 {input.bam2} --IPA {input.passfiltergene} --gtf {input.gtf} --anno {input.anno} --result {output.result} --result2 {output.result2}
        """


rule filterpass2:
    input:
        alldata="results/IPA/expression_filter/{sample1}_vs_{sample2}_passfiltergene.txt"
    output:
        passfilter="results/IPA/filterfinal/{sample1}_vs_{sample2}/passfilter.txt"
    shell:
        """
        Rscript ../scripts/valcano_IPAfinal.R {input.alldata} {wildcards.sample1} {wildcards.sample2} {output.passfilter}
        """

rule mergetxt2:
    input:
        passfilter=expand("results/IPA/filterfinal/{diffname}/passfilter.txt",diffname=diff)
    output:
        results="results/IPA/filterfinal/results.txt"
    shell:
        """
        cat {input.passfilter} > {output.results}
        """

rule integrate2:
    input:
        "results/IPA/filterfinal/results.txt"
    output:
        "results/IPA/filterfinal/interate.txt"
    shell:
        """
        python ../scripts/intergrate_IPA.py --merge {input} --result {output}
        """

rule ortho_to_zebrafish:
    input:
        "results/IPA/expression_filter/IPA.txt"
    output:
        mapfile = "results/IPA/expression_filter/IPA_zebrafish.txt",
    shell:
        """
        python ../../2_alternative_splicing_Whippet/scripts/psl_v2.py --psl {input} --result {output.mapfile} --ortho ../../2_alternative_splicing_Whippet/Ortho/Monopterus_albus__v__Zebrafish.tsv --s1 Monopterus_albus --s2 Zebrafish
        """

rule ortho_to_zebrafish_add_blat:
    input:
        mapfile = "results/IPA/expression_filter/IPA_zebrafish.txt",
        blatfile = "../../2_alternative_splicing_Whippet/Ortho/blat_best_map_to_zebrafish.txt"
    output:
        mapfile2 = "results/IPA/expression_filter/IPA_result_zebrafish_blat.txt",
        cut = "results/IPA/expression_filter/IPA_zebrafish_blat.txt"
    shell:
        """
        python ../../2_alternative_splicing_Whippet/scripts/addblat.py --orthoresult {input.mapfile} --blatresult {input.blatfile} --result {output.mapfile2}
        cat {output.mapfile2} | awk '{{print $2}}' | sort | uniq > {output.cut}
        """


SAMPLE=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]
diff2 = []
for i in range(0,len(SAMPLE)-1):
    diff2.append(SAMPLE[i] + '_vs_' + SAMPLE[i+1])
#print(len(diff2))
GO = ["female_gamete_generation","sex_change","sex_determination","sex_differentiation","spermatogenesis","sex_reversal","oogenesis","GO_sex_related"]


rule all:
    input:
        expand("results/tf_filter/{diffname}_filter.txt",diffname = diff2),
        "results/merge/sex_gene.txt"
         

rule get_tf:
    input:
        "../core_scripts/results/ananse_influence/{diffname}_influence.txt"
    output:
        "results/tf_filter/{diffname}_filter.txt"
    shell:
        """
        cat {input} | awk  '$2>=0.8{{print $0}}' > {output}
        """

rule merge:
    input:
        expand("results/tf_filter/{diffname}_filter.txt",diffname = diff2)
    output:
        file1 = "results/merge/merge.txt",
        file2 = "results/merge/merge_sort.txt",
    shell:
        """
        python ../scripts/merge.py --influence {input} --result {output}
        sort -k10,10nr {output.file1} > {output.file2}
        """

rule ortho_to_zebrafish:
    input:
        "results/merge/merge.txt",
    output:
        file1 = "results/merge/merge_name.txt",
        mapfile = "results/merge/merge_name_zebrafish.txt",
    shell:
        """
        cat {input} | cut -f 1 > {output.file1}
        python ../core_scripts/Ortho/psl_v2.py --psl {output.file1} --result {output.mapfile} --ortho ../core_scripts/Ortho/Monopterus_albus__v__Zebrafish.tsv --s1 Monopterus_albus --s2 Zebrafish
        """

rule ortho_to_zebrafish_add_blat:
    input:
        mapfile = "results/merge/merge_name_zebrafish.txt",
        blatfile = "../core_scripts/Ortho/blat_best_map_to_zebrafish.txt"
    output:
        mapfile2 = "results/merge/merge_name_zebrafish_blat.txt",
    shell:
        """
        python ../core_scripts/Ortho/addblat.py --orthoresult {input.mapfile} --blatresult {input.blatfile} --result {output.mapfile2}
        """

rule add_GO:
    input:
        mapfile = "results/merge/merge_name_zebrafish_blat.txt",
        go = expand("../core_scripts/GO/{name2}.txt",name2=GO),
    output:
        "results/merge/merge_annotation.txt"
    shell:
        """
        python ../core_scripts/scripts/GO_annoation.py --tf {input.mapfile} --go {input.go} --result {output}
        """

rule get_sex_gene:
    input:
        "results/merge/merge_annotation.txt"
    output:
        file1 = "results/merge/sex_gene.txt",
        file2 = "results/merge/sex_gene_modified.txt"
    shell:
        """
        cat {input}|grep "GO_sex_related"|cut -f 9,10 > {output.file1}
        python ../core_scripts/scripts/sex_gene_modified.py --file1 {output.file1} --result {output.file2}
        """


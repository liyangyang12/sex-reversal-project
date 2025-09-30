SAMPLE=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]
diff2 = []
for i in range(0,len(SAMPLE)-1):
    diff2.append(SAMPLE[i] + '_vs_' + SAMPLE[i+1])
print(len(diff2))
GO = ["female_gamete_generation","sex_change","sex_determination","sex_differentiation","spermatogenesis","sex_reversal","oogenesis","GO_sex_related"]

rule all:
    input:
        expand("results/ananse_influence_parse/{diffname}_annotation.txt",diffname = diff2),
        

rule get_name:
    input:
        "results/ananse_influence/{diffname}_influence.txt"
    output:
        "results/ananse_influence_parse/{diffname}_name.txt"
    shell:
        """
        cat {input} | cut -f 1 > {output}
        """

rule ortho_to_zebrafish:
    input:
        "results/ananse_influence_parse/{diffname}_name.txt",
    output:
        mapfile = "results/ananse_influence_parse/{diffname}_name_zebrafish.txt",
    shell:
        """
        python Ortho/psl_v2.py --psl {input} --result {output.mapfile} --ortho Ortho/Monopterus_albus__v__Zebrafish.tsv --s1 Monopterus_albus --s2 Zebrafish
        """

rule ortho_to_zebrafish_add_blat:
    input:
        mapfile = "results/ananse_influence_parse/{diffname}_name_zebrafish.txt",
        blatfile = "Ortho/blat_best_map_to_zebrafish.txt"
    output:
        mapfile2 = "results/ananse_influence_parse/{diffname}_name_zebrafish_blat.txt",
    shell:
        """
        python Ortho/addblat.py --orthoresult {input.mapfile} --blatresult {input.blatfile} --result {output.mapfile2}
        """

rule add_GO:
    input:
        mapfile = "results/ananse_influence_parse/{diffname}_name_zebrafish_blat.txt",
        go = expand("GO/{name2}.txt",name2=GO),
    output:
        "results/ananse_influence_parse/{diffname}_annotation.txt"
    shell:
        """
        python scripts/GO_annoation.py --tf {input.mapfile} --go {input.go} --result {output}
        """
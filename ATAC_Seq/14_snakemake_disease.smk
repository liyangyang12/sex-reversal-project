rule all:
    input:
        "results/merge/blat_best_map_to_human.txt",
        "results/disease/merge_name_human_blat.txt",
        "results/disease/disease_merge.txt"


rule get_best_blat:
    input:    
        "../core_scripts/Ortho/Monopterus_albus_Homo_sapiens_annote.result1"
    output:
        "results/merge/blat_best_map_to_human.txt"
    shell:
        """
        python ../core_scripts/Ortho/parse_blat_psl.py --blat {input} --result {output}
        """


rule ortho_to_human:
    input:
        "results/merge/merge_sort.txt",
    output:
        file1 = "results/disease/merge_name.txt",
        mapfile = "results/disease/merge_name_human.txt",
    shell:
        """
        cat {input} | cut -f 1 > {output.file1}
        python ../core_scripts/Ortho/psl_v2.py --psl {output.file1} --result {output.mapfile} --ortho ../core_scripts/Ortho/Monopterus_albus__v__Human.tsv --s1 Monopterus_albus --s2 Human
        """

rule ortho_to_human_add_blat:
    input:
        mapfile = "results/disease/merge_name_human.txt",
        blatfile = "results/merge/blat_best_map_to_human.txt"
    output:
        mapfile2 = "results/disease/merge_name_human_blat.txt",
    shell:
        """
        python ../core_scripts/Ortho/addblat.py --orthoresult {input.mapfile} --blatresult {input.blatfile} --result {output.mapfile2}
        """

rule disease_merge:
    input:
        "results/disease/disease_score.txt"
    output:
        "results/disease/disease_merge.txt"
    shell:
        """
        python merge.py --disease {input} --result {output}
        """
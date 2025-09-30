rule all:
    input:
        "spliceosome/spliceosome_all_zebrafish_blat_filter.txt"


#cat spliceosome.txt |grep -v "GeneID" |cut -f 3 >spliceosome_all.txt
rule ortho_to_malbus:
    input:
        "spliceosome/spliceosome_all.txt"
    output:
        mapfile = "spliceosome/spliceosome_all_zebrafish.txt",
    shell:
        """
        python spliceosome/psl_v2.py --psl {input} --result {output.mapfile} --ortho spliceosome/Monopterus_albus__v__Human.tsv --s1 Human --s2 Monopterus_albus
        """

rule ortho_to_malbus_add_blat:
    input:
        mapfile = "spliceosome/spliceosome_all_zebrafish.txt",
        blatfile = "spliceosome/blat_best_map_from_malbus_to_Homo_sapiens.txt"
    output:
        mapfile2 = "spliceosome/spliceosome_all_zebrafish_blat.txt",
    shell:
        """
        python spliceosome/addblat.py --orthoresult {input.mapfile} --blatresult {input.blatfile} --result {output.mapfile2}
        """

rule filter:
    input:
        "spliceosome/spliceosome_all_zebrafish_blat.txt"
    output:
        "spliceosome/spliceosome_all_zebrafish_blat_filter.txt"
    shell:
        """
        python spliceosome/filter.py --orthoresult {input} --result {output}
        """
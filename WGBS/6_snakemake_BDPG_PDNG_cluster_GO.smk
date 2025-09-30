cluster = ["c1","c2","c3"]
rule all:
    input:
        expand("results/BDPG_PDNG_cluster/{cluster}_zebrafish_blat.txt",cluster=cluster)



rule ortho_to_zebrafish:
    input:
        "results/BDPG_PDNG_cluster/{cluster}.txt"
    output:
        mapfile = "results/BDPG_PDNG_cluster/{cluster}_zebrafish.txt",
    shell:
        """
        python ../../perl/2_alternative_splicing_Whippet/scripts/psl_v2.py --psl {input} --result {output.mapfile} --ortho ../../perl/2_alternative_splicing_Whippet/Ortho/Monopterus_albus__v__Zebrafish.tsv --s1 Monopterus_albus --s2 Zebrafish
        """

rule ortho_to_zebrafish_add_blat:
    input:
        mapfile ="results/BDPG_PDNG_cluster/{cluster}_zebrafish.txt",
        blatfile = "../../perl/2_alternative_splicing_Whippet/Ortho/blat_best_map_to_zebrafish.txt"
    output:
        mapfile2 = "results/BDPG_PDNG_cluster/{cluster}_result_zebrafish_blat.txt",
        cut = "results/BDPG_PDNG_cluster/{cluster}_zebrafish_blat.txt"
    shell:
        """
        python ../../perl/2_alternative_splicing_Whippet/scripts/addblat.py --orthoresult {input.mapfile} --blatresult {input.blatfile} --result {output.mapfile2}
        cat {output.mapfile2} | awk '{{print $2}}' | sort | uniq > {output.cut}
        """

loc = ["promoter","genebody"]
cor = ["positive","negative"]
rule all:
    input:
        expand("results/Cor_GO/Cor_{loc}_{cor}_zebrafish_blat.txt",loc=loc,cor=cor)


rule get_gene:
    input:
        gene1 = "results/cor_DMR_methlevel_exp/{loc}_{cor}.txt",
    output:
        "results/Cor_GO/Cor_{loc}_{cor}.txt"
    shell:
        """
        cat {input.gene1} | grep -v "NA" | cut -f 1  | sort |uniq > {output}
        """

rule ortho_to_zebrafish:
    input:
        "results/Cor_GO/Cor_{loc}_{cor}.txt"
    output:
        mapfile = "results/Cor_GO/Cor_{loc}_{cor}_zebrafish.txt",
    shell:
        """
        python ../../perl/2_alternative_splicing_Whippet/scripts/psl_v2.py --psl {input} --result {output.mapfile} --ortho ../../perl/2_alternative_splicing_Whippet/Ortho/Monopterus_albus__v__Zebrafish.tsv --s1 Monopterus_albus --s2 Zebrafish
        """

rule ortho_to_zebrafish_add_blat:
    input:
        mapfile ="results/Cor_GO/Cor_{loc}_{cor}_zebrafish.txt",
        blatfile = "../../perl/2_alternative_splicing_Whippet/Ortho/blat_best_map_to_zebrafish.txt"
    output:
        mapfile2 = "results/Cor_GO/Cor_{loc}_{cor}_result_zebrafish_blat.txt",
        cut = "results/Cor_GO/Cor_{loc}_{cor}_zebrafish_blat.txt"
    shell:
        """
        python ../../perl/2_alternative_splicing_Whippet/scripts/addblat.py --orthoresult {input.mapfile} --blatresult {input.blatfile} --result {output.mapfile2}
        cat {output.mapfile2} | awk '{{print $2}}' | sort | uniq > {output.cut}
        """

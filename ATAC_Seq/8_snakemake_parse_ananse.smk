SAMPLE=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]
GO = ["female_gamete_generation","sex_change","sex_determination","sex_differentiation","spermatogenesis","sex_reversal","oogenesis"]

rule all:
    input:
        expand("results/ananse_network_filter/{sample}_sort_network.tsv",sample = SAMPLE),
        expand("results/ananse_network_filter/{sample}_top_network.tsv",sample = SAMPLE),
        expand("results/tf_degree/{sample}_name_zebrafish_blat.txt",sample = SAMPLE),
        expand("results/tf_degree/{sample}_annotation.txt",sample = SAMPLE)


rule filter_network:
    input:
        "results/ananse_network/{sample}_network.tsv"
    output:
        "results/ananse_network_filter/{sample}_sort_network.tsv"
    shell:
        """
        sort -k 2r {input}  > {output}
        """

rule get_gene:
    input:
        "results/ananse_network_filter/{sample}_sort_network.tsv"
    output:
        "results/ananse_network_filter/{sample}_top_network.tsv"
    shell:
        """
        #head -n 100000 {input} >{output}
        cat {input} | awk  '$2>=0.8{{print $0}}' > {output}
        """

rule tf_degree:
    input:
        "results/ananse_network_filter/{sample}_top_network.tsv"
    output:
        count = "results/tf_degree/{sample}_count.txt",
        network = "results/tf_degree/{sample}_network.txt",
        sortcount = "results/tf_degree/{sample}_sortcount.txt",
        name = "results/tf_degree/{sample}_name.txt",
    shell:
        """
        python scripts/TF_count.py --network {input} --resultcount {output.count} --resultnetwork {output.network}
        sort -k2,2nr {output.count} > {output.sortcount}
        cat {output.sortcount} | cut -f 1 > {output.name}
        """

rule ortho_to_zebrafish:
    input:
        "results/tf_degree/{sample}_name.txt",
    output:
        mapfile = "results/tf_degree/{sample}_name_zebrafish.txt",
    shell:
        """
        python Ortho/psl_v2.py --psl {input} --result {output.mapfile} --ortho Ortho/Monopterus_albus__v__Zebrafish.tsv --s1 Monopterus_albus --s2 Zebrafish
        """

rule ortho_to_zebrafish_add_blat:
    input:
        mapfile = "results/tf_degree/{sample}_name_zebrafish.txt",
        blatfile = "Ortho/blat_best_map_to_zebrafish.txt"
    output:
        mapfile2 = "results/tf_degree/{sample}_name_zebrafish_blat.txt",
    shell:
        """
        python Ortho/addblat.py --orthoresult {input.mapfile} --blatresult {input.blatfile} --result {output.mapfile2}
        """

rule add_GO:
    input:
        mapfile = "results/tf_degree/{sample}_name_zebrafish_blat.txt",
        go = expand("GO/{name2}.txt",name2=GO),
    output:
        "results/tf_degree/{sample}_annotation.txt"
    shell:
        """
        python scripts/GO_annoation_v1.py --tf {input.mapfile} --go {input.go} --result {output}
        """

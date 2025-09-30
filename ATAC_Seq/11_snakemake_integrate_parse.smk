SAMPLE=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]
GO = ["female_gamete_generation","sex_change","sex_determination","sex_differentiation","spermatogenesis","sex_reversal","oogenesis"]

#SAMPLE2=["VO","MO"]
#SAMPLE2=["ISE","ISM"]

SAMPLE2=["ES","MS"]
SAMPLE3 = ["VO","MO","ISE","ISM","ES","MS"]
diff2 = ["VO_vs_MO","ISE_vs_ISM","ES_vs_MS"]

#print(len(diff2))

rule all:
    input:
        expand("results/extract_top_tf_network/{diffname}.txt",diffname = diff2),
        expand("results/extract_top_tf_network/{sample2}_network.txt",sample2 = SAMPLE2),
        expand("results/extract_top_tf_network/{sample2}_sortcount.txt",sample2 = SAMPLE3)

rule extracttf:
    input:
        influence = "../core_scripts/results/ananse_influence/{diffname}_influence.txt",
    output:
        filte = "results/extract_top_tf_network/{diffname}.txt"
    shell:
        """
        cat {input.influence} |sed -n '2,16p' |cut -f 1 > {output.filte}
        """

rule extract_network:
    input:
        #filte = "results/extract_top_tf_network/VO_vs_MO.txt",
        #filte = "results/extract_top_tf_network/ISE_vs_ISM.txt",
        filte = "results/extract_top_tf_network/ES_vs_MS.txt",
        net = "../core_scripts/results/network_integrate/{sample2}_GO.txt"
    output:
        enet = "results/extract_top_tf_network/{sample2}_network.txt",
    shell:
        """
        python ../core_scripts/scripts/filter_network.py --filte {input.filte} --net {input.net} --result {output.enet}
        """
        
rule extract_network_degree:
    input:
        "results/extract_top_tf_network/{sample2}_network.txt"
    output:
        count = "results/extract_top_tf_network/{sample2}_count.txt",
        sortcount = "results/extract_top_tf_network/{sample2}_sortcount.txt",
    shell:
        """
        python ../core_scripts/scripts/TF_count_integrate.py --network {input} --resultcount {output.count}
        sort -k2,2nr {output.count} > {output.sortcount}
        """


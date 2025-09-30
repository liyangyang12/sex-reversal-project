SAMPLE=["MO","ISM","MS"]
species = ["Homo_sapiens","Mus_musculus"]


rule all:
    input:
        #expand("results/best_map/blat_best_map_from_{species}_to_malbus.txt",species = species),
        expand("results/TRRUST_prepare/{sample}_and_{species}_only.txt",sample=SAMPLE,species = species),
        expand("results/TRRUST_compare/{sample}_and_{species}_only.txt",sample=SAMPLE,species = species)
        

rule get_best_blat:
    input:    
        "blat/Monopterus_albus_{species}_annote.result1"
    output:
        "results/best_map/blat_best_map_from_{species}_to_malbus.txt"
    shell:
        """
        python scripts/parse_blat_psl.py --blat {input} --result {output}
        """


rule modify:
    input:
        enet = "results/extract_top_tf_network/{sample}_network.txt",
        bestmap = "../core_scripts/results/best_map/blat_best_map_from_{species}_to_malbus.txt",
        Ortho = "../core_scripts/Ortho/Orthogroups_modified.tsv"
    output:
        result1 = "results/TRRUST_prepare/{sample}_and_{species}_only.txt",
        result2 = "results/TRRUST_prepare/{sample}_and_{species}_multi.txt",
    shell:
        """
        python ../core_scripts/scripts/modified_TRRUST.py --network {input.enet} --ortho {input.Ortho} --bestmap {input.bestmap} --species {wildcards.species} --result {output.result1} --multi {output.result2}
        """


rule compare:
    input:
        trrust = "../core_scripts/TRRUST/{species}.tsv",
        net = "results/TRRUST_prepare/{sample}_and_{species}_only.txt",
    output:
        compare = "results/TRRUST_compare/{sample}_and_{species}_only.txt"
    shell:
        """
        python ../core_scripts/scripts/compare.py --trrust {input.trrust} --net {input.net} --result {output.compare}
        """

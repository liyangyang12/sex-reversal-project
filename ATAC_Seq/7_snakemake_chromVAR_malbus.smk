rule all:
        input:
            "results/malbus/zscore_malbus.txt",
            "results/malbus/zscore_DEG_malbus.txt"

rule ortho_to_zebrafish:
    input:
        "zscore.txt"
    output:
        mapfile = "results/malbus/zscore_malbus.txt",
    shell:
        """
        python ../ortho/psl2.py --psl {input} --result {output.mapfile} --ortho ../ortho/Orthogroups.tsv --s1 Monopterus_albus --s2 Danio_rerio --s3 Homo_sapiens --s4 Mus_musculus
        """




rule ortho_to_zebrafish_DEG:
    input:
        "zscore_DEG.txt"
    output:
        mapfile = "results/malbus/zscore_DEG_malbus.txt",
    shell:
        """
        python ../ortho/psl2.py --psl {input} --result {output.mapfile} --ortho ../ortho/Orthogroups.tsv --s1 Monopterus_albus --s2 Danio_rerio --s3 Homo_sapiens --s4 Mus_musculus
        """


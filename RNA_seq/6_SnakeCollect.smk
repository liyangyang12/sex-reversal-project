#!/usr/bin/env snakemake

names = [
    #"stringtie_stringtie", 
    "stringtie_taco",
    # "cufflink_cuffmerge",
    # "cufflink_taco"
]

rule all:
    input:
        expand("results/collect/{name}.gtf.gz", name=names),


def get_input_gtf(name):
    if name == "stringtie_taco":
        return "results/stringtie/taco/stringtie_taco.sorted.gtf.gz"
    #elif name == "stringtie_taco":
    #    return "results/stringtie/taco/stringtie_taco.sorted.gtf.gz"
    # elif name == "cufflink_cuffmerge":
    #     return "results/cufflinks/merged.filtered.gtf"
    # elif name == "cufflink_taco":
    #     return "results/cufflinks/taco/merged/assembly.gtf"
    else:
        print(name)
        assert False


rule collect:
    input:
        gtf = lambda wildcards: get_input_gtf(wildcards.name)
    output:
        gtf = "results/collect/{name}.gtf.gz",
        tbi = "results/collect/{name}.gtf.gz.tbi"
    shell:
        """
        cp {input.gtf} {output.gtf}
        tabix -p gff {output.gtf}
        """


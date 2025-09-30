SAMPLE2=["PVO","VO","MO","ISE","ISM","ISL","ES","MS","LS"]

regiom = ["fiveUTR","exon","intron","threeUTR","cpgIsland"]
#regiom = ["cpgIsland"]

rule all:
    input:
        expand("results/CGI/computematrix2/{region}_results.gz",region=regiom),
        expand("results/CGI/computematrix2/pdf/{region}_results.pdf",region=regiom)


rule cutbed:
    input:
        "regionbed/{region}.bed"
    output:
        "regionbed/{region}_cut.bed"
    shell:
        """
        cut -f 3,4,5,7 {input} > {output}
        """


rule computematrix:
    input:
        bw =expand("../results/methylpy/bw/{name}_meth.bw",name = SAMPLE2),
        bed = "regionbed/{region}_cut.bed",
    output:
        gz = "results/CGI/computematrix2/{region}_results.gz",
        tab = "results/CGI/computematrix2/{region}_results.tab"
    shell:
        """
        computeMatrix scale-regions -p 30 -R {input.bed} -S {input.bw} -b 0 -a 0 --skipZeros -o {output.gz} --outFileNameMatrix {output.tab}
        """


rule plot:
    input:
        "results/CGI/computematrix2/{region}_results.gz",
    output:
        "results/CGI/computematrix2/pdf/{region}_results.pdf"
    shell:
        """
        plotProfile -m {input} -out {output} --perGroup --plotHeight 10 --plotWidth 5 --yMin 0 --yMax 1  --plotFileFormat pdf --colors "#a23124" "#cb3228" "#e86c44" "#f4aa62" "#ecc264" "#a9d3e2" "#72a6c8" "#4470aa" "#1d2e62"
        """


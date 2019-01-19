import os
import tempfile

#
# A little complicated to find the temp dir
#
SSD_TMP_DIR = "/data/scratch/ssd"
if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
elif os.path.exists(SSD_TMP_DIR):
    TMPDIR = SSD_TMP_DIR
else:
    TMPDIR = tempfile.gettempdir()

configfile: "phasedsv.json"
SD = os.path.dirname(workflow.snakefile)


rule all:
    input:
        phasedSites     = "phased_regions.bed",
        phasedRegions   = "PhasedRegions.bed",
        regionsbed      = "Regions.bed",
        regionsasm      = "Regions.asm"

# vcffilter -g "( GT = 1|0 | GT = 0|1 ) "  | awk '{{ if ($NF != ".") print;}}' | \

rule MakePhasedSites:
    input:
        vcf=config["vcf"]
    output:
        sites="phased_regions.bed"
    params:
        sample=config["sample"]
    shell:"""
vcftools --gzvcf {input.vcf} --indv {params.sample} --recode  --stdout | \
awk '{{ if (length($4) ==1 && length($5) ==1 && $substr($0,0,1) == "#" || substr($10,0,3) == "0|1" || substr($10,0,3) == "1|0") print; }}' | \
 bedtools cluster -i stdin -d 5000 | \
 bedtools groupby -g 11 -c 11,1,2,2 -o count,first,first,last | \
awk '{{ print $3"\\t"$4"\\t"$5"\\t"$5-$4"\\t"$2;}}' > {output.sites}
"""

rule MergePhasedSites:
    input:
        sites="phased_regions.bed"
    output:
        merged="PhasedRegions.bed"
    params:
        ref=config["ref"]
    shell:"""
cat {input.sites} | \
  awk '{{ if ($4 >= 1000 && $5 > 2) print;}}' | \
  bedtools slop -g {params.ref}.fai -b 7000 | \
  bedtools sort | \
  awk 'BEGIN{{ pc="NONE"; prev=-1000000;}} {{ if ($1 == $pc && $2 - $prev < 10000) {{ $2=prev-1;}} prev=$3; pc=$1; print $1"\\t"$2"\\t"$3;}}' | \
   awk '{{ if ($3-$2 > 20000) print;}}' | \
   bedtools merge > PhasedRegions.bed
"""

rule SplitByZygosity:
    input:
        phased="PhasedRegions.bed"
    output:
        autozygous="Autozygous.bed"
    params:
        sd=SD
    shell:"""
bedtools subtract -a {params.sd}/HighComplexity.bed -b {input.phased} > {output.autozygous}
"""

rule MakeRegions:
    input:
        het="PhasedRegions.bed",
        auto="Autozygous.bed"
    output:
        regions="Regions.bed"
    params:
        window=config["window"],
        step=config["step"],
        ref=config["ref"]
    shell:"""
bedtools slop -b 5000 -g {params.ref}.fai -i {input.het} | \
bedtools makewindows -b stdin -w {params.window} -s {params.step} | \
awk '{{ if ($3!= prev) {{ print;}} prev=$3;}}' | \
awk '{{ print $0"\\tHET";}}' > {input.het}.rgn

bedtools makewindows -b {input.auto} -w {params.window} -s {params.step} | \
awk '{{ if ($3 != prev) {{ print;}} prev=$3;}}' | \
 awk '{{ print $0"\\tHOM";}}' > {input.auto}.rgn
cat {input.het}.rgn {input.auto}.rgn | bedtools sort > {output.regions}
"""

rule MakeAsm:
    input:
        bed = "Regions.bed"
    output:
        asm = "Regions.asm"
    shell:"""
cat {input.bed} | awk '{{ print $1"."$2"-"$3"/"$4;}}' > {output.asm}
"""


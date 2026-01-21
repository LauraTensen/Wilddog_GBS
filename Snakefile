
################################################################################################
#-------------------------------mapping to reference genome: Snakemake
# DO NOT DOWNLOAD TOOLS OR CONDA

#-------------------------------Edit config
#nano config_I.yaml

#ASSEMBLY: /folder/CamFam4.fa
#OUTDIR: mapping_out
#PATH_WITH_FILES: /folder/gbs_reads

#------------------------------- Make Snakefile
#nano Snakefile_I

#################################################setup#################################################
configfile: "config_I.yaml"
import os
from snakemake.utils import makedirs


#################################################set vars##############################################

# Detect samples as subdirectories

SAMPLE_DIR = config["PATH_WITH_FILES"]
SAMPLES = [d for d in os.listdir(SAMPLE_DIR)
           if os.path.isdir(os.path.join(SAMPLE_DIR, d))]

REF = config["ASSEMBLY"]

if "OUTDIR" in config:
    print("\nSaving to " + config["OUTDIR"] + "\n")
    workdir: config["OUTDIR"]

makedirs("logs_slurm")


####################################################rule all##########################################

rule all:
    input:
        # put all needed output here


####################################################align#############################################

# map to reference
rule bwa_mem2:
    input:
        assembly=REF,
        ref_amb="/folder/CamFam4.fa.amb",
        ref_ann="/folder/CamFam4.fa.ann",
        ref_bwt="/folder/CamFam4.fa.bwt.2bit.64",
        ref_pac="/folder/CamFam4.fa.pac",
        ref_123="/folder/CamFam4.fa.0123",
        reads=lambda wildcards: os.path.join(SAMPLE_DIR, wildcards.sample, f"{wildcards.sample}_R1.fastq.gz"),
        reads2=lambda wildcards: os.path.join(SAMPLE_DIR, wildcards.sample, f"{wildcards.sample}_R2.fastq.gz")
    output:
       temp("mapped_reads/{sample}.bam")
    resources:
        cpus=8,
        mem_mb=32000
    shell:
        """
        rg="{wildcards.sample}"
        module load tools ngs # always put before loading tools
        module load bwa-mem2/2.2.1
        module load samtools/1.20
        echo $rg
        bwa-mem2 mem -t {resources.cpus} -R "@RG\\tID:$rg\\tSM:$rg" {input.assembly} {input.reads} {input.reads2} | samtools view -b - > {output}
        """

# sort aligned reads
rule samtools_sort:
    input:
        rules.bwa_mem2.output
    output:
        "processed_reads/{sample}.sorted.bam"
    message:
        "Rule {rule} processing"
    group:
        'group'
    shell:
        """
        module load tools ngs # always put before loading tools
        module load samtools/1.20
        
        # -m 2G: maximum memorty 2G
        # -@ 7: use 7 threads 

        samtools sort -m 2G -@ 7 -O bam {input} > {output}
        """

# index aligned reads
rule samtools_index:
    input:
        rules.samtools_sort.output
    output:
        "processed_reads/{sample}.sorted.bam.bai"
    message:
        "Rule {rule} processing"
    group:
        "group"
    shell:
        """
        module load tools ngs # always put before loading tools
        module load samtools/1.20

        # -@ 16: use 16 threads 

        samtools index -@ 16 {input}
        """


##############################################quality control#####################################################
# report mapping stats
rule qualimap_report:
    input:
        check=rules.samtools_index.output,
        bam=rules.samtools_sort.output
    output:
        "mapping_stats/qualimap/{sample}/genome_results.txt"
    params:
        outdir = "mapping_stats/qualimap/{sample}/"
    message:
        "Rule {rule} processing"
    group:
        "group"
    shell:
        """
        module load tools ngs # always put before loading tools
        module load gcc/12.2.0 #for qualimap
        module load intel/perflibs/2020_update4 #for qualimap
        module load R/4.4.0 #for qualimap
        module load java/1.8.0 #for qualimap
        module load qualimap/2.3
        unset DISPLAY
        qualimap bamqc -bam {input.bam} --java-mem-size=4G -nt 8 -outformat PDF -outdir {params.outdir}
        """


# summarise all stats in 1 file, fast
rule qualimap_summary:
    input:
        expand("mapping_stats/qualimap/{sample}/genome_results.txt", sample = SAMPLES)
    output:
        "mapping_stats/sample_quality_summary.tsv"
    message:
        'Rule {rule} processing'
    group:
        "group"
    params:
        scripts_dir = os.path.join(workflow.basedir, "scripts/")
    shell:
        """
        sh {params.scripts_dir}create_qualimap_summary.sh
        """


################################################variant calling#########################################

# index ref with samtools faid
rule samtools_faidx:
    input:
        REF
    message:
        "Rule {rule} processing"
    output:
        REF + ".fai"
    group:
        "variants"
    shell:
        """
        module load tools ngs # always put before loading tools
        module load samtools/1.20
        samtools faidx {input}
        """

# call SNPs
rule bcftools_variants:
    input:
        assembly = REF,
        index_assembly = rules.samtools_faidx.output,
        bam_reads=expand("processed_reads/{sample}.sorted.bam", sample=SAMPLES),
        bai_reads=expand("processed_reads/{sample}.sorted.bam.bai", sample=SAMPLES)
    output:
        temp("SNPs/wilddog.vcf.gz")
    message:
        "Rule {rule} processing"
    group:
        "variants"
    shell:
        """
        module load tools ngs # always put before loading tools
        module load samtools/1.20
        module load bcftools/1.9

        # -a: output csv with allelic depth, number of high-quality bases, Phred-scale p-value
        # -q: minimum mapping quality = 20
        # -Q: minimum base quality = 10
        # -I: no indel calling
        # -Ou: output format uncompressed BCF
        # -f: faidx-indexed reference file in FASTA 

        bcftools mpileup -a AD,DP,SP -q 20 -Q 10 -I -Ou -f {input.assembly} {input.bam_reads} | \

        # -f GQ,GP: FORMAT fields output genotype quality & genotype probability
        # -m: model for multiallelic and rare-variant calling
        # -Oz: output compressed VCF
        # -o: specify output file
        
        bcftools call -f GQ,GP -m -Oz -o {output}
        """

# Extract SNPs only
rule bcftools_view:
    input:
        rules.bcftools_variants.output,  
    output:
       "SNPs/wilddog_snps.vcf.gz"  
    message:
        "Rule {rule} processing"
    shell:
        """
        module load tools ngs # always put before loading tools
        module load bcftools/1.9

        # -m2 -M2: only biallelic snps
        # -v: only snps
        # -Oz: output compressed vcf

        bcftools view -m2 -M2 -v snps -Oz -o {output} {input}
        """

# index vcf file
rule bcftools_index:
    input:
        rules.bcftools_view.output
    message:
        "Rule {rule} processing"
    output:
        "SNPs/wilddog.vcf.gz.csi"
    group:
        "variants"
    shell:
        """
        module load tools ngs # always put before loading tools
        module load bcftools/1.9
        bcftools index {input}
        """

# get statistics of variant
rule bcftools_stats:
    input:
        vcf = rules.bcftools_view.output,  
        idx = rules.bcftools_index.output 
    output:
        stats = "SNPs/variant.stats" 
    message:
        "Rule {rule} processing"
    shell:
        """
        module load tools ngs # always put before loading tools
        module load bcftools/1.9
        bcftools stats {input.vcf} > {output.stats}
        """

# stats for filtering
rule vcftools_analysis:
    input:
        vcf = rules.bcftools_view.output,  
        vcf_index = rules.bcftools_index.output 
    output:
        frq = "SNPs/stats/wilddog.frq",
        idepth = "SNPs/stats/wilddog.idepth",
        ldepth_mean = "SNPs/stats/wilddog.ldepth.mean",
        lqual = "SNPs/stats/wilddog.lqual",
        imiss = "SNPs/stats/wilddog.imiss",
        lmiss = "SNPs/stats/wilddog.lmiss"
    shell:
        """
        module load tools ngs
        module load perl/5.36.1
        module load vcftools/0.1.15

        # Allele frequency calculation
        vcftools --gzvcf {input.vcf} --freq2 --out SNPs/stats/wilddog --max-alleles 2

        # Mean depth calculation
        vcftools --gzvcf {input.vcf} --depth --out SNPs/stats/wilddog

        # Mean depth per site calculation
        vcftools --gzvcf {input.vcf} --site-mean-depth --out SNPs/stats/wilddog

        # Site quality calculation
        vcftools --gzvcf {input.vcf} --site-quality --out SNPs/stats/wilddog

        # Missing data per individual calculation
        vcftools --gzvcf {input.vcf} --missing-indv --out SNPs/stats/wilddog

        # Missing data per site calculation
        vcftools --gzvcf {input.vcf} --missing-site --out SNPs/stats/wilddog
        """

# filter out some SNPs
rule vcf_filtering
    input:
        rules.bcftools_view.output
    message:
        "Rule {rule} processing"
    output:
        "SNPs/wilddog_filtered.vcf.gz"
    group:
        "variants"
    shell:
        """
        module load tools ngs # always put before loading tools
        module load perl/5.36.1
        module load vcftools/0.1.15

        # --gzvcf: specify compressed vcf input
        # --maf 0.2: include only sites with minor allele freq <0.2
        # --max-missing: exclude sites with >0.2 missing data

        vcftools --gzvcf {input} --max-missing 0.6 --maf 0.05 --minQ 20 \
        --min-meanDP 2 --max-meanDP 100 --recode --stdout | gzip -c > {output}
        """


################################################data prep#########################################
# convert to plink format: fast
rule plink_vcf:
    input:
        rules.vcf_filtering.output
    output:
        bed = "plink/snps_dataset.bed",
        bim = "plink/snps_dataset.bim",
        fam = "plink/snps_dataset.fam"
    message:
        "Rule {rule} processing"
    group:
        "variants"
    shell:
        """
        module load tools ngs
        module load plink2/1.90beta7

        # --vcf:
        # --make-bed: 
        # --set-missing-var-ids @:#: 
        # --out SNPs: 
        # --allow-extra-chr: 

        plink --vcf {input} --make-bed --set-missing-var-ids @:# \
        --out plink --allow-extra-chr
        """

# remove variants in LD
rule plink_pruning:
    input:
        rules.plink_vcf.output
    output:
        pruned_in = "plink/ld_pruned.prune.in",
        pruned_out = "plink/ld_pruned.prune.out"
    log:
        "logs_slurm/ld_pruning.log"
    message:
        "Rule {rule} processing"
    shell:
        """
        module load tools ngs
        module load plink2/1.90beta7

        # --bfile
        # --indep-pairwise 30 15 0.2
        # --allow-extra-chr

        plink --bfile plink/snps_dataset --indep-pairwise 30 15 0.2 --allow-extra-chr --out plink/ld_pruned
        """

# extract the variants not in LD
rule plink_extract:
        rules.plink_vcf.output,
        pruned_in = "plink/ld_pruned.prune.in"
    output:
        filtered_bed = "plink/snps_dataset_ld_filtered.bed",
        filtered_bim = "plink/snps_dataset_ld_filtered.bim",
        filtered_fam = "plink/snps_dataset_ld_filtered.fam",
        filtered_vcf =  "vcf_pruned/snps_dataset_ld_filtered_noFID.vcf"
    log:
        "logs_slurm/extract_pruned_snps.log"
    message:
        "Rule {rule} processing"
    shell:
        """
        module load tools ngs
        module load plink2/1.90beta7

        plink --bfile plink/snps_dataset --extract {input.pruned_in} --make-bed --allow-extra-chr --out plink/snps_dataset_ld_filtered
        plink --bfile plink/snps_dataset_ld_filtered --recode vcf --allow-extra-chr --out vcf_pruned/snps_dataset_ld_filtered
        awk 'BEGIN {OFS="\t"} /^#CHROM/ {for (i=10; i<=NF; i++) $i = substr($i, index($i, "_")+1); print; next} {print}' vcf_pruned/snps_dataset_ld_filtered.vcf > vcf_pruned/snps_dataset_ld_filtered_noFID.vcf
        """

# get vcf file per population 
rule vcftools_split:
    input:
        "vcf_pruned/snps_dataset_ld_filtered_noFID.vcf"
    output:
        expand("vcf_pruned/snps_{pair}.recode.vcf", pair=["MTP", "KNP", "FRM", "KNP_FRM", "MTP_FRM", "MTP_KNP"])
    message:
        "Rule {rule} processing"
    shell:
        """
        module load tools ngs # always put before loading tools
        module load perl/5.36.1
        module load vcftools/0.1.15

        vcftools --vcf {input} \
        --keep /folder/pop_list/pop_MTP.txt \
        --recode --out vcf_pruned/snps_MTP

        vcftools --vcf {input} \
        --keep /folder/pop_list/pop_KNP.txt \
        --recode --out vcf_pruned/snps_KNP

        vcftools --vcf {input} \
        --keep /folder/pop_list/pop_FRM.txt \
        --recode --out vcf_pruned/snps_FRM

        vcftools --vcf {input} \
        --keep /folder/pop_list/pairpop_KNP_FRM.txt \
        --recode --out vcf_pruned/snps_KNP_FRM

        vcftools --vcf {input} \
        --keep /folder/pop_list/pairpop_MTP_FRM.txt \
        --recode --out vcf_pruned/snps_MTP_FRM

        vcftools --vcf {input} \
        --keep /folder/pop_list/pairpop_MTP_KNP.txt \
        --recode --out vcf_pruned/snps_MTP_KNP
        """

# get raw & map plink files
rule plink_extract:
        "plink/snps_dataset_ld_filtered"
    output:
        filtered_bed = "plink/snps_dataset_ld_filtered.raw",
        filtered_bim = "plink/snps_dataset_ld_filtered.map",
    message:
        "Rule {rule} processing"
    shell:
        """
        module load tools ngs
        module load plink2/1.90beta7

        plink --bfile plink/snps_dataset_ld_filtered --recode --allow-extra-chr --out plink/snps_dataset_ld_filtered
        plink --bfile plink/snps_dataset_ld_filtered --recode A --allow-extra-chr --out plink/snps_dataset_ld_filtered

        """

################################################clean up#########################################

# email when done or error log 

EMAIL = "mail@outlook.com"

onsuccess:
    shell("mail -s 'DONE' {EMAIL} < {log}")

onerror:
   shell("mail -s 'ERROR' {EMAIL} < {log}")

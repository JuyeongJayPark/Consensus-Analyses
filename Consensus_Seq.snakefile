#import os 
#
#current_path = os.getcwd()
#
#samples = []
#with open("WholeGenomeReSequencing.configure.v.0.0.3.txt") as config_file:
#    for line in config_file:
#        if line.startswith("SAMPLE="):
#            sample = line.split("SAMPLE=")[1].split(',')[0].strip()
#            samples.append(sample)
#        elif line.startswith("REFERENCE="):
#            ref_path = line.split("REFERENCE=")[1].strip()
#
#print(f"Sample_list : {samples}")
#print(f"Reference : {ref_path}")

configfile: "config.yaml"

samples = config["samples"].split(" ")
print(f"Samples: {samples}")

rule all:
    input:
        expand("OutPut/09.GenotypeVariant/{sample}/{sample}.vcf.gz", sample = samples),
        expand("OutPut/09.GenotypeVariant/{sample}/{sample}.vcf.gz.tbi", sample = samples),
        expand("OutPut/23.Consensus/{sample}/{sample}.gatk.consensus.fa", sample = samples),
        expand("OutPut/24.Bcftools_VCF/{sample}/{sample}.vcf.gz", sample = samples),
        expand("OutPut/24.Bcftools_VCF/{sample}/{sample}.vcf.gz.tbi", sample = samples),
        expand("OutPut/23.Consensus/{sample}/{sample}.bcftools.consensus.fa", sample = samples),
        expand("OutPut/25.Base_Count/{sample}/{sample}.basecall.xls", sample = samples),
        expand("OutPut/27.FTP/{sample}/02_vcf_file/{sample}.gatk.vcf", sample = samples),
        expand("OutPut/27.FTP/{sample}/02_vcf_file/{sample}.bcftools.vcf", sample = samples),
        expand("OutPut/27.FTP/{sample}/05_consensus_file/{sample}.gatk.consensus.fa", sample = samples),
        expand("OutPut/27.FTP/{sample}/05_consensus_file/{sample}.bcftools.consensus.fa", sample = samples),
        expand("OutPut/27.FTP/{sample}/06_base_count/{sample}.basecall.xls", sample = samples)

rule gatk_compress_vcf:
    input:
        vcf = "OutPut/09.GenotypeVariant/{sample}/{sample}.vcf"
    output: 
        vcf_gz = "OutPut/09.GenotypeVariant/{sample}/{sample}.vcf.gz"
    shell:
        "bgzip -c {input.vcf} > {output.vcf_gz}"

rule gatk_vcf_index:
    input:
        vcf_gz = "OutPut/09.GenotypeVariant/{sample}/{sample}.vcf.gz"
    output: 
        tbi = "OutPut/09.GenotypeVariant/{sample}/{sample}.vcf.gz.tbi"
    shell:
        "tabix -p vcf -f {input.vcf_gz} > {output.tbi}"

rule gatk_consensus_seq:
    input:
        vcf_gz = "OutPut/09.GenotypeVariant/{sample}/{sample}.vcf.gz",
        tbi = "OutPut/09.GenotypeVariant/{sample}/{sample}.vcf.gz.tbi"
    output:
        consensus_fa = "OutPut/23.Consensus/{sample}/{sample}.gatk.consensus.fa"
    params:
        reference = config["ref"]
    shell:
        "cat {params.reference} | bcftools consensus {input.vcf_gz} > {output.consensus_fa}"

rule bcftools_call:
    input:
        bam = "OutPut/05.De-duplication/{sample}/{sample}.dedup.bam"
    output:
        vcf_gz = "OutPut/24.Bcftools_VCF/{sample}/{sample}.vcf.gz"
    params:
        reference = config["ref"]
    shell:
        "bcftools mpileup -Oz -a DP,AD,ADF,ADR -f {params.reference} {input.bam} | bcftools call -mv -Oz -o {output.vcf_gz}"

rule bcftools_vcf_index:
    input:
        vcf_gz = "OutPut/24.Bcftools_VCF/{sample}/{sample}.vcf.gz"
    output: 
        tbi = "OutPut/24.Bcftools_VCF/{sample}/{sample}.vcf.gz.tbi"
    shell:
        "tabix -p vcf -f {input.vcf_gz} > {output.tbi}"

rule bcftools_consensus_seq:
    input:
        vcf_gz = "OutPut/24.Bcftools_VCF/{sample}/{sample}.vcf.gz",
        tbi = "OutPut/24.Bcftools_VCF/{sample}/{sample}.vcf.gz.tbi"
    output:
        consensus_fa = "OutPut/23.Consensus/{sample}/{sample}.bcftools.consensus.fa"
    params:
        reference = config["ref"]
    shell:
        "cat {params.reference} | bcftools consensus {input.vcf_gz} > {output.consensus_fa}"

rule calculate_base_count:
    input:
        bam = "OutPut/05.De-duplication/{sample}/{sample}.dedup.bam"
    output:
        xls = "OutPut/25.Base_Count/{sample}/{sample}.basecall.xls"
    params:
        reference = config["ref"]
    shell:
        "/usr/bin/python2.7 /TBI/Share/NonHumanTeam/Script/WGRS/script/BamToBaseCount.ver1.1.py {params.reference} {input.bam} {output.xls}"

rule upload_ftp:
    input:
        gatk_vcf = "OutPut/09.GenotypeVariant/{sample}/{sample}.vcf",
        bcftools_vcf = "OutPut/24.Bcftools_VCF/{sample}/{sample}.vcf.gz",
        gatk_consensus_seq = "OutPut/23.Consensus/{sample}/{sample}.gatk.consensus.fa",
        bcftools_consensus_seq = "OutPut/23.Consensus/{sample}/{sample}.bcftools.consensus.fa",
        base_count_xls = "OutPut/25.Base_Count/{sample}/{sample}.basecall.xls"
    output:
        ln_gatk_vcf = "OutPut/27.FTP/{sample}/02_vcf_file/{sample}.gatk.vcf",
        ln_bcftools_vcf = "OutPut/27.FTP/{sample}/02_vcf_file/{sample}.bcftools.vcf.gz",
        ln_gatk_consensus_seq = "OutPut/27.FTP/{sample}/05_consensus_file/{sample}.gatk.consensus.fa",
        ln_bcftools_consensus_seq = "OutPut/27.FTP/{sample}/05_consensus_file/{sample}.bcftools.consensus.fa",
        ln_base_count_xls = "OutPut/27.FTP/{sample}/06_base_count/{sample}.basecall.xls"
    shell:
        "ln -s {config[workdir]}/{input.gatk_vcf} {config[workdir]}/{output.ln_gatk_vcf} && "
        "ln -s {config[workdir]}/{input.bcftools_vcf} {config[workdir]}/{output.ln_bcftools_vcf} && "
        "ln -s {config[workdir]}/{input.gatk_consensus_seq} {config[workdir]}/{output.ln_gatk_consensus_seq} && "
        "ln -s {config[workdir]}/{input.bcftools_consensus_seq} {config[workdir]}/{output.ln_bcftools_consensus_seq} && "
        "ln -s {config[workdir]}/{input.base_count_xls} {config[workdir]}/{output.ln_base_count_xls}"


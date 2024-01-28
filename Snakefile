import os
configfile: "config.yaml"

SAMPLES = config['samples']
CHROMOSOMES = [str(i) for i in range(1, 23)] + ['MT','Y','X']
output_dir = config['output_dir']
raw_dir = config['rawdir']
alignment_dir = output_dir + '/alignment'
log_dir = output_dir + '/logs'
benchmark_dir = log_dir + '/benchmarks'
alignment_dir = output_dir + '/alignment'
qc_dir = output_dir + '/qc'
vcf_dir = output_dir + '/vcf'
stats_dir = output_dir + '/stats'
ref_fasta = config['ref']
ref_fasta_without_extension = os.path.splitext(ref_fasta)[0]
Tertiary = config['Tertiary']
num_threads = config['computing_threads']

results_fastqc = expand(qc_dir + '/fastqc/{sample}_sorted_MarkDup_RG_fastqc.html', sample = SAMPLES)
results_qualimap = expand(qc_dir + '/qualimap/{sample}/qualimapReport.html', sample = SAMPLES)
results_VCF = expand(output_dir + '/VCF/{sample}_align_sort_MarkDup_RG.vcf', sample = SAMPLES)
results_report = expand(output_dir + '/Clincal_Report/{sample}_report.html', sample = SAMPLES)


if Tertiary == 1:
    rule all:
        input:
            results_report +
            results_fastqc +
            results_qualimap 
else:
    rule all:
        input:
            results_fastqc +
            results_qualimap 
		
 
rule fastp:
    """
    Trim reads and perform QC using fastp with raw fastq
    """
    input:
        R1 = raw_dir + '/{sample}_R1.fastq.gz',
        R2 = raw_dir + '/{sample}_R2.fastq.gz'
    output:
        R1 = qc_dir + '/fastp/{sample}_R1_fastp.fastq.gz',
        R2 = qc_dir + '/fastp/{sample}_R2_fastp.fastq.gz',
        stats = qc_dir + '/fastp/{sample}_fastp.html'
    log:
        log_dir + '/{sample}.fastp.log'
    threads:
	    num_threads
    benchmark:
        benchmark_dir + '/fastp/{sample}.tsv'
    shell:
        'fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} -h {output.stats} --thread {threads} &> {log}'

rule ref_index:
    """
    Create index of the reference
    """
    input:
        ref_fasta
    output:
        ref_fasta + '.fai'
    benchmark:
        benchmark_dir + '/samtools_index.tsv'
    shell:
        """
        samtools faidx {input} -o {output}
        """

rule mm2:
    """
    mapping and alignment using mm2
    """
    input:
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
        ref = ref_fasta
    output:
        alignment_dir + '/alignment/{sample}.sam'
    log:
        log_dir + '/{sample}.mm2.log'
    threads:
        num_threads
    benchmark:
        benchmark_dir + '/mm2/{sample}.tsv'
    shell:
        """
        minimap2 -ax sr -t {threads} {input.ref} {input.R1} {input.R2} > {output}
        """
		
rule sort:
    """
    Conversion and position sorting
    """
    input:
        rules.mm2.output
    output:
        alignment_dir + '/alignment/{sample}_sorted.bam'
    log:
        log_dir + '/{sample}.samtools_sort.log'
    threads:
	    num_threads
    benchmark:
        benchmark_dir + '/samtools_sort/{sample}.tsv'
    shell:
        """
        samtools sort {input} -o {output} -@ {threads} -m 1G &> {log}
        """

rule markDuplicates:
    """
    Conversion and position sorting
    """
    input:
        rules.sort.output
    output:
        bam = alignment_dir + '/alignment/{sample}_sorted_MarkDup.bam',
        metrics = output_dir + '/stats/{sample}_marked_dup_metrics.txt',
        temp = directory(output_dir + '/temp/{sample}/')
    log:
        log_dir + '/{sample}.markDuplicates.log'
    threads:
	    num_threads
    benchmark:
        benchmark_dir + '/markDuplicates/{sample}.tsv'
    shell:
        """
        gatk MarkDuplicates I={input} O={output.bam} M={output.metrics} TMP_DIR={output.temp} &> {log}
        """

rule fastqc:
    """
    Run QC for each fastq
    """
    input:
        alignment_dir + '/alignment/{sample}_sorted_MarkDup_RG.bam'
    output:
        qc_dir + '/fastqc/{sample}_sorted_MarkDup_RG_fastqc.html'
    log:
        log_dir + '/{sample}_fastqc.log'
    threads:
	    num_threads
    benchmark:
        benchmark_dir + '/fastqc/{sample}.tsv'
    shell:
        'fastqc -o $(dirname {output}) {input} &> {log}'
		

rule qualimap:
    """
    Calculate bam statistics
    """
    input:
        bam = alignment_dir + '/alignment/{sample}_sorted_MarkDup_RG.bam'
    output:
        qc_dir + '/qualimap/{sample}/qualimapReport.html'
    log:
        log_dir + '/{sample}_qualimap.log'
    benchmark:
        benchmark_dir + '/qualimap/{sample}.tsv'
    shell:
        """
        qualimap bamqc -bam {input.bam} -outdir {qc_dir}/qualimap/{wildcards.sample}/ --java-mem-size=30000M &> {log}
        """
		
rule GATK_create_dict:
    """
    Create dictionary of the reference
    """
    input:
        ref_fasta
    output:
        ref_fasta_without_extension + '.dict'
    log:
        log_dir + '/GATK_CreateSequenceDictionary.log'
    benchmark:
        benchmark_dir + '/GATK_CreateSequenceDictionary.tsv'
    shell:
        """
        gatk CreateSequenceDictionary -R {input} -O {output} &> {log}
        """

rule GATK_Read_Group:
    """
    Add Read groups
    """
    input:
        rules.markDuplicates.output.bam
    output:
        alignment_dir + '/alignment/{sample}_sorted_MarkDup_RG.bam'
    log:
        log_dir + '/{sample}_GATK_Read_Group.log'
    benchmark:
        benchmark_dir + '/GATK_Read_Group/{sample}.tsv'
    shell:
        """
        gatk AddOrReplaceReadGroups I={input} O={output} RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20 &> {log}
        """
		
rule bam_index:
    """
    Create Bam index
    """
    input:
       rules.GATK_Read_Group.output
    output:
        alignment_dir + '/alignment/{sample}_sorted_MarkDup_RG.bam.csi'
    log:
        log_dir + '/{sample}_samtools_index.log'
    threads:
        num_threads	
    benchmark:
        benchmark_dir + '/samtools_index/{sample}.tsv'
    shell:
        """
        samtools index {input} -@ {threads} -m 1G &> {log}
        """
		
rule Haplotypecaller:
    """
    Perform Variant calling using Haplotypecaller
    """
    input:
        bam = rules.GATK_Read_Group.output,
        ref = ref_fasta,
        index_ref = rules.ref_index.output,
        dict = rules.GATK_create_dict.output,
        index_bam = rules.bam_index.output
    output:
        vcf = output_dir + '/VCF/{sample}/{sample}_align_sort_MarkDup_RG_{chromosome}.vcf'
    params:
        sample = SAMPLES,
        chromosome = CHROMOSOMES
    threads: 4
    shell:
        """
        gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf} --min-base-quality-score 30 -L {wildcards.chromosome}
        """
        
rule Merge_VCF:
    input:
        vcf_files = lambda wildcards: expand(output_dir + "/VCF/{sample}/{sample}_align_sort_MarkDup_RG_{chromosome}.vcf",
                                             sample=wildcards.sample, 
                                             chromosome=CHROMOSOMES)
    output:
        vcf = output_dir + "/VCF/{sample}/{sample}_merged.vcf"
    shell:
        """
        vcf-concat {input.vcf_files} > {output.vcf}
        """

if Tertiary == 1:
    rule VEP:
        """
        Perform variant annotation using Ensembl VEP
        """
        input:
            rules.Merge_VCF.output
        output:
            output_dir + '/VEP/{sample}_VEP_output.vcf'
        benchmark:
            benchmark_dir + '/VEP/{sample}.tsv'
        log:
            log_dir + '/{sample}_VEP.log'
        shell:
            """
            vep --dir_cache VEP_cache --check_existing --offline -i {input} -o {output} --vcf &> {log}
            """

    rule vcfanno:
        input:
            VCF = rules.VEP.output
        output:
            output_dir + '/Clinvar/{sample}_VEP_clinvar.vcf'
        shell:
            """
            vcfanno clinvar/clinvar.toml {input.VCF} > {output} 
            """
	
    rule filter:
        input:
            rules.vcfanno.output
        output:
            output_dir + '/Clinvar/{sample}_VEP_clinvar_filtered.vcf'
        shell:
            "bcftools view -i 'INFO/CLNDN != \"\" & INFO/CLNDN != \"not_provided\" & INFO/CLNDN != \"not_specified\"' {input} > {output}"
			
    rule report:
        input:
            vcf_file = rules.filter.output
        output:
            report_html = output_dir + '/Clincal_Report/{sample}_report.html'
        shell:
            "python scripts/report.py {input.vcf_file} {output.report_html}"

else:
    print('Tertiary analysis is set to 0 in config file, skipping VEP and clinvar annotation')

if Tertiary == 1:
    onsuccess:
        shell('multiqc {output_dir} -o {output_dir} && rm -rf {output_dir}/temp/')
else:
    onsuccess:
        shell('multiqc {output_dir} -o {output_dir} && rm -rf {output_dir}/temp/')
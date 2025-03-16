from os.path import join

# Configuration
CONFIG = {
    'genome_dir': '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome_star',  # Directory for STAR genome index
    'genome_fasta': '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa',
    'gtf': '/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_H2AK119Ub_cross_V5/COMMON_DATA/gencode.v43.basic.annotation.nochr.gtf',
    'max_threads': 32,  # Maximum threads for heavy processes
    'medium_threads': 16,  # Medium thread count for moderate processes
    'light_threads': 4,  # Light thread count for simple processes
    'fastq_dir': '/beegfs/scratch/ric.broccoli/kubacki.michal/Azenta_projects/90-1127887346/00_fastq'  # Directory with input fastq files
}

# Get sample names from fastq files (1-12)
SAMPLES = [str(i) for i in range(1, 13)]

# Output directories
RESULTS_DIR = 'results'
LOGS_DIR = 'logs'

# Final output files for the complete pipeline
rule all:
    input:
        # FastQC reports
        expand(join(RESULTS_DIR, 'fastqc', '{sample}_R1_001_fastqc.html'), sample=SAMPLES),
        expand(join(RESULTS_DIR, 'fastqc', '{sample}_R2_001_fastqc.html'), sample=SAMPLES),
        # Trimmed FastQC reports
        expand(join(RESULTS_DIR, 'fastqc', '{sample}_R1_trimmed_fastqc.html'), sample=SAMPLES),
        expand(join(RESULTS_DIR, 'fastqc', '{sample}_R2_trimmed_fastqc.html'), sample=SAMPLES),
        # STAR alignments
        expand(join(RESULTS_DIR, 'star', '{sample}', '{sample}_Aligned.sortedByCoord.out.bam'), sample=SAMPLES),
        expand(join(RESULTS_DIR, 'star', '{sample}', '{sample}_Aligned.sortedByCoord.out.bam.bai'), sample=SAMPLES),
        # featureCounts results
        join(RESULTS_DIR, 'counts', 'all_samples_counts.txt'),
        # MultiQC report
        join(RESULTS_DIR, 'multiqc', 'multiqc_report.html')

# FastQC on raw reads
rule fastqc_raw:
    input:
        r1 = join(CONFIG['fastq_dir'], '{sample}_R1_001.fastq.gz'),
        r2 = join(CONFIG['fastq_dir'], '{sample}_R2_001.fastq.gz')
    output:
        html_r1 = join(RESULTS_DIR, 'fastqc', '{sample}_R1_001_fastqc.html'),
        html_r2 = join(RESULTS_DIR, 'fastqc', '{sample}_R2_001_fastqc.html'),
        zip_r1 = join(RESULTS_DIR, 'fastqc', '{sample}_R1_001_fastqc.zip'),
        zip_r2 = join(RESULTS_DIR, 'fastqc', '{sample}_R2_001_fastqc.zip')
    log:
        join(LOGS_DIR, 'fastqc', '{sample}_raw.log')
    threads: CONFIG['light_threads']
    resources:
        mem_mb = 4000
    shell:
        """
        fastqc -t {threads} \
            -o $(dirname {output.html_r1}) \
            {input.r1} {input.r2} 2> {log}
        """

# Trim adapters and low quality bases with Trimmomatic
rule trimmomatic:
    input:
        r1 = join(CONFIG['fastq_dir'], '{sample}_R1_001.fastq.gz'),
        r2 = join(CONFIG['fastq_dir'], '{sample}_R2_001.fastq.gz')
    output:
        r1 = join(RESULTS_DIR, 'trimmed', '{sample}_R1_trimmed.fastq.gz'),
        r2 = join(RESULTS_DIR, 'trimmed', '{sample}_R2_trimmed.fastq.gz'),
        r1_unpaired = join(RESULTS_DIR, 'trimmed', '{sample}_R1_unpaired.fastq.gz'),
        r2_unpaired = join(RESULTS_DIR, 'trimmed', '{sample}_R2_unpaired.fastq.gz')
    log:
        join(LOGS_DIR, 'trimmomatic', '{sample}.log')
    threads: CONFIG['medium_threads']
    resources:
        mem_mb = 8000
    shell:
        """
        trimmomatic PE -threads {threads} \
            {input.r1} {input.r2} \
            {output.r1} {output.r1_unpaired} \
            {output.r2} {output.r2_unpaired} \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
            LEADING:3 TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36 2> {log}
        """

# FastQC on trimmed reads
rule fastqc_trimmed:
    input:
        r1 = rules.trimmomatic.output.r1,
        r2 = rules.trimmomatic.output.r2
    output:
        html_r1 = join(RESULTS_DIR, 'fastqc', '{sample}_R1_trimmed_fastqc.html'),
        html_r2 = join(RESULTS_DIR, 'fastqc', '{sample}_R2_trimmed_fastqc.html'),
        zip_r1 = join(RESULTS_DIR, 'fastqc', '{sample}_R1_trimmed_fastqc.zip'),
        zip_r2 = join(RESULTS_DIR, 'fastqc', '{sample}_R2_trimmed_fastqc.zip')
    log:
        join(LOGS_DIR, 'fastqc', '{sample}_trimmed.log')
    threads: CONFIG['light_threads']
    resources:
        mem_mb = 4000
    shell:
        """
        fastqc -t {threads} \
            -o $(dirname {output.html_r1}) \
            {input.r1} {input.r2} 2> {log}
        """

# Align reads with STAR
rule star_align:
    input:
        r1 = rules.trimmomatic.output.r1,
        r2 = rules.trimmomatic.output.r2,
        index = CONFIG['genome_dir']
    output:
        bam = join(RESULTS_DIR, 'star', '{sample}', '{sample}_Aligned.sortedByCoord.out.bam'),
        log = join(RESULTS_DIR, 'star', '{sample}', '{sample}_Log.final.out')
    params:
        prefix = join(RESULTS_DIR, 'star', '{sample}', '{sample}_')
    log:
        join(LOGS_DIR, 'star', '{sample}.log')
    threads: CONFIG['max_threads']
    resources:
        mem_mb = 40000
    shell:
        """
        STAR --runThreadN {threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes Standard \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 2> {log}
        """

# Index BAM files
rule index_bam:
    input:
        rules.star_align.output.bam
    output:
        join(RESULTS_DIR, 'star', '{sample}', '{sample}_Aligned.sortedByCoord.out.bam.bai')
    log:
        join(LOGS_DIR, 'samtools', '{sample}_index.log')
    threads: CONFIG['light_threads']
    resources:
        mem_mb = 4000
    shell:
        """
        samtools index -@ {threads} {input} 2> {log}
        """

# Count reads with featureCounts
rule feature_counts:
    input:
        bams = expand(join(RESULTS_DIR, 'star', '{sample}', '{sample}_Aligned.sortedByCoord.out.bam'), sample=SAMPLES),
        gtf = CONFIG['gtf']
    output:
        counts = join(RESULTS_DIR, 'counts', 'all_samples_counts.txt')
    log:
        join(LOGS_DIR, 'featureCounts', 'all_samples.log')
    threads: CONFIG['max_threads']
    resources:
        mem_mb = 16000
    shell:
        """
        featureCounts -T {threads} \
            -p -t exon -g gene_id \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bams} 2> {log}
        """

# Generate MultiQC report
rule multiqc:
    input:
        fastqc_raw = expand(join(RESULTS_DIR, 'fastqc', '{sample}_R1_001_fastqc.zip'), sample=SAMPLES),
        fastqc_trimmed = expand(join(RESULTS_DIR, 'fastqc', '{sample}_R1_trimmed_fastqc.zip'), sample=SAMPLES),
        star_logs = expand(join(RESULTS_DIR, 'star', '{sample}', '{sample}_Log.final.out'), sample=SAMPLES),
        counts = join(RESULTS_DIR, 'counts', 'all_samples_counts.txt')
    output:
        report = join(RESULTS_DIR, 'multiqc', 'multiqc_report.html')
    log:
        join(LOGS_DIR, 'multiqc', 'multiqc.log')
    params:
        outdir = join(RESULTS_DIR, 'multiqc')
    shell:
        """
        multiqc -f \
            -o {params.outdir} \
            {RESULTS_DIR}/fastqc \
            {RESULTS_DIR}/star \
            {RESULTS_DIR}/counts 2> {log}
        """

import os, sys
from os.path import basename, join
from glob import glob

WD = config['wd']
FASTQ_DIR = config['orig_fastqs']

DADA2_DIR=config['dada2_sm_dir']
FILTERED_FASTQ_DIR = join(WD, config['filtered_fastqs'])
ERROR_MODEL_DIR = join(WD, config['error_model_dir'])
SEQTABS_DIR = join(WD, config['seqtables_dir'])
MERGED_SEQTAB_DIR = join(WD, config['merged_seqtab_dir'])
TAX_ASSIGNED_DIR = join(WD, config['taxonomy_dir'])
ANALYSES_DIR = join(WD, config['analysis_dir'])

SAMPLES = list(set(glob_wildcards('{fastq_dir}/{{sample}}_{{pair}}.fastq.gz'.format(fastq_dir=FASTQ_DIR)).sample))
r1_filt = expand(os.path.join(FILTERED_FASTQ_DIR, "{sample}_R1_filt.fastq.gz"), sample=SAMPLES)
r2_filt = expand(os.path.join(FILTERED_FASTQ_DIR, "{sample}_R2_filt.fastq.gz"), sample=SAMPLES)


################################################################################
localrules multiqc, dada2_to_phyloseq

rule all:
	input:
		expand(os.path.join(WD,  "0_qc_reports/{sample}_R{read}_fastqc.html"), sample=SAMPLES, read=['1', '2']),
		os.path.join(WD, "0_qc_reports/multiqc_report.html"), 
		expand(os.path.join(ANALYSES_DIR, "dada2_phyloseq_{tax_method}_{min_sample_reads}_ffun_{min_samples}_{min_amplicon_reads}.Rds"), tax_method=config["taxonomy"], min_sample_reads = config["phyloseq"]["min_sample_reads"], min_samples = config["phyloseq"]["min_samples"], min_amplicon_reads = config["phyloseq"]["min_amplicon_reads"])


################################################################################
rule pre_fastqc:
	input: os.path.join(FASTQ_DIR, "{sample}_R{read}.fastq.gz")
	output: os.path.join(WD, "0_qc_reports/{sample}_R{read}_fastqc.html")
	threads: 1
	log: os.path.join(WD, "slurm_logs/preFastqc_{sample}_R{read}")
	shell: """
	   mkdir -p {WD}/0_qc_reports/
	   fastqc {input} --outdir {WD}/0_qc_reports/
	"""
		
################################################################################
rule multiqc:
	input: expand(os.path.join(WD, "0_qc_reports/{sample}_R{read}_fastqc.html"), sample=SAMPLES, read=['1', '2'])
	output: os.path.join(WD, "0_qc_reports/multiqc_report.html")
	threads: 1
	shell: """
	   cd {WD}/0_qc_reports/
	   multiqc {WD}/0_qc_reports/
	"""		

################################################################################
rule filter_and_trim:
	input:
		r1 = os.path.join(FASTQ_DIR, "{sample}_R1.fastq.gz"),
		r2 = os.path.join(FASTQ_DIR, "{sample}_R2.fastq.gz")
	output:
		r1 = os.path.join(FILTERED_FASTQ_DIR, "{sample}_R1_filt.fastq.gz"),
		r2 = os.path.join(FILTERED_FASTQ_DIR, "{sample}_R2_filt.fastq.gz")
	threads: 1
	resources:
		mem=10,
		time=1
	params:
		r1_trunc = config['filter_and_trim']['r1_trunc'],
		r2_trunc = config['filter_and_trim']['r2_trunc']
	log: os.path.join(WD, "slurm_logs/filterTrim_{sample}")
	shell: """
		mkdir -p {FILTERED_FASTQ_DIR}
		Rscript {DADA2_DIR}/scripts/filter_and_trim.R \
			{input.r1} {input.r2} \
			{params.r1_trunc} {params.r2_trunc} \
			{FILTERED_FASTQ_DIR}
	"""

################################################################################
rule learn_errors:
	input: r1_filt, r2_filt
	output:
		r1_error = os.path.join(ERROR_MODEL_DIR, "r1_error_rates.Rds"),
		r2_error = os.path.join(ERROR_MODEL_DIR, "r2_error_rates.Rds")
	threads: 1
	resources:
		mem=30,
		time=24
	log: os.path.join(WD, "slurm_logs/learnErrors.txt")
	shell: """
		mkdir -p {ERROR_MODEL_DIR}
		Rscript {DADA2_DIR}/scripts/error_model.R \
			{ERROR_MODEL_DIR} {input}
	"""

################################################################################
rule sample_inference:
	input:
		error1 = rules.learn_errors.output.r1_error,
		error2 = rules.learn_errors.output.r2_error,
		r1 = os.path.join(FILTERED_FASTQ_DIR, "{sample}_R1_filt.fastq.gz"),
		r2 = os.path.join(FILTERED_FASTQ_DIR, "{sample}_R2_filt.fastq.gz")
	output:
		seq_all  = os.path.join(SEQTABS_DIR, "{sample}_seqCount.Rds"),
		seq_filt = os.path.join(SEQTABS_DIR, "{sample}_seqCount_noChim.Rds")
	threads: 2
	resources:
		mem=20,
		time=2
	log: os.path.join(WD, "slurm_logs/sampleInference_{sample}")
	shell: """
		mkdir -p {SEQTABS_DIR}
		Rscript {DADA2_DIR}/scripts/sample_inference.R \
			{input.error1} {input.error2} {input.r1} {input.r2} {SEQTABS_DIR}
	"""

################################################################################
rule merge_seqtabs:
	input: expand(os.path.join(SEQTABS_DIR, "{sample}_seqCount_noChim.Rds"), sample=SAMPLES)
	output:
		rds = os.path.join(MERGED_SEQTAB_DIR, "dada2_seqtab.Rds"),
		tsv = os.path.join(MERGED_SEQTAB_DIR, "dada2_seqtab.tsv")
	threads: 2
	log: os.path.join(WD, "slurm_logs/mergeSeqtab.txt")
	shell: """
		mkdir -p {MERGED_SEQTAB_DIR}
        	Rscript {DADA2_DIR}/scripts/merge_seqtabs.R {MERGED_SEQTAB_DIR} {input}
	"""

################################################################################
rule assign_taxonomy:
	input: rules.merge_seqtabs.output.rds
	output: os.path.join(TAX_ASSIGNED_DIR, "dada2_{tax_method}_taxonomy.Rds")
	threads: 1
	log: os.path.join(WD, "slurm_logs/assignTax_{tax_method}")
	shell: """
		mkdir -p {TAX_ASSIGNED_DIR}
        Rscript {DADA2_DIR}/scripts/assign_taxonomies.R \
			{input} {wildcards.tax_method} \
			{DADA2_DIR}/taxonomy_databases/ {TAX_ASSIGNED_DIR}
	"""

##############################################################################
rule dada2_to_phyloseq:
	input:
		samp_data = config["phyloseq"]["sample_data_table"],
		amp_table = rules.merge_seqtabs.output.rds,
		tax_table = rules.assign_taxonomy.output
	output:
		#os.path.join(ANALYSES_DIR, "dada2_phyloseq_{tax_method}.Rds"),
		os.path.join(ANALYSES_DIR, "dada2_phyloseq_{tax_method}_{min_sample_reads}_ffun_{min_samples}_{min_amplicon_reads}.Rds")
	threads: 1
	log: os.path.join(WD, "slurm_logs/dada2Phyloseq_{tax_method}")
	shell: """
		mkdir -p {ANALYSES_DIR}
		Rscript {DADA2_DIR}/scripts/dada2phyloseq.R \
			{input.samp_data} {input.amp_table} {input.tax_table} {wildcards.tax_method}\
			{wildcards.min_sample_reads} {wildcards.min_samples} {wildcards.min_amplicon_reads} \
			{ANALYSES_DIR}
	"""

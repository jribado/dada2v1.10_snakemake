# load libraries
library(dada2); packageVersion("dada2")

################################################################################
# specify command line options
args <- commandArgs(TRUE)
fwd_fastq = args[1]
rev_fastq = args[2]
forward_trunc = as.numeric(args[3])
reverse_trunc = as.numeric(args[4])
output_dir = args[5]

################################################################################
# parse files
samp_prefix <- strsplit(basename(fwd_fastq), "_[12]|_R[12]")[[1]][1]
fwd_path <- paste0(output_dir, "/", samp_prefix, "_R1_filt.fastq.gz")
rev_path <- paste0(output_dir, "/", samp_prefix, "_R2_filt.fastq.gz")
print(fwd_path)

print(paste("Filtering fastq files", fwd_fastq, "and", rev_fastq))
filterAndTrim(fwd_fastq, fwd_path, rev_fastq, rev_path,
          truncLen=c(forward_trunc, reverse_trunc),
          maxEE=c(2, 6),
          truncQ=2,
          maxN=0,
          rm.phix=TRUE,
          compress=TRUE,
          verbose=TRUE,
	  multithread=FALSE)

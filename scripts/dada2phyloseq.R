# load libraries
require(data.table); require(dplyr)
require(dada2)
require(phyloseq)
require(genefilter)

# specify input parameters
args = commandArgs(trailingOnly=TRUE)
samp_table <- read.delim(args[1])
amp_table <- readRDS(args[2])
tax_table <- readRDS(args[3])
tax_method <- args[4]
min_reads <- as.numeric(args[5])
min_samps <- as.numeric(args[6])
min_amp_reads <- as.numeric(args[7])
out_path <- args[8]

###############################################################################
# make phyloseq objects
rownames(samp_table) <- samp_table[, 1] ## set rownames
samp_table <- samp_table[, -1]

ps <- phyloseq(sample_data(samp_table),
							 otu_table(amp_table, taxa_are_rows=F),
               tax_table(tax_table))
ps
saveRDS(ps, paste0(out_path, "/dada2_phyloseq_", tax_method, ".Rds"))


###############################################################################
# filter phyloseq object
# remove samples with less than specified number of minimum reads
ps_min_reads <- subset_samples(ps, sample_sums(ps) > min_reads)
# remove noisy taxa
ffun <- filterfun(kOverA(min_samps, min_amp_reads))
ps_filtered  <- filter_taxa(ps_min_reads, ffun, TRUE)
saveRDS(ps_filtered,
	paste0(out_path, "/dada2_phyloseq_", tax_method, "_", min_reads, "_ffun_", min_samps, "_", min_amp_reads, ".Rds"))
	ps_filtered

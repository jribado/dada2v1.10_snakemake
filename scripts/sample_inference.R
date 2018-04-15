# load libraries
library(dada2); packageVersion("dada2")

# specify input parameters
args = commandArgs(trailingOnly=TRUE)
errF <- readRDS(args[1])
errR <- readRDS(args[2])
fastq_fwd  <- args[3]
fastq_rev  <- args[4]
seqtab_out <- args[5]


################################################################################
# run sample inference
print(paste("Processing:", fastq_fwd, "and", fastq_rev))
print("Running dada2 derepFastq() on forward reads...")
derepF <- derepFastq(fastq_fwd, verbose=TRUE)
print("Running dada2 dada() on forward reads...")
ddF <- dada(derepF, err=errF, multithread=TRUE, verbose=T)

print("Running dada2 derepFastq() on reverse reads...")
derepR <- derepFastq(fastq_rev, verbose=TRUE)
print("Running dada2 dada() on reverse reads...")
ddR <- dada(derepR, err=errR, multithread=TRUE, verbose=T)

################################################################################
# merge sample pairs
sample_prefix <- strsplit(basename(fastq_fwd), '[_R]')[[1]][1]
mergers <- vector("list", 1)
names(mergers) <- sample_prefix
merger <- mergePairs(ddF, derepF, ddR, derepR)
mergers[[sample_prefix]] <- merger

################################################################################
# construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
print(paste("Saving results to", paste0(seqtab_out, "/", sample_prefix, "_seqCount.Rds")))
saveRDS(seqtab, paste0(seqtab_out, "/", sample_prefix, "_seqCount.Rds"))
cat("Unique sequences prior to chimera removal:", ncol(seqtab), " for", nrow(seqtab), " samples.")

# check whether sample even has chimeras. Will output the full sequence table in the absence of chimeras. Only an issue for control blank sample/samples with low DNA yield that did not amplify correctly and will be filtered later in the pipeline.
if(ncol(seqtab) > 0){
	seqtab_nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
	saveRDS(seqtab_nochim, paste0(seqtab_out, "/", sample_prefix, "_seqCount_noChim.Rds"))
	message("Unique sequences after to chimera removal:", ncol(seqtab_nochim), " for", nrow(seqtab), " samples.")
} else {
	saveRDS(seqtab, paste0(seqtab_out, "/", sample_prefix, "_seqCount_noChim.Rds"))
	message("No chimeras removed. ")
}

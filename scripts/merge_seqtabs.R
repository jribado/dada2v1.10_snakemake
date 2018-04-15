# load libraries
library(dada2); packageVersion("dada2")
library(readr); packageVersion("readr")

# specify input parameters
args = commandArgs(trailingOnly=TRUE)
seqtable_out <- args[1]
seqtab_paths <- args[2:length(args)]

################################################################################
# load all sequence tables
load_seqtables <- function(seqtab_paths){
  seqtables <- list()
  i = 1
  for (path in seqtab_paths){
    print(paste('Loading seqtable', i , 'of', length(seqtab_paths)))

    seqtab_in <- readRDS(path, refhook = NULL)
    if (class(seqtab_in) != "matrix" | nrow(seqtab_in) == 0 | ncol(seqtab_in) == 0){
      print(paste("Warning, sequence table is incomplete, ignoring:", path))
    }else{
      seqtables[[i]] <- seqtab_in
      i = i+1
    }
  }
  return(seqtables)
}

################################################################################
# merge sequence tables
seqtabs <- load_seqtables(seqtab_paths)
merged_seqtab <- do.call(mergeSequenceTables, seqtabs)
print("Dimensions of merged sequence table:")

################################################################################
# save sequence tables
print(table(nchar(getSequences(merged_seqtab))))
saveRDS(merged_seqtab, paste0(seqtable_out, "/dada2_seqtab.Rds"))
merged_seqtab_tsv <- cbind(data.frame(ID=rownames(data.frame(merged_seqtab))),
                           data.frame(merged_seqtab))
write_tsv(merged_seqtab_tsv, paste0(seqtable_out, "/dada2_seqtab.tsv"))

# load required libraries
library(dada2); packageVersion("dada2")

# specify input parameters
args = commandArgs(trailingOnly=TRUE)
error_output <- args[1]
files        <- args[2:length(args)]
fwd_fastqs   <- files[1:(length(files)/2)]
rev_fastqs   <- files[((length(files)/2)+1):length(files)]

# Learn forward and reverse error rates
print("Learning forward error rates from 1 million reads...")
errF <- learnErrors(fwd_fastqs, nread=1e6, multithread=TRUE)
print("Saving forward error rates...")
saveRDS(errF, paste0(error_output, "/r1_error_rates.Rds"))

print("Learning reverse error rates from 1 million reads...")
errR <- learnErrors(rev_fastqs, nread=1e6, multithread=TRUE)
print("Saving reverse error rates...")
saveRDS(errR, paste0(error_output, "/r2_error_rates.Rds"))

# load libraries
library(dada2); packageVersion("dada2")
library(readr); packageVersion("readr")

# specify input parameters
args = commandArgs(trailingOnly=TRUE)
seqtab     <- readRDS(args[1])
tax_method <- args[2]
tax_ref_path <- args[3]
taxa_out   <- args[4]

# specify training and test sets
if(tax_method=="rdp"){
	train_set   <- paste0(tax_ref_path, "/rdp_train_set_14.fa.gz")
	species_set <- paste0(tax_ref_path, "/rdp_species_assignment_14.fa.gz")
} else if(tax_method=="silva"){
	train_set   <- paste0(tax_ref_path, "/silva_nr_v123_train_set.fa.gz")
	species_set <- paste0(tax_ref_path, "/silva_species_assignment_v123.fa.gz")
} else{
	print("Invald option: choose 'rdp' or 'silva'. Other taxonomic assignment methods are unavaliable.")
}

# run taxonomic classification
message("Assigning initial taxonomies using ", tax_method, " training set...")
taxa <- assignTaxonomy(seqtab, train_set, verbose=TRUE, multithread=TRUE)
message("Assigning species using ", tax_method, " training set...")
taxa <- addSpecies(taxa, species_set, allowMultiple=TRUE, verbose=TRUE)
taxa_df <- data.frame(taxa)
taxa_df$Sequence <- rownames(taxa)

# save data
message("Saving data from ", tax_method, " taxonomy assignments...")
saveRDS(taxa, paste0(taxa_out, "/dada2_", tax_method, "_taxonomy.Rds"))
write_tsv(taxa_df, paste0(taxa_out, "/dada2_", tax_method, "_taxonomy.tsv"))

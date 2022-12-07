library(plyr)
library(tidyr)
library(readr)
library(jsonlite)
library(tidyverse)
suppressMessages(library(DECIPHER))
suppressMessages(library(Biostrings))
library(Rsamtools)
library(GenomicAlignments)
library(qgraph)
library(wordspace)
library(data.table)
library(argparser)

parser <- arg_parser("Commandline parser")
parser <- add_argument(
  parser,
  "--run",
  help = "Run"
)
parser <- add_argument(
  parser,
  "--bam_file",
  help = "File containing the NGS reference data"
)
parser <- add_argument(
  parser,
  "--umi_cutoff_R9",
  help = "cutoff value for filtering reads"
)
parser <- add_argument(
  parser,
  "--umi_cutoff_V14",
  help = "cutoff value for filtering reads"
)


argv <- parse_args(parser)
run <- argv$run
nanostat_summary <- argv$nanostat_summary
mutserve_summary <- argv$mutserve_summary
ngs_data <- argv$ngs_data
corresponding_positions <- argv$corresponding_positions
umi_cutoff <- ifelse(
  str_detect(run, "V14"),
  argv$umi_cutoff_V14,
  argv$umi_cutoff_R9
)

clip_sequences <- function(i) {
  soft_clip <- unlist(explodeCigarOpLengths(bam_file[[1]]$cigar[i], ops="S"))
  return(subseq(bam_file[[1]]$seq[i], as.numeric(soft_clip[1]), -as.numeric(soft_clip[2])))
}



input_dir <- "cladogram"
run <- "run12_V14"
sample <- "barcode20"
# run <- "run12_V14"
# sample <- "SAPHIR_4624_5104"
# run <- "run10_V14"
# sample <- "A_B_950_50_5104"
STR_start <- 2472
STR_end <- 2506 ### adapted to 2506 instead of 2505
STR_range_start <- 2450
STR_range_end <- 2570

mutation_classification <- read_csv(mutation_classification_path)

bam_file <- scanBam(bam_file_path)


sequences <- DNAStringSet(bam_file[[1]]$seq)
names(sequences) <- bam_file[[1]]$qname

writeXStringSet(sequences, paste(input_dir, sample, "sequences.fasta", sep = "/"))
sequences_clipped <- unlist(DNAStringSetList(sapply(seq_along(sequences), clip_sequences)))
names(sequences_clipped) <- names(sequences)
writeXStringSet(sequences_clipped, paste(input_dir, sample, "clipped.fasta", sep = "/"))
sequences_clipped_unique <- DNAStringSet(unique(sequences_clipped))
# names(sequences_clipped_unique) <- seq_along(sequences_clipped_unique)
# writeXStringSet(sequences_clipped_unique, paste(input_dir, sample, "clipped_unique.fasta", sep = "/"))

sequences_clipped_aligned <- 
  readDNAStringSet(
    filepath = paste(input_dir, sample, "clipped_multiple_alignment.fasta", sep = "/"), 
    format = "fasta"
    )
sequences_clipped <- unlist(DNAStringSetList(sapply(seq_along(sequences_clipped_aligned), clip_sequences)))
names(sequences_clipped) <- names(sequences_clipped_aligned)

sequences_clipped_STR_removed <- remove_STR_region(sequences_clipped_aligned)

sequence_table <- as.data.frame(sequences_clipped_STR_removed) %>% 
  rename(sequence = x)
n_positions <- max(str_count(sequence_table$sequence)) + 1
clusters = row.names(sequence_table)
n_clusters <- nrow(sequence_table)
variant_cutoff <- 0.005
lower_limit <- ceiling(n_clusters * variant_cutoff)
upper_limit <- floor(n_clusters - lower_limit)

sequence_table_transposed <- sequence_table %>%
  mutate( cluster = clusters) %>% 
  separate(sequence, sep = "", into = as.character(seq_along(1:n_positions))) %>% 
  pivot_longer(cols = as.character(seq_along(1:n_positions)), names_to = "position", values_to = "nucleotide")

n_distinct_nucleotides_per_pos <- sequence_table_transposed %>% 
  group_by(position, nucleotide) %>% 
  summarise(n_nucleotides = n()) %>% 
  filter(n_nucleotides < upper_limit & n_nucleotides > lower_limit) %>%
  group_by(position) %>% 
  summarise(n_variants = n_distinct(nucleotide)) %>% 
  filter(n_variants > 1)

sequence_table_haplotypes <- sequence_table_transposed %>% 
  inner_join(n_distinct_nucleotides_per_pos, by = "position") %>% 
  group_by(cluster) %>% 
  summarise(haplotype = paste(nucleotide, collapse = ""),
            positions = paste(position, collapse = ","))
positions_string <- unique(sequence_table_haplotypes$positions)

sequence_table_haplotypes_unique <- sequence_table_haplotypes %>% 
  group_by(haplotype) %>% 
  summarize(haplotype_occurences = n(), 
            cluster = paste(cluster, collapse = ", "), 
            positions = positions_string) %>% 
  filter(haplotype_occurences > 4)


sequence_table_haplotypes_unique_adist <- sequence_table_haplotypes_unique %>% 
  mutate(dist = adist(haplotype))

writeXStringSet(toDNAStringSet(sequence_table_haplotypes_unique), paste(input_dir, sample, "haplotype.fasta", sep = "/"))
write_tsv(sequence_table_haplotypes_unique, paste(input_dir, sample, "haplotype.tsv", sep = "/"))

haplotypes_TypeA <- sequence_table_haplotypes_unique %>% 
  filter(grepl("CTCCCAGGAAACGC", haplotype)) 

n_positions_typeA <- max(str_count(haplotypes_TypeA$haplotype)) + 1
clusters_typeA <- row.names(haplotypes_TypeA)
n_clusters_typeA <- nrow(haplotypes_TypeA)
lower_limit_typeA <- ceiling(n_clusters_typeA * variant_cutoff)
upper_limit_typeA <- floor(n_clusters_typeA - lower_limit_typeA)

sequence_table_transposed_typeA <- haplotypes_TypeA %>%
  separate(haplotype, sep = "", into = as.character(seq_along(1:n_positions_typeA))) %>% 
  pivot_longer(cols = as.character(seq_along(1:n_positions_typeA)), names_to = "position", values_to = "nucleotide")

n_distinct_nucleotides_per_pos_typeA <- sequence_table_transposed_typeA %>% 
  group_by(position) %>% 
  summarise(n_variants = n_distinct(nucleotide)) %>% 
  filter(n_variants > 1) %>% 
  select(position)

sequence_table_haplotypes_typeA <- sequence_table_transposed_typeA %>% 
  inner_join(n_distinct_nucleotides_per_pos_typeA, by = "position") %>% 
  group_by(cluster) %>% 
  summarise(haplotype = paste(nucleotide, collapse = ""),
            positions = paste(position, collapse = ","))

sequence_table_haplotypes_unique_typeA <- sequence_table_haplotypes_typeA %>% 
  distinct(haplotype, .keep_all = TRUE)

# haplotypes_TypeA_edist <- haplotypes_TypeA %>%
#   transmute(adist(haplotype)) %>% 
#   rename(haplotypes_TypeA$cluster)
# names(haplotypes_TypeA_edist) <- haplotypes_TypeA$cluster


writeXStringSet(toDNAStringSet(sequence_table_haplotypes_unique_typeA), paste(input_dir, sample, "haplotype_TypeA.fasta", sep = "/"))
write_tsv(sequence_table_haplotypes_unique_typeA, paste(input_dir, sample, "haplotypes_typeA.tsv", sep = "/"))



# fas <- readDNAStringSet("cladogram/barcode22/clipped_unique_multiple_alignment.fasta")
# masked_DNA <- MaskAlignment(fas, showPlot = TRUE)
# masked_DNA_ranges <- MaskAlignment(fas, type = "ranges")
# DNA_replace_masked_regions <- replaceAt(fas, masked_DNA_ranges)
# sequences_clipped_unique <- DNAStringSet(unique(DNA_replace_masked_regions))
# writeXStringSet(sequences_clipped_unique, paste(input_dir, sample, "removed_masked_unique.fasta", sep = "/"))
# values <- MaskAlignment(fas, type="values")
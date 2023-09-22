library(stringr)
library(argparser)
library(tidyr)
library(purrr)
library(readr)
library(jsonlite)
library(dplyr)
library(Biostrings)
parser <- arg_parser("Commandline parser")
parser <- add_argument(
  parser,
  "--sample_sheet",
  help = "Sample sheet to join barcodes and samples"
)
parser <- add_argument(
  parser,
  "--run",
  help = "Run name"
)
parser <- add_argument(
  parser,
  "--barcode",
  help = "Barcode name"
)
parser <- add_argument(
  parser,
  "--aligned_fasta",
  help = "MAFFT aligned fasta file"
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
sample_sheet <- argv$sample_sheet
barcode <- argv$barcode
run <- argv$run
aligned_fasta <- argv$aligned_fasta
variant_cutoff <- ifelse(
  str_detect(run, "V14"),
  as.numeric(argv$umi_cutoff_V14),
  as.numeric(argv$umi_cutoff_R9)
)

# sample_sheet <- "~/UMI_LPA_KIV2/run12_V14/lib/Barcode_Sample_overview.js"
# barcode <- "barcode20"
# run <- "run12_V14"
# aligned_fasta <- "~/UMI_LPA_KIV2/cladogram/barcode20/clipped_multiple_alignment.fasta"
# variant_cutoff <- 0.005

find_STR_start <- function(sequence){
  matches <- matchPattern("-", sequence)
  STR_start_matches <- start(matches)[between(start(matches), STR_range_start, STR_range_end)]
  if(!is_empty(STR_start_matches)) {
    return(min(STR_start_matches))
  }
  return(0)
}
find_STR_end <- function(sequence){
  matches <- matchPattern("-", sequence)
  STR_end_matches <- end(matches)[between(end(matches), STR_range_start, STR_range_end)]
  if(!is_empty(STR_end_matches)) {
    return(max(STR_end_matches))
  }
  return(0)
}

remove_STR_region <- function(sequences){
  STR_ranges <- list()
  STR_ranges[["start"]] <- sapply(sequences, find_STR_start)
  STR_ranges[["end"]] <- sapply(sequences, find_STR_end)
  min_STR_start <- min(STR_ranges$start)
  max_STR_end <- max(STR_ranges$end)

  if( min_STR_start > 0 & max_STR_end > 0){
    STR_IRange <- IRanges(start = min_STR_start, end = max_STR_end)
    sequences_STR_removed <- replaceAt(sequences, STR_IRange)
    return (sequences_STR_removed)  
  }
  return(sequences)
  
}
toDNAStringSet <- function(haplotype_table){
  haplotypes <- DNAStringSet(as.vector(haplotype_table$haplotype))
  names(haplotypes) <- haplotype_table$cluster
  return( haplotypes )
}
get_sample_name <- function(barcode, sample_sheet){
  barcode_nanopore <- paste0("NB", str_sub(barcode, start = -2))
  sample_barcode_overview <- fromJSON(sample_sheet)
  sample_frame <- sample_barcode_overview %>% filter(Barcode == barcode_nanopore)
  return(as.character(sample_frame["Sample"]))
}

STR_range_start <- 2450
STR_range_end <- 2570
sample_name <- get_sample_name(barcode, sample_sheet)
# sample <- str_sub(sample_name, end = -6)
fragment <- as.numeric(str_sub(sample_name, start = -4))

sequences_clipped <- 
  readDNAStringSet(
    filepath = aligned_fasta, 
    format = "fasta"
  )

if(fragment == 5104){
  sequences_clipped <- remove_STR_region(sequences_clipped)
}
sequence_table <- as.data.frame(sequences_clipped)
names(sequence_table) <- "sequence"
n_positions <- max(str_count(sequence_table$sequence)) + 1
clusters = row.names(sequence_table)
n_clusters <- nrow(sequence_table)
lower_limit <- floor(n_clusters * variant_cutoff)
upper_limit <- ceiling(n_clusters - lower_limit)

sequence_table_transposed <- sequence_table %>%
  mutate( cluster = clusters) %>% 
  separate(sequence, sep = "", into = as.character(seq_along(1:n_positions))) %>% 
  pivot_longer(cols = as.character(seq_along(1:n_positions)), names_to = "position", values_to = "nucleotide")

### filter noise with variant_cutoff dependant upper_ and lower_limit
### summarize data per position and keep only positions containing more than one variant 
n_distinct_nucleotides_per_pos <- sequence_table_transposed %>% 
  group_by(position, nucleotide) %>% 
  summarise(n_nucleotides = n()) %>% 
  filter(n_nucleotides <= upper_limit & n_nucleotides >= lower_limit) %>%
  group_by(position) %>% 
  summarise(n_variants = n_distinct(nucleotide)) %>% 
  filter(n_variants > 1)

if(nrow(n_distinct_nucleotides_per_pos) > 0){
  sequence_table_transposed <- sequence_table_transposed %>% 
    inner_join(n_distinct_nucleotides_per_pos, by = "position")
}
sequence_table_haplotypes <- sequence_table_transposed %>% 
  group_by(cluster) %>% 
  summarise(haplotype = paste(nucleotide, collapse = ""),
            positions = paste(position, collapse = ","))
positions_string <- unique(sequence_table_haplotypes$positions)

sequence_table_haplotypes_unique <- sequence_table_haplotypes %>% 
  group_by(haplotype) %>% 
  summarize(haplotype_occurences = n(), 
            cluster = paste(cluster, collapse = ", "), 
            positions = positions_string) %>% 
  filter(haplotype_occurences > lower_limit)

writeXStringSet(toDNAStringSet(sequence_table_haplotypes_unique), paste(sample_name, "haplotype.fasta", sep = "_"))
write_tsv(sequence_table_haplotypes_unique, paste(sample_name, "haplotype.tsv", sep = "_"))
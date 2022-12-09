library(Rsamtools)
library(GenomicAlignments)
library(argparser)

parser <- arg_parser("Commandline parser")
parser <- add_argument(
  parser,
  "--bam_file",
  help = "File containing the NGS reference data"
)

argv <- parse_args(parser)
bam_file_path <- argv$bam_file

# bam_file_path <- "~/UMI_LPA_KIV2/run12_V14/ont_pl/barcode20/align/final/final.bam"

clip_sequences <- function(i) {
  soft_clip <- unlist(explodeCigarOpLengths(bam_file[[1]]$cigar[i], ops="S"))
  sequence <- bam_file[[1]]$seq[i]
  return(subseq(sequence, as.numeric(soft_clip[1]), -as.numeric(soft_clip[2])))
}

remove_short_sequences <- function(sequences){
  median <- median(width(sequences))
  shortest_sequence <- floor(median / 2)
  return( sequences[width(sequences) > shortest_sequence] )
}

bam_file <- scanBam(bam_file_path)
sequences <- DNAStringSet(bam_file[[1]]$seq)
names(sequences) <- bam_file[[1]]$qname
sequences_filtered <- remove_short_sequences(sequences)

writeXStringSet(sequences_filtered, "sequences.fasta")
sequences_clipped <- unlist(DNAStringSetList(sapply(seq_along(sequences_filtered), clip_sequences)))
names(sequences_clipped) <- names(sequences_filtered)
writeXStringSet(sequences_clipped, "clipped.fasta")

suppressMessages(library(DECIPHER))
suppressMessages(library(Biostrings))
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

clip_sequences <- function(i) {
  soft_clip <- unlist(explodeCigarOpLengths(bam_file[[1]]$cigar[i], ops="S"))
  return(subseq(bam_file[[1]]$seq[i], as.numeric(soft_clip[1]), -as.numeric(soft_clip[2])))
}

bam_file <- scanBam(bam_file_path)
sequences <- DNAStringSet(bam_file[[1]]$seq)
names(sequences) <- bam_file[[1]]$qname

writeXStringSet(sequences, "sequences.fasta")
sequences_clipped <- unlist(DNAStringSetList(sapply(seq_along(sequences), clip_sequences)))
names(sequences_clipped) <- names(sequences)
writeXStringSet(sequences_clipped, "clipped.fasta")
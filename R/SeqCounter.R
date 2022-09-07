# SeqCounter
#' Counting sequence occurrences in arbitrary subsets of the genome
#'
#' This function, given search sequences, reference genome and genomic regions,
#' returns the overall count of the number of exact matches for each sequence in
#' the genomic regions.
#' 
#' @usage SeqCounter(genome,regions,patterns)
#' @param genome The reference genome, BSgenome object
#' @param regions The genomic regions, GRanges object
#' @param patterns The search sequences, XStringSet object
#' @return A named vector with the overall count for each pattern
#' @author Luca Procaccio\cr Politecnico di Milano\cr Mantainer: Luca Procaccio\cr
#' E-mail: <luca.procaccio@@mail.polimi.it>
#' 
#' @examples 
#' 
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(Biostrings)
#' 
#' hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
#' promoters <- promoters(txdb, upstream=2000, downstream=200)
#' promoters <- trim(promoters)
#' patterns <- DNAStringSet(c("TATAAAA", "TATATAT", "TATAAAT"))
#' 
#' SeqCounter(hg38,promoters,patterns)
#' 
#' @importFrom GenomicRanges GRanges
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings vcountPDict
#' @importFrom methods is
#' 
#' @export

SeqCounter <- function(genome, regions, patterns) {
  
  # Check missing arguments
  if (missing(genome)) { stop("Genome not found") }
  if (missing(regions)) { stop("Regions not found") }
  if (missing(patterns)) { stop("Patterns not found") }
  
  # Check correctness of class
  if (!is(genome, "BSgenome")) { 
    stop("'genome' must be a BSgenome object") }
  if (!is(regions, "GRanges")) { 
    stop("'regions' must be a GRanges object") }
  if (!is(patterns,"DNAStringSet")) { 
    stop("'patterns' must be a DNAStringSet object") }
  
  # Get the corresponding sequences from the reference
  region_seqs <- getSeq(genome,regions)
  
  # Count occurences between patterns and region sequences
  occurences <- vcountPDict(patterns, region_seqs)
  
  # Get sum by rows
  summed_occurences <- apply(occurences, 1, sum)
  names(summed_occurences) <- patterns
  
  # Return a named vector with the overall count for each pattern
  return(summed_occurences)
  
}
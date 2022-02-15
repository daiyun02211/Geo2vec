#' @title Sample motifs on a subset of genome.
#'
#' @description A function used to extract and sample given motifs from a user defined sub-regions of the genome (Inherited from github.com/ZW-xjtlu/RegionPropertiesFeatures)
#'
#' @param motif A character string indicating the query sequence or motifs; Vague mapping rules in \code{\link{IUPAC_CODE_MAP}} is supported when \code{fixed} = FALSE.
#'
#' @param region The \code{\link{GRanges}} object defines the subset region of the genome.
#'
#' @param sequence The \code{\link{BSgenome}} object containing the sequence of the genome.
#'
#' @param fixed FALSE to support the vague mapping rule of the character string, default is FALSE.
#'
#' @param N Number of ranges sampled, by default, it returns all the matched ranges without sub-sampling.
#'
#' @param replace Whether sample with replacement or not, default is FALSE.
#'
#' @return A \code{GRanges} object contains the (sampled) mapped regions of the query sequence on the given subset of the genome.
#'
#' @examples
#' 
#' ## Sample DRACH motif on a given range:
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' exons_gr <- GRanges(seqnames = c('chr10', 'chr19'),
#'                     IRanges(start = c(73813868, 3976270),
#'                             end = c(73813908, 3976310)),
#'                     strand = c('-', '-'))
#' motif_gr <- sampleSequence('DRACH', exons_gr, BSgenome.Hsapiens.UCSC.hg38, fixed = F)
#' 
#' 
#' @importFrom GenomicRanges reduce GRanges
#' @importFrom GenomicFeatures mapFromTranscripts
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Biostrings vmatchPattern
#' @export
sampleSequence <- function(motif, region, bsgenome, N = NULL, fixed = FALSE, replace = FALSE){
  stopifnot(is(region, "GRangesList")|is(region, "GRanges"))
  if(is(region, "GRangesList")) region <- unlist(region)
  region <- reduce(region)
  
  region_dnass <- getSeq(x=bsgenome, 
                         names=seqnames(region), 
                         start=start(region), 
                         end=end(region),
                         strand=strand(region), 
                         as.character=FALSE)
  indx <- paste0("reg_", seq_along(region))
  regions_GRL <- split(region, indx)
  regions_GRL <- regions_GRL[indx]
  rm(indx)
  vmp <- vmatchPattern(motif, region_dnass, fixed = fixed)
  rm(region_dnass)
  vmp_gr <- GRanges(seqnames = rep(names(regions_GRL), elementNROWS(vmp)), ranges = unlist(vmp))
  rm(vmp)
  motif_on_regions <- mapFromTranscripts(vmp_gr, regions_GRL)
  rm(vmp_gr, regions_GRL)
  mcols(motif_on_regions) = NULL
  if(is.null(N)) N <- length(motif_on_regions)
  if(replace == FALSE) N <- min(N, length(motif_on_regions))
  indx2 <- sample.int(length(motif_on_regions), N, replace = replace)
  motif_on_regions <- motif_on_regions[indx2]
  rm(indx2)
  seqlengths(motif_on_regions) <- seqlengths(region)
  return(motif_on_regions)
}
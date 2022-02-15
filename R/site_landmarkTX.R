#' Generate landmarkTX for single-nucleotide level data.
#' 
#' @description This function accepts GRanges data with width 1 and generate landmarkTX for each input.
#' 
#' @param x A \code{GRanges} object 
#' @param txdb A transcript database, currently only TxDb and EnsDb are supported.
#' @param long_tx A \code{logical} object indicating whether to consider only the longest transcript for each site.
#' @param mRNA A \code{character} object indicating whether to consider only protein coding transcripts.
#' @return A \code{\link{GRanges}} object whose metacolumns are the generated landmarkTX.
#' @examples 
#' 
#' ## Generate landmarkTX for example data.
#' library(EnsDb.Hsapiens.v86)
#' input <- import.bed(system.file("extdata", "example.bed", package = "Geo2vec"))
#' encoding <- site_landmarkTX(input, EnsDb.Hsapiens.v86, long_tx=T, mRNA=T)
#' 
#' @import magrittr
#' @import GenomicFeatures
#' @export
site_landmarkTX <- function(x, txdb, long_tx=TRUE, mRNA=TRUE){
  db_type <- class(txdb)[1]
  if (!(db_type %in% c('TxDb', 'EnsDb'))){
    stop("Error: currently only TxDb and EnsDb are supported.")
  } 
  mcols(x) <- NULL
  if (db_type == 'EnsDb'){
    seqlevelsStyle(x) <- 'Ensembl'
    if (mRNA){
      filter <- AnnotationFilterList(SeqNameFilter(seqlevels(x)),
                                     TxBiotypeFilter('protein_coding'))
      txdb <- addFilter(txdb, filter)
    }
  }
  if (db_type == 'EnsDb'){
    tx <- transcripts(txdb)
    exbytx <- exonsBy(txdb, by='tx')
  }else{
    seqlevelsStyle(x) <- 'UCSC'
    tx <- transcripts(txdb)
    exbytx <- exonsBy(txdb, by='tx', use.names=T)
    if (mRNA){
      cdsbytx <- cdsBy(txdb, by='tx', use.names=TRUE)
      coding_names <- names(cdsbytx) %>% unique()
      tx <- tx[coding_names]
      exbytx <- exbytx[coding_names]
    }
  } 
  
  exbytx <- unlist(exbytx)
  # Currently landmark only applies to single base sites
  if (mean(overlapsAny(x, exbytx)) != 1){
    stop("Error: Some targets can not be mapped to exons! Currently landmarkTX only supports sites from exon.")
  }
  
  xbyex <- findOverlaps(x, exbytx)
  if (long_tx){
    xbyex <- mapLongTX(xbyex, exbytx, tx, exon_only=TRUE)
  }
  
  target_name <- names(exbytx[xbyex %>% subjectHits()])
  map_name <- paste0(xbyex %>% queryHits(), '-', target_name)
  
  x <- x[queryHits(xbyex)]
  names(x) <- map_name
  
  # Dist to site-containing transcript
  tx <- transcripts(txdb)
  names(tx) <- tx$tx_name
  map2tx <- pmapToTranscripts(x, tx[target_name])
  x$tx_u5 <- start(map2tx) - 1
  x$tx_u3 <- width(tx[target_name]) - end(map2tx)
  
  # Dist to site-containing exon
  map2ex <- pmapToTranscripts(x, exbytx[subjectHits(xbyex)])
  x$ex_u5 <- start(map2ex) - 1
  x$ex_u3 <- width(exbytx[subjectHits(xbyex)]) - end(map2ex)
  
  # Although an exon annotation may contain both CDS and UTR
  # the distance to UTR will provide CDS/UTR boundary info.
  
  # Dist to the border between CDS and UTRs
  if (db_type == 'EnsDb'){
    cdsbytx <- cdsBy(txdb, by='tx')
  }else{
    cdsbytx <- cdsBy(txdb, by='tx', use.names=T)
  }
  # range(cdsbytx) is slower
  cds_min <- lapply(start(cdsbytx), min) %>% unlist()
  cds_max <- lapply(end(cdsbytx), max) %>% unlist()
  cds_ranges <- GRanges(seqnames(tx[cdsbytx %>% names()]),
                        IRanges(cds_min, cds_max),
                        strand(tx[cdsbytx %>% names()]))
  cds_starts <- resize(cds_ranges, width=1, fix='start')
  cds_ends <- resize(cds_ranges, width=1, fix='end')
  
  mcols(tx) <- NULL
  tx_starts <- resize(tx, width=1, fix='start')
  tx_ends <- resize(tx, width=1, fix='end')
  tx_starts[names(cdsbytx)] <- cds_starts
  tx_ends[names(cdsbytx)] <- cds_ends
  
  strand_effect <- rep(1, length(x))

  strand_effect[which(as.vector(strand(x)) == '-')] <- -1
  
  x$cds_utr5 <- (start(resize(x, width=1, fix='start')) -
                   start(tx_starts[target_name])) * strand_effect
  
  x$cds_utr3 <- (start(tx_ends[target_name]) -
                   start(resize(x, width=1, fix='end'))) * strand_effect
  
  return(x)
}
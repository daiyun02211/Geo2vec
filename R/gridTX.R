#' Generate gridTX for single-nucleotide level data.
#' 
#' @description This function accepts GRanges data with width 1 and generate gridTX for each input.
#' 
#' @param x A \code{GRanges} object 
#' @param txdb A transcript database, currently only TxDb and EnsDb are supported.
#' @param ngrid An integer indicating the number of grids, i.e., how many pieces will the transcript be divided into.
#' @param exon_only A \code{logical} object indicating whether to consider only sites from exonic regions.
#' @param long_tx A \code{logical} object indicating whether to consider only the longest transcript for each site.
#' @param mRNA A \code{character} object indicating whether to consider only protein coding transcripts.
#' @param unify_strand A \code{logical} object indicating whether to unify the direction from 5'end to 3'end.
#' @return A \code{\link{GRanges}} object whose metacolumns are the generated gridTX.
#' @examples 
#' 
#' ## Generate gridTX for example data.
#' library(EnsDb.Hsapiens.v86)
#' input <- import.bed(system.file("extdata", "example.bed", package = "Geo2vec"))
#' encoding <- gridTX(input+10, EnsDb.Hsapiens.v86, ngrid=40, exon_only=T, long_tx=T, mRNA=T)
#' 
#' @import magrittr
#' @export
gridTX <- function(x, txdb, ngrid=40, exon_only=TRUE, long_tx=TRUE, mRNA=TRUE, unify_strand=TRUE){
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
    cdsbytx <- cdsBy(txdb, by='tx')
    inbytx <- intronsByTranscript(txdb)
    utr5bytx <- fiveUTRsByTranscript(txdb)
    utr3bytx <- threeUTRsByTranscript(txdb)
    coding_names <- c(names(cdsbytx), names(utr5bytx), names(utr3bytx))
    coding_names <- coding_names %>% unique()
  }else{
    seqlevelsStyle(x) <- 'UCSC'
    tx <- transcripts(txdb)
    names(tx) <- tx$tx_name
    exbytx <- exonsBy(txdb, by='tx', use.names=TRUE)
    cdsbytx <- cdsBy(txdb, by='tx', use.names=TRUE)
    inbytx <- intronsByTranscript(txdb, use.names=TRUE)
    utr5bytx <- fiveUTRsByTranscript(txdb, use.names=TRUE)
    utr3bytx <- threeUTRsByTranscript(txdb, use.names=TRUE)
    coding_names <- c(names(cdsbytx), names(utr5bytx), names(utr3bytx))
    coding_names <- coding_names %>% unique()
    if (mRNA){
      tx <- tx[coding_names]
      exbytx <- exbytx[coding_names]
    }
  }
  
  x_start <- GRanges(seqnames=seqnames(x),
                     IRanges(start=start(x)),
                     strand=strand(x))
  
  if (!exon_only){
    xbytx <- findOverlaps(x_start, tx)
    xbytx <- xbytx[end(x[xbytx %>% queryHits()]) <=
                     end(tx[names(tx[xbytx %>% subjectHits()])])]    
  } else {
    xbytx <- findOverlaps(x_start, exbytx)
    x_end <- GRanges(seqnames=seqnames(x),
                     IRanges(start=end(x)),
                     strand=strand(x))
    xebyex <- findOverlaps(x_end, exbytx)
    
    sname <- paste0(xbytx %>% queryHits(), '-', 
                    names(exbytx[xbytx %>% subjectHits()]))
    ename <- paste0(xebyex %>% queryHits(), '-',
                    names(exbytx[xebyex %>% subjectHits()]))
    xbytx <- xbytx[sname %in% intersect(sname, ename)]
  }
  
  if (long_tx){
    xbytx <- mapLongTX(xbytx, exbytx, tx, exon_only=exon_only)
  } 
  
  if (exon_only) {
    hit_names <- names(exbytx)[subjectHits(xbytx)]
  } else {
    hit_names <- names(tx)[subjectHits(xbytx)]
  }
  
  mapped_names <- hit_names %>% unique()
  mapped_ex <- exbytx[mapped_names] %>% unlist()
  mapped_cds <- cdsbytx[names(cdsbytx) %in% mapped_names] %>% unlist()
  mapped_utr5 <- utr5bytx[names(utr5bytx) %in% mapped_names] %>% unlist()
  mapped_utr3 <- utr3bytx[names(utr3bytx) %in% mapped_names] %>% unlist()
  mapped_in <- inbytx[names(inbytx) %in% mapped_names] %>% unlist()
  
  mapped_tx <- tx[names(tx) %in% mapped_names]
  tx_starts <- start(mapped_tx) %>% rep(each=ngrid)
  stride1 <- (rep(width(mapped_tx), each=ngrid) * 
                rep(c(0:(ngrid-1)) * (1/ngrid), length(mapped_tx))) %>% round()
  stride2 <- (rep(width(mapped_tx), each=ngrid) * 
                rep(c(1:ngrid) * (1/ngrid), length(mapped_tx))) %>% round()
  
  equitx <- GRanges(seqnames=seqnames(mapped_tx) %>% rep(each=ngrid),
                    IRanges(start=(tx_starts+stride1), 
                            end=(tx_starts+stride2-1)),
                    strand=strand(mapped_tx) %>% rep(each=ngrid))
  names(equitx) <- names(mapped_tx) %>% rep(each=ngrid)
  
  equitx$length <- width(equitx)
  
  pexon <- mapGrid(equitx, mapped_ex)
  equitx$exon <- pexon
  
  pintron <- mapGrid(equitx, mapped_in)
  equitx$intron <- pintron
  
  pcds <- mapGrid(equitx, mapped_cds)
  equitx$cds <- pcds
  
  putr5 <- mapGrid(equitx, mapped_utr5)
  equitx$utr5 <- putr5
  
  putr3 <- mapGrid(equitx, mapped_utr3)
  equitx$utr3 <- putr3
  
  equitx <- equitx %>% split(names(equitx))
  equitx <- equitx[hit_names]
  
  mapped_x <- x[xbytx %>% queryHits()]
  names(mapped_x) <- paste0(queryHits(xbytx), '-',
                            names(exbytx)[subjectHits(xbytx)])
  
  equitx <- equitx %>% unlist()
  names(equitx) <- rep(paste0(queryHits(xbytx), '-',
                              names(exbytx)[subjectHits(xbytx)]), each=ngrid)
  
  hits <- findOverlaps(mapped_x, equitx)
  true_hits <- hits[which(names(mapped_x[queryHits(hits)]) ==
                            names(equitx[subjectHits(hits)]))]
  
  target_map <- rep(0, length(equitx))
  target_map[subjectHits(true_hits)] <- 1
  equitx$target <- target_map
  equitx <- equitx %>% split(names(equitx))
  
  if (unify_strand){
    strands <- strand(equitx %>% unlist()) %>% as.vector()
    strands <- strands[seq(1, length(strands), ngrid)]
    minus_idx <- which(strands == "-")
    
    minus_grid <- equitx[minus_idx]
    names(minus_grid) <- NULL
    minus_grid <- minus_grid %>% unlist() %>% sort(decreasing = T) 
    minus_grid <- minus_grid %>% split(names(minus_grid))
    
    equitx <- c(equitx[-minus_idx], minus_grid)
  }
  return(equitx)
}


mapGrid <- function(x, target){
  hits <- findOverlaps(x, target)
  true_hits <- hits[which(names(x[queryHits(hits)]) ==
                            names(target[subjectHits(hits)]))]
  overlaps <- pintersect(x[queryHits(true_hits)], target[subjectHits(true_hits)])
  
  percent <- width(overlaps) / width(x[queryHits(true_hits)])
  percent <- percent %>% split(queryHits(true_hits))
  percent <- lapply(percent, function (x) sum(x)) %>% unlist()
  
  ptarget <- rep(0, length(x))
  ptarget[unique(queryHits(true_hits))] <- percent
  return(ptarget)
}



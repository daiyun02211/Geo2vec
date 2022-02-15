#' Generate chunkTX for single-nucleotide level data.
#' 
#' @description This function accepts GRanges data with width 1 and generate chunkTX for each input.
#' 
#' @param x A \code{GRanges} object 
#' @param txdb A transcript database, currently only TxDb and EnsDb are supported.
#' @param exon_only A \code{logical} object indicating whether to consider only sites from exonic regions.
#' @param long_tx A \code{logical} object indicating whether to consider only the longest transcript for each site.
#' @param mRNA A \code{character} object indicating whether to consider only protein coding transcripts.
#' @return A \code{\link{GRanges}} object whose metacolumns are the generated chunkTX.
#' @examples 
#' 
#' ## Generate chunkTX for example data.
#' library(EnsDb.Hsapiens.v86)
#' input <- import.bed(system.file("extdata", "example.bed", package = "Geo2vec"))
#' encoding <- site_chunkTX(input, EnsDb.Hsapiens.v86, exon_only=T, long_tx=T, mRNA=T)
#' 
#' @import magrittr
#' @export
site_chunkTX <- function(x, txdb, exon_only=TRUE, long_tx=TRUE, mRNA=TRUE){
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
  
  if (exon_only){
    if (mean(overlapsAny(x, exbytx)) != 1){
      if (!mRNA){
        stop("Error: Some targets can not be mapped to exons!")
      }else{
        stop("Error: Some targets can not be mapped to coding exons!")
      }
    }
    xbytx <- findOverlaps(x, exbytx)
  } else {
    xbytx <- findOverlaps(x, tx)
  }
  
  if (long_tx){
    xbytx <- mapLongTX(xbytx, exbytx, tx, exon_only=exon_only)
    if (exon_only){
      long_hits_names <- paste0(xbytx %>% queryHits(), '-',
                                names(exbytx[xbytx %>% subjectHits()]))
    } else {
      long_hits_names <- paste0(xbytx %>% queryHits(), '-',
                                names(tx[xbytx %>% subjectHits()]))
    }
  } else {
    long_hits_names <- NULL
  }
  
  if (exon_only) {
    hit_names <- names(exbytx)[subjectHits(xbytx)]
  } else {
    hit_names <- names(tx)[subjectHits(xbytx)]
  }
  
  mapped_names <- hit_names %>% unique()
  mapped_ex <- exbytx[mapped_names]
  mapped_in <- inbytx[names(inbytx) %in% mapped_names] 
  mapped_cds <- cdsbytx[names(cdsbytx) %in% mapped_names] %>% unlist()
  mcols(mapped_cds) <- NULL
  mapped_utr5 <- utr5bytx[names(utr5bytx) %in% mapped_names] %>% unlist()
  mcols(mapped_utr5) <- NULL
  mapped_utr3 <- utr3bytx[names(utr3bytx) %in% mapped_names] %>% unlist()
  mcols(mapped_utr3) <- NULL
  
  mapped_cds$target <- 0
  mapped_cds$exon <- 1
  mapped_cds$intron <- 0
  mapped_cds$cds <- 1
  mapped_cds$utr3 <- 0
  mapped_cds$utr5 <- 0
  mapped_cds$length <- width(mapped_cds)
  
  mapped_utr5$target <- 0
  mapped_utr5$exon <- 1
  mapped_utr5$intron <- 0
  mapped_utr5$cds <- 0
  mapped_utr5$utr3 <- 0
  mapped_utr5$utr5 <- 1
  mapped_utr5$length <- width(mapped_utr5)
  
  mapped_utr3$target <- 0
  mapped_utr3$exon <- 1
  mapped_utr3$intron <- 0
  mapped_utr3$cds <- 0
  mapped_utr3$utr3 <- 1
  mapped_utr3$utr5 <- 0
  mapped_utr3$length <- width(mapped_utr3)
  
  if (db_type == 'TxDb'){
    selected_in <- mapped_in[names(mapped_in) %in% coding_names] %>% unlist()
  } else {
    if (!mRNA){
      selected_in <- mapped_in[names(mapped_in) %in% coding_names] %>% unlist()
    }else{
      selected_in <- mapped_in %>% unlist()
    }
  }
  
  mcols(selected_in) <- NULL
  selected_in$target <- 0
  selected_in$exon <- 0
  selected_in$intron <- 1
  selected_in$cds <- 0
  selected_in$utr3 <- 0
  selected_in$utr5 <- 0
  selected_in$length <- width(selected_in)
  
  t1_family <- mapChunk(x, mapped_cds, 'cds', mapped_ex, selected_in,
                        mapped_cds, mapped_utr3, mapped_utr5,
                        long_tx, long_hits_names)
  t2_family <- mapChunk(x, mapped_utr3, 'utr3', mapped_ex, selected_in,
                        mapped_cds, mapped_utr3, mapped_utr5,
                        long_tx, long_hits_names)
  t3_family <- mapChunk(x, mapped_utr5, 'utr5', mapped_ex, selected_in,
                        mapped_cds, mapped_utr3, mapped_utr5,
                        long_tx, long_hits_names)
  
  if (!exon_only){
    t0_family <- mapChunk(x, selected_in, 'coding_in', mapped_ex, selected_in,
                          mapped_cds, mapped_utr3, mapped_utr5,
                          long_tx, long_hits_names)
  }
  
  if (!mRNA){
    mapped_ex <- mapped_ex[!(names(mapped_ex) %in% coding_names)]
    mapped_ex <- mapped_ex %>% unlist()
    mcols(mapped_ex) <- NULL
    
    mapped_in <- mapped_in[!(names(mapped_in) %in% coding_names)]
    mapped_in <- mapped_in %>% unlist()
    mcols(mapped_in) <- NULL
    
    mapped_ex$target <- 0
    mapped_ex$exon <- 1
    mapped_ex$intron <- 0
    mapped_ex$cds <- 0
    mapped_ex$utr3 <- 0
    mapped_ex$utr5 <- 0
    mapped_ex$length <- width(mapped_ex)
    
    mapped_in$target <- 0
    mapped_in$exon <- 0
    mapped_in$intron <- 1
    mapped_in$cds <- 0
    mapped_in$utr3 <- 0
    mapped_in$utr5 <- 0
    mapped_in$length <- width(mapped_in)
    
    t4_family <- mapChunk(x, mapped_ex, 'exon', mapped_ex, mapped_in,
                          mapped_cds, mapped_utr3, mapped_utr5,
                          long_tx, long_hits_names)
    
    if (!exon_only){
      t5_family <- mapChunk(x, mapped_in, 'intron', mapped_ex, mapped_in,
                            mapped_cds, mapped_utr3, mapped_utr5,
                            long_tx, long_hits_names)
      out <- c(t0_family, t1_family, t2_family, t3_family, t4_family, t5_family)
    }else{
      out <- c(t1_family, t2_family, t3_family, t4_family)
    }
  }else{
    if (exon_only) out <- c(t1_family, t2_family, t3_family)
    else out <- c(t0_family, t1_family, t2_family, t3_family)
  }
  return(out)
}


mapChunk <- function(x, target, t_type, exon, intron, cds, utr3, utr5,
                     long_tx, long_names=NULL){
  if (!(t_type %in% c('exon', 'coding_in', 'intron', 'cds', 'utr3', 'utr5'))){
    stop("Error: Type should be one of exon, coding_in, intron, cds, utr3 and utr5")
  }
  xbytarget <- findOverlaps(x, target)
  m_name <- paste0(xbytarget %>% queryHits(), '-', 
                   names(target[xbytarget %>% subjectHits()]))
  
  if (long_tx){
    long_idx <- m_name %in% long_names
    m_x <- x[queryHits(xbytarget)[long_idx]]
    m_t <- target[subjectHits(xbytarget)[long_idx]]
    m_name <- m_name[long_idx]
  } else {
    m_x <- x[xbytarget %>% queryHits()]
    m_t <- target[xbytarget %>% subjectHits()]
  }
  
  names(m_x) <- m_name
  
  m_t_s1 <- GRanges(seqnames <- seqnames(m_t),
                    IRanges(start(m_t), 
                            ifelse((start(m_x)-1) > start(m_t), 
                                   start(m_x)-1, start(m_t))),
                    strand <- strand(m_t))
  names(m_t_s1) <- m_name
  m_t_s1 <- m_t_s1[width(m_t_s1) > 0]
  
  m_t_s2 <- GRanges(seqnames <- seqnames(m_t),
                    IRanges(ifelse((start(m_x)+1) <= end(m_t),
                                   start(m_x)+1, end(m_t)),
                            end(m_t)),
                    strand <- strand(m_t))
  names(m_t_s2) <- m_name
  m_t_s2 <- m_t_s2[width(m_t_s2) > 0]
  
  if (length(m_t_s1)){
    m_t_s1$target <- 0
    m_t_s1$exon <- as.integer(!(t_type =='intron' | t_type == 'coding_in'))
    m_t_s1$intron <- as.integer(t_type =='intron' | t_type == 'coding_in')
    m_t_s1$cds <- as.integer(t_type == 'cds')
    m_t_s1$utr3 <- as.integer(t_type == 'utr3')
    m_t_s1$utr5 <- as.integer(t_type == 'utr5')
    m_t_s1$length <- width(m_t_s1)
  }
  
  if (length(m_t_s2)){
    m_t_s2$target <- 0
    m_t_s2$exon <- as.integer(!(t_type =='intron' | t_type == 'coding_in'))
    m_t_s2$intron <- as.integer(t_type =='intron' | t_type == 'coding_in')
    m_t_s2$cds <- as.integer(t_type == 'cds')
    m_t_s2$utr3 <- as.integer(t_type == 'utr3')
    m_t_s2$utr5 <- as.integer(t_type == 'utr5')
    m_t_s2$length <- width(m_t_s2)
  }
  
  if (length(m_x)){
    m_x$target <- 1
    m_x$exon <- as.integer(!(t_type =='intron' | t_type == 'coding_in'))
    m_x$intron <- as.integer(t_type =='intron' | t_type == 'coding_in')
    m_x$cds <- as.integer(t_type == 'cds')
    m_x$utr3 <- as.integer(t_type == 'utr3')
    m_x$utr5 <- as.integer(t_type == 'utr5')
    m_x$length <- width(m_x)
  }
  
  if (t_type %in% c('coding_in', 'cds', 'utr5', 'utr3')){
    family <- c(cds, utr5, utr3, intron)
  } else {
    family <- c(exon, intron)
  }
  
  fname <- names(family)
  names(family) <- NULL
  
  family <- split(family, fname)[names(m_t)]
  names(family) <- m_name

  family <- GRangesList(family)
  family <- unlist(family)
  
  names(m_t) <- m_name
  tmp_hits <- findOverlaps(family, m_t)
  true_hits <- names(family[queryHits(tmp_hits)]) == names(m_t[subjectHits(tmp_hits)])
  
  family <- family[-queryHits(tmp_hits)[true_hits]]
  
  family <- c(family, m_x, m_t_s1, m_t_s2)
  family <- family %>% sort()
  family <- family %>% split(names(family))
  
  return(family)
}


mapLongTX <- function(hits, ex, tx, exon_only=TRUE){
  dup_idx <- which(queryHits(hits) %in% queryHits(hits)[queryHits(hits) %>% duplicated()])
  dup_hits <- as.data.frame(hits[dup_idx, ])
  if (exon_only){
    dup_tx_names <- ex[dup_hits$subjectHits] %>% names()
  } else {
    dup_tx_names <- tx[dup_hits$subjectHits] %>% names()
  }
  tx_width <- tx %>% width()
  tx_names <- tx$tx_name
  names(tx_width) <- tx_names
  dup_hits$tx_width <- tx_width[dup_tx_names]
  rownames(dup_hits) <- dup_idx
  tx_width_vec <- dup_hits$tx_width
  names(tx_width_vec) <- rownames(dup_hits)
  map_idx <- tapply(tx_width_vec, dup_hits$queryHits,
                    function(x) names(x)[which.max(x)])
  hits <- hits[-1 * as.integer(rownames(dup_hits)[!rownames(dup_hits) %in% as.character(map_idx)]),]
  return(hits)
}



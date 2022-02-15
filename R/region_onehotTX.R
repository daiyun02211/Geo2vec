#' Generate onehotTX for region level data.
#' 
#' @description This function accepts \code{GRanges} data with width larger than 1 and generate onehotTX for each input.
#' 
#' @param x A \code{GRanges} object 
#' @param txdb A transcript database, currently only TxDb and EnsDb are supported.
#' @param window An integer indicates the width of flanking region which equals to 2*window + region width.
#' @param exon_only A \code{logical} object indicating whether to consider only regions whose starts and ends are from exons.
#' @param long_tx A \code{logical} object indicating whether to consider only the longest transcript for each region.
#' @param mRNA A \code{character} object indicating whether to consider only protein coding transcripts.
#' @return A \code{\link{GRanges}} object whose metacolumns are the generated onehotTX.
#' @examples 
#' 
#' ## Generate onehotTX for example data.
#' library(EnsDb.Hsapiens.v86)
#' input <- import.bed(system.file("extdata", "example.bed", package = "Geo2vec"))
#' encoding <- region_onehotTX(input+10, EnsDb.Hsapiens.v86, window=0, exon_only=T, long_tx=T, mRNA=T)
#' 
#' @import magrittr
#' @export
region_onehotTX <- function(x, txdb, window=0, exon_only=TRUE, long_tx=TRUE, mRNA=TRUE){
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
  
  x_flank <- x + window
  
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
  
  t1_family <- region_mapOH(x, x_flank, mapped_cds, 'cds', 
                            tx, mapped_ex, selected_in,
                            mapped_cds, mapped_utr3, mapped_utr5,
                            long_tx, long_hits_names, exon_only)
  t2_family <- region_mapOH(x, x_flank, mapped_utr3, 'utr3', 
                            tx, mapped_ex, selected_in,
                            mapped_cds, mapped_utr3, mapped_utr5,
                            long_tx, long_hits_names, exon_only)
  t3_family <- region_mapOH(x, x_flank, mapped_utr5, 'utr5', 
                            tx, mapped_ex, selected_in,
                            mapped_cds, mapped_utr3, mapped_utr5,
                            long_tx, long_hits_names, exon_only)
  
  if (!exon_only){
    t0_family <- region_mapOH(x, x_flank, selected_in, 'coding_in', 
                              tx, mapped_ex, selected_in,
                              mapped_cds, mapped_utr3, mapped_utr5,
                              long_tx, long_hits_names, exon_only)
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
    
    t4_family <- region_mapOH(x, x_flank, mapped_ex, 'exon',
                              tx, mapped_ex, mapped_in,
                              mapped_cds, mapped_utr3, mapped_utr5,
                              long_tx, long_hits_names, exon_only)
    
    if (!exon_only){
      t5_family <- region_mapOH(x, x_flank, mapped_in, 'intron', 
                                tx, mapped_ex, mapped_in,
                                mapped_cds, mapped_utr3, mapped_utr5,
                                long_tx, long_hits_names, exon_only)
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


region_mapOH <- function(x, x_flank, target, t_type, tx,
                         exon, intron, cds, utr3, utr5,
                         long_tx, long_names=NULL, exon_only=TRUE){
  if (!(t_type %in% c('exon', 'coding_in', 'intron', 'cds', 'utr3', 'utr5'))){
    stop("Error: Type should be one of exon, coding_in, intron, cds, utr3 and utr5")
  }
  x_start <- GRanges(seqnames=seqnames(x),
                     IRanges(start=start(x)),
                     strand=strand(x))
  xbytarget <- findOverlaps(x_start, target)
  
  if (!exon_only){
    xbytarget <- xbytarget[end(x[xbytarget %>% queryHits()]) <=
                             end(tx[names(target[xbytarget %>% subjectHits()])])]
    m_name <- paste0(xbytarget %>% queryHits(), '-', 
                      names(target[xbytarget %>% subjectHits()]))
  } else {
    x_end <- GRanges(seqnames=seqnames(x),
                     IRanges(start=end(x)),
                     strand=strand(x))
    xebyex <- findOverlaps(x_end, exon)
    
    m_name <- paste0(xbytarget %>% queryHits(), '-', 
                      names(target[xbytarget %>% subjectHits()]))
    ename <- paste0(xebyex %>% queryHits(), '-',
                    names(exon[xebyex %>% subjectHits()]))
    
    xbytarget <- xbytarget[m_name %in% intersect(m_name, ename)]
    m_name <- m_name[m_name %in% intersect(m_name, ename)]
  }
  
  if (long_tx){
    long_idx <- m_name %in% long_names
    m_x <- x[queryHits(xbytarget)[long_idx]]
    m_fx <- x_flank[queryHits(xbytarget)[long_idx]]
    m_t <- target[subjectHits(xbytarget)[long_idx]]
    m_name <- m_name[long_idx]
  } else {
    m_x <- x[xbytarget %>% queryHits()]
    m_fx <- x_flank[queryHits(xbytarget)]
    m_t <- target[xbytarget %>% subjectHits()]
  }
  
  m_tx <- tx[names(m_t)]
  left_out_idx <- which(start(m_fx) < start(m_tx))
  start(m_fx[left_out_idx]) <- start(m_tx[left_out_idx])
  right_out_idx <- which(end(m_fx) > end(m_tx))
  end(m_fx[right_out_idx]) <- end(m_tx[right_out_idx])
  
  # we may add mcols padding here
  
  names(m_x) <- m_name
  
  m_t_s1 <- GRanges(seqnames <- seqnames(m_t),
                    IRanges(start=ifelse(start(m_fx) > start(m_t),
                                         start(m_fx), start(m_t)),
                            end=start(m_x)-1),
                    strand <- strand(m_t))
  names(m_t_s1) <- m_name
  m_t_s1 <- m_t_s1[width(m_t_s1) > 0]
  
  m_t_s2 <- GRanges(seqnames <- seqnames(m_t),
                    IRanges(start=start(m_x)+1,
                            end=ifelse(end(m_fx) < end(m_t),
                                       end(m_fx), end(m_t))),
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
  
  left_ridx <- which(start(m_t) > start(m_fx))
  right_ridx <- which(end(m_t) < end(m_fx))
  
  if (length(left_ridx) > 0){
    left_remain <- GRanges(seqnames=seqnames(m_fx[left_ridx]),
                           IRanges(start(m_fx[left_ridx]),
                                   start(m_t[left_ridx])-1),
                           strand=strand(m_fx[left_ridx]))
    names(left_remain) <- m_name[left_ridx]
    
    all_sm_t <- GRanges()
    all_sm_t_in <- GRanges()
    smid_x <- GRanges()
    
    if (t_type %in% c('coding_in', 'cds', 'utr5', 'utr3')){
      targets <- list(cds, utr5, utr3, intron)
      tnames <- c('cds', 'utr5', 'utr3', 'intron')
    } else {
      targets <- list(exon, intron)
      tnames <- c('exon', 'intron')
    }
    
    for (i in seq_along(targets)){
      target <- targets[[i]]
      tname <- tnames[i]
      
      x <- GRanges(seqnames=seqnames(left_remain),
                   IRanges(start=start(left_remain)),
                   strand=strand(left_remain))
      names(x) <- names(left_remain)
      xbytarget <- findOverlaps(x, target)
      
      if (length(xbytarget) == 0) next
      
      smt_name <- names(target[xbytarget %>% subjectHits()])
      smx_name <- names(x[xbytarget %>% queryHits()]) %>%
        strsplit(split='-') %>% unlist()
      
      smx_name <- smx_name[seq(2, length(smx_name), 2)]
      xbytarget <- xbytarget[which(smx_name == smt_name)]
      
      if (length(xbytarget) == 0) next
      
      sm_name <- names(x[xbytarget %>% queryHits()])
      sm_x <- left_remain[xbytarget %>% queryHits()]
      sm_t <- target[xbytarget %>% subjectHits()]
      names(sm_t) <- sm_name
      
      sm_t_in <- GRanges(seqnames <- seqnames(sm_t),
                         IRanges(start(sm_x), end(sm_t)),
                         strand <- strand(sm_t))
      names(sm_t_in) <- sm_name
      sm_t_in <- sm_t_in[width(sm_t_in) > 0]
      
      if (length(sm_t_in) > 0){
        sm_t_in$target <- 0
        sm_t_in$exon <- as.integer(tname != 'intron')
        sm_t_in$intron <- as.integer(tname == 'intron')
        sm_t_in$cds <- as.integer(tname == 'cds')
        sm_t_in$utr3 <- as.integer(tname == 'utr3')
        sm_t_in$utr5 <- as.integer(tname == 'utr5')
        sm_t_in$length <- width(sm_t_in)
      }
      
      remain_idx <- which(end(sm_x) > end(sm_t))
      smid_remain <- GRanges(seqnames=seqnames(sm_x[remain_idx]),
                             IRanges(start=end(sm_t[remain_idx])+1,
                                     end=end(sm_x[remain_idx])),
                             strand=strand(sm_x[remain_idx]))
      names(smid_remain) <- names(sm_x[remain_idx])
      smid_x <- c(smid_x, smid_remain)
      all_sm_t <- c(all_sm_t, sm_t)
      all_sm_t_in <- c(all_sm_t_in, sm_t_in)
      left_remain <- left_remain[-which(names(left_remain) %in% sm_name)]
    }
    
    if (length(left_remain) != 0){
      stop('Some of input can not be fully mapped to TXs')
    }
    
    smm_t_in <- GRanges()
    if (length(smid_x) > 0){
      if (t_type %in% c('coding_in', 'cds', 'utr5', 'utr3')){
        target <- c(cds, utr5, utr3, intron)
      } else {
        target <- c(exon, intron)
      }
      
      regionbymid <- findOverlaps(target, smid_x)
      mr_name <- names(target[regionbymid %>% queryHits()])
      mmid_name <- names(smid_x[regionbymid %>% subjectHits()]) %>%
        strsplit(split='-') %>% unlist()
      mmid_name <- mmid_name[seq(2, length(mmid_name), 2)]
      regionbymid <- regionbymid[which(mr_name == mmid_name)]
      
      mr_in <- target[regionbymid %>% queryHits()]
      names(mr_in) <- names(smid_x[regionbymid %>% subjectHits])
      smm_t_in <- c(smm_t_in, mr_in)
    }
  } else {
    all_sm_t <- GRanges()
    all_sm_t_in <- GRanges()
    smm_t_in <- GRanges()
  }
  
  if (length(right_ridx) > 0){
    right_remain <- GRanges(seqnames=seqnames(m_fx[right_ridx]),
                            IRanges(end(m_t[right_ridx]+1),
                                    end(m_fx[right_ridx])),
                            strand=strand(m_fx[right_ridx]))
    names(right_remain) <- m_name[right_ridx]
    
    all_em_t <- GRanges()
    all_em_t_in <- GRanges()
    emid_x <- GRanges()
    
    if (t_type %in% c('coding_in', 'cds', 'utr5', 'utr3')){
      targets <- list(cds, utr5, utr3, intron)
      tnames <- c('cds', 'utr5', 'utr3', 'intron')
    } else {
      targets <- list(exon, intron)
      tnames <- c('exon', 'intron')
    }
    
    for (i in seq_along(targets)){
      target <- targets[[i]]
      tname <- tnames[i]
      
      x <- GRanges(seqnames=seqnames(right_remain),
                   IRanges(start=end(right_remain), width=1),
                   strand=strand(right_remain))
      names(x) <- names(right_remain)
      xbytarget <- findOverlaps(x, target)
      
      if (length(xbytarget) == 0) next
      
      emt_name <- names(target[xbytarget %>% subjectHits()])
      emx_name <- names(x[xbytarget %>% queryHits()]) %>%
        strsplit(split='-') %>% unlist()
      
      emx_name <- emx_name[seq(2, length(emx_name), 2)]
      xbytarget <- xbytarget[which(emx_name == emt_name)]
      
      if (length(xbytarget) == 0) next
      
      em_name <- names(x[xbytarget %>% queryHits()])
      em_x <- right_remain[xbytarget %>% queryHits()]
      em_t <- target[xbytarget %>% subjectHits()]
      names(em_t) <- em_name
      
      em_t_in <- GRanges(seqnames <- seqnames(em_t),
                         IRanges(start(em_t), end(em_x)),
                         strand <- strand(em_t))
      names(em_t_in) <- em_name
      em_t_in <- em_t_in[width(em_t_in) > 0]
      
      if (length(em_t_in) > 0){
        em_t_in$target <- 0
        em_t_in$exon <- as.integer(tname != 'intron')
        em_t_in$intron <- as.integer(tname == 'intron')
        em_t_in$cds <- as.integer(tname == 'cds')
        em_t_in$utr3 <- as.integer(tname == 'utr3')
        em_t_in$utr5 <- as.integer(tname == 'utr5')
        em_t_in$length <- width(em_t_in)
      }
      
      remain_idx <- which(start(em_x) < start(em_t))
      emid_remain <- GRanges(seqnames=seqnames(em_x[remain_idx]),
                             IRanges(start=start(em_x[remain_idx]),
                                     end=start(em_t[remain_idx])-1),
                             strand=strand(em_x[remain_idx]))
      names(emid_remain) <- names(em_x[remain_idx])
      emid_x <- c(emid_x, emid_remain)
      all_em_t <- c(all_em_t, em_t)
      all_em_t_in <- c(all_em_t_in, em_t_in)
      right_remain <- right_remain[-which(names(right_remain) %in% em_name)]
    }
    
    if (length(right_remain) != 0){
      stop('Some of input can not be fully mapped to TXs')
    }
    
    emm_t_in <- GRanges()
    if (length(emid_x) > 0){
      if (t_type %in% c('coding_in', 'cds', 'utr5', 'utr3')){
        target <- c(cds, utr5, utr3, intron)
      } else {
        target <- c(exon, intron)
      }
      
      regionbymid <- findOverlaps(target, emid_x)
      mr_name <- names(target[regionbymid %>% queryHits()])
      mmid_name <- names(emid_x[regionbymid %>% subjectHits()]) %>%
        strsplit(split='-') %>% unlist()
      mmid_name <- mmid_name[seq(2, length(mmid_name), 2)]
      regionbymid <- regionbymid[which(mr_name == mmid_name)]
      
      mr_in <- target[regionbymid %>% queryHits()]
      names(mr_in) <- names(emid_x[regionbymid %>% subjectHits])
      emm_t_in <- c(emm_t_in, mr_in)
    }
  } else {
    all_em_t <- GRanges()
    all_em_t_in <- GRanges()
    emm_t_in <- GRanges()
  }
  
  family <- c(m_x, m_t_s1, m_t_s2,
              smm_t_in, all_sm_t_in,
              emm_t_in, all_em_t_in)
  family <- family %>% sort()
  family <- family %>% split(names(family))
  
  return(family)
  
}


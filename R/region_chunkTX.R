#' Generate chunkTX for region level data.
#' 
#' @description This function accepts \code{GRanges} data with width larger than 1 and generate chunkTX for each input.
#' 
#' @param x A \code{GRanges} object. 
#' @param txdb A transcript database, currently only TxDb and EnsDb are supported.
#' @param exon_only A \code{logical} object indicating whether to consider only regions whose starts and ends are from exons.
#' @param long_tx A \code{logical} object indicating whether to consider only the longest transcript for each region.
#' @param mRNA A \code{character} object indicating whether to consider only protein coding transcripts.
#' @return A \code{\link{GRanges}} object whose metacolumns are the generated chunkTX.
#' @examples 
#' 
#' ## Generate chunkTX for example data.
#' library(EnsDb.Hsapiens.v86)
#' input <- import.bed(system.file("extdata", "example.bed", package = "Geo2vec"))
#' encoding <- region_chunkTX(input+10, EnsDb.Hsapiens.v86, exon_only=T, long_tx=T, mRNA=T)
#' 
#' @import magrittr
#' @export
region_chunkTX <- function(x, txdb, exon_only=TRUE, long_tx=TRUE, mRNA=TRUE){
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
  } else {
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
  
  if (exon_only) {
    hit_names <- names(exbytx)[subjectHits(xbytx)]
  } else {
    hit_names <- names(tx)[subjectHits(xbytx)]
  }
  
  mapped_names <- hit_names %>% unique()
  mapped_cds <- cdsbytx[names(cdsbytx) %in% mapped_names] %>% unlist()
  mcols(mapped_cds) <- NULL
  mapped_utr5 <- utr5bytx[names(utr5bytx) %in% mapped_names] %>% unlist()
  mcols(mapped_utr5) <- NULL
  mapped_utr3 <- utr3bytx[names(utr3bytx) %in% mapped_names] %>% unlist()
  mcols(mapped_utr3) <- NULL
  
  mapped_ex <- exbytx[names(exbytx) %in% mapped_names]
  mapped_in <- inbytx[names(inbytx) %in% mapped_names]
  
  selected_ex <- mapped_ex[names(mapped_ex) %in% coding_names] %>% unlist()
  mcols(selected_ex) <- NULL
  
  selected_in <- mapped_in[names(mapped_in) %in% coding_names] %>% unlist()
  mcols(selected_in) <- NULL
  
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
  
  selected_in$target <- 0
  selected_in$exon <- 0
  selected_in$intron <- 1
  selected_in$cds <- 0
  selected_in$utr3 <- 0
  selected_in$utr5 <- 0
  selected_in$length <- width(selected_in)
  
  t1_family <- region_mapChunk(x, mapped_cds, 'cds', tx,
                               selected_ex, selected_in,
                               mapped_cds, mapped_utr3, mapped_utr5,
                               long_tx, long_hits_names, exon_only)
  t2_family <- region_mapChunk(x, mapped_utr3, 'utr3', tx,
                               selected_ex, selected_in,
                               mapped_cds, mapped_utr3, mapped_utr5,
                               long_tx, long_hits_names, exon_only)
  t3_family <- region_mapChunk(x, mapped_utr5, 'utr5', tx, 
                               selected_ex, selected_in, 
                               mapped_cds, mapped_utr3, mapped_utr5, 
                               long_tx, long_hits_names, exon_only)
  
  # whether the Intronic boundaries of the region shall be considered 
  if (!exon_only){
    t0_family <- region_mapChunk(x, selected_in, 'coding_in', tx, 
                                 selected_ex, selected_in, 
                                 mapped_cds, mapped_utr3, mapped_utr5,
                                 long_tx, long_hits_names, exon_only)
  }
  
  # whether the transcript without CDS/UTR annotation shall be considered
  if (!mRNA){
    mapped_ex <- mapped_ex[!(names(mapped_ex) %in% coding_names)]
    mapped_ex <- mapped_ex %>% unlist()
    mcols(mapped_ex) <- NULL
    
    mapped_ex$target <- 0
    mapped_ex$exon <- 1
    mapped_ex$intron <- 0
    mapped_ex$cds <- 0
    mapped_ex$utr3 <- 0
    mapped_ex$utr5 <- 0
    mapped_ex$length <- width(mapped_ex)
    
    mapped_in <- mapped_in[!(names(mapped_in) %in% coding_names)]
    mapped_in <- mapped_in %>% unlist()
    mcols(mapped_in) <- NULL
    
    mapped_in$target <- 0
    mapped_in$exon <- 0
    mapped_in$intron <- 1
    mapped_in$cds <- 0
    mapped_in$utr3 <- 0
    mapped_in$utr5 <- 0
    mapped_in$length <- width(mapped_in)
    
    t4_family <- region_mapChunk(x, mapped_ex, 'exon', tx, 
                                 mapped_ex, mapped_in,
                                 mapped_cds, mapped_utr3, mapped_utr5,
                                 long_tx, long_hits_names, exon_only)
    
    if (!exon_only){
      t5_family <- region_mapChunk(x, mapped_in, 'intron', tx, 
                                   mapped_ex, mapped_in, 
                                   mapped_cds, mapped_utr3, mapped_utr5,
                                   long_tx, long_hits_names, exon_only)
      out <- c(t0_family, t1_family, t2_family, t3_family, t4_family, t5_family)
    }else{
      out <- c(t1_family, t2_family, t3_family, t4_family)
    }
  } else{
    if (exon_only) out <- c(t1_family, t2_family, t3_family)
    else out <- c(t0_family, t1_family, t2_family, t3_family)
  }
  return(out)
}


region_mapChunk <- function(x, target, t_type, tx, exon, intron,
                            cds, utr3, utr5, long_tx, long_names=NULL,
                            exon_only=TRUE){
  if (!(t_type %in% c('exon', 'coding_in', 'intron', 'cds', 'utr3', 'utr5'))){
    stop("Error: Type should be one of exon, coding_in, intron, cds, utr3 and utr5")
  }
  x_start <- GRanges(seqnames=seqnames(x),
                     IRanges(start=start(x)),
                     strand=strand(x))
  xbytarget <- findOverlaps(x_start, target)
  
  if (length(xbytarget) == 0){
    return(c())
  }else{
    if (!exon_only){
      xbytarget <- xbytarget[end(x[xbytarget %>% queryHits()]) <=
                               end(tx[names(target[xbytarget %>% subjectHits()])])]
      sm_name <- paste0(xbytarget %>% queryHits(), '-', 
                        names(target[xbytarget %>% subjectHits()]))
    } else {
      x_end <- GRanges(seqnames=seqnames(x),
                       IRanges(start=end(x)),
                       strand=strand(x))
      xebyex <- findOverlaps(x_end, exon)
      
      sm_name <- paste0(xbytarget %>% queryHits(), '-', 
                        names(target[xbytarget %>% subjectHits()]))
      ename <- paste0(xebyex %>% queryHits(), '-',
                      names(exon[xebyex %>% subjectHits()]))
      
      xbytarget <- xbytarget[sm_name %in% intersect(sm_name, ename)]
      sm_name <- sm_name[sm_name %in% intersect(sm_name, ename)]
    }
    
    if (long_tx){
      long_idx <- sm_name %in% long_names
      sm_x <- x[queryHits(xbytarget)[long_idx]]
      sm_t <- target[subjectHits(xbytarget)[long_idx]]
      sm_name <- sm_name[long_idx]
    } else {
      sm_x <- x[xbytarget %>% queryHits()]
      sm_t <- target[xbytarget %>% subjectHits()]
    }
    
    sm_t_out1 <- GRanges(seqnames <- seqnames(sm_t),
                         IRanges(start(sm_t),
                                 ifelse((start(sm_x)-1) > start(sm_t), 
                                        start(sm_x)-1, start(sm_t))),
                         strand <- strand(sm_t))
    names(sm_t_out1) <- sm_name
    sm_t_out1 <- sm_t_out1[width(sm_t_out1) > 0]
    
    if (length(sm_t_out1) > 0){
      sm_t_out1$target <- 0
      sm_t_out1$exon <- as.integer(!(t_type =='intron' | t_type == 'coding_in'))
      sm_t_out1$intron <- as.integer(t_type =='intron' | t_type == 'coding_in')
      sm_t_out1$cds <- as.integer(t_type == 'cds')
      sm_t_out1$utr3 <- as.integer(t_type == 'utr3')
      sm_t_out1$utr5 <- as.integer(t_type == 'utr5')
      sm_t_out1$length <- width(sm_t_out1)
    }
    
    sm_t_in <- GRanges(seqnames <- seqnames(sm_t),
                       IRanges(start(sm_x),
                               ifelse(end(sm_t)<end(sm_x), end(sm_t), end(sm_x))),
                       strand <- strand(sm_t))
    names(sm_t_in) <- sm_name
    sm_t_in <- sm_t_in[width(sm_t_in) > 0]
    
    sm_t_in$target <- 1
    sm_t_in$exon <- as.integer(!(t_type =='intron' | t_type == 'coding_in'))
    sm_t_in$intron <- as.integer(t_type =='intron' | t_type == 'coding_in')
    sm_t_in$cds <- as.integer(t_type == 'cds')
    sm_t_in$utr3 <- as.integer(t_type == 'utr3')
    sm_t_in$utr5 <- as.integer(t_type == 'utr5')
    sm_t_in$length <- width(sm_t_in)
    
    sm_t_out2 <- GRanges(seqnames <- seqnames(sm_t),
                         IRanges(ifelse(end(sm_x)<end(sm_t), end(sm_x)+1, end(sm_t)+1),
                                 end(sm_t)),
                         strand <- strand(sm_t))
    names(sm_t_out2) <- sm_name
    sm_t_out2 <- sm_t_out2[width(sm_t_out2) > 0]
    
    if (length(sm_t_out2) > 0){
      sm_t_out2$target <- 0
      sm_t_out2$exon <- as.integer(!(t_type =='intron' | t_type == 'coding_in'))
      sm_t_out2$intron <- as.integer(t_type =='intron' | t_type == 'coding_in')
      sm_t_out2$cds <- as.integer(t_type == 'cds')
      sm_t_out2$utr3 <- as.integer(t_type == 'utr3')
      sm_t_out2$utr5 <- as.integer(t_type == 'utr5')
      sm_t_out2$length <- width(sm_t_out2)
    }
    
    remaining <- sm_x[which(end(sm_t) < end(sm_x))]
    
    if (length(remaining) > 0){
      remaining <- GRanges(seqnames=seqnames(remaining),
                           IRanges(end(sm_t[which(end(sm_t)<end(sm_x))])+1,
                                   end(remaining)),
                           strand=strand(remaining))
      names(remaining) <- sm_name[which(end(sm_t)<end(sm_x))]
      
      all_em_t <- GRanges()
      all_em_t_in <- GRanges()
      all_em_t_out <- GRanges()
      mid_x <- GRanges()
      
      if (t_type %in% c('cds', 'utr5', 'utr3', 'coding_in')){
        if (!exon_only){
          targets <- list(cds, utr5, utr3, intron)
          tnames <- c('cds', 'utr5', 'utr3', 'intron')
        } else {
          targets <- list(cds, utr5, utr3)
          tnames <- c('cds', 'utr5', 'utr3')
        }
      } else {
        if (!exon_only){
          targets <- list(exon, intron)
          tnames <- c('exon', 'intron')
        } else {
          targets <- list(exon)
          tnames <- c('exon')
        }
      }
      
      for (i in seq_along(targets)){
        target <- targets[[i]]
        tname <- tnames[i]
        
        x <- GRanges(seqnames=seqnames(remaining),
                     IRanges(start=end(remaining)),
                     strand=strand(remaining))
        names(x) <- names(remaining)
        xbytarget <- findOverlaps(x, target)
        
        if (length(xbytarget) == 0) next
        
        emt_name <- names(target[xbytarget %>% subjectHits()])
        emx_name <- names(x[xbytarget %>% queryHits()]) %>%
          strsplit(split='-') %>% unlist()
        
        emx_name <- emx_name[seq(2, length(emx_name), 2)]
        xbytarget <- xbytarget[which(emx_name == emt_name)]
        
        if (length(xbytarget) == 0) next
        
        em_name <- names(x[xbytarget %>% queryHits()])
        
        em_x <- remaining[xbytarget %>% queryHits()]
        em_t <- target[xbytarget %>% subjectHits()]
        
        names(em_t) <- em_name
        
        em_t_in <- GRanges(seqnames <- seqnames(em_t),
                           IRanges(start(em_t), end(em_x)),
                           strand <- strand(em_t))
        names(em_t_in) <- em_name
        em_t_in <- em_t_in[width(em_t_in) > 0]
        
        if (length(em_t_in) > 0){
          em_t_in$target <- 1
          em_t_in$exon <- as.integer(tname != 'intron')
          em_t_in$intron <- as.integer(tname == 'intron')
          em_t_in$cds <- as.integer(tname == 'cds')
          em_t_in$utr3 <- as.integer(tname == 'utr3')
          em_t_in$utr5 <- as.integer(tname == 'utr5')
          em_t_in$length <- width(em_t_in)
        }
        
        em_t_out <- GRanges(seqnames <- seqnames(em_t),
                            IRanges(ifelse(end(em_x)<end(em_t), end(em_x)+1, end(em_t)),
                                    end(em_t)),
                            strand <- strand(em_t))
        names(em_t_out) <- em_name
        em_t_out <- em_t_out[width(em_t_out) > 0]
        
        if (length(em_t_out) > 0){
          em_t_out$target <- 0
          em_t_out$exon <- as.integer(tname != 'intron')
          em_t_out$intron <- as.integer(tname == 'intron')
          em_t_out$cds <- as.integer(tname == 'cds')
          em_t_out$utr3 <- as.integer(tname == 'utr3')
          em_t_out$utr5 <- as.integer(tname == 'utr5')
          em_t_out$length <- width(em_t_out)
        }
        
        remain_idx <- which(start(em_x) < start(em_t))
        mid_remain <- GRanges(seqnames=seqnames(em_x[remain_idx]),
                              IRanges(start=start(em_x[remain_idx]),
                                      end=start(em_t[remain_idx])-1),
                              strand=strand(em_x[remain_idx]))
        names(mid_remain) <- names(em_x[remain_idx])
        mid_x <- c(mid_x, mid_remain)
        all_em_t <- c(all_em_t, em_t)
        all_em_t_in <- c(all_em_t_in, em_t_in)
        all_em_t_out <- c(all_em_t_out, em_t_out)
        remaining <- remaining[-which(names(remaining) %in% em_name)]
      }
      
      if (length(remaining) != 0) stop('Some of input can not be fully mapped to TXs')
      
      mm_t_in <- GRanges()
      
      if (length(mid_x) > 0){
        if (t_type %in% c('cds', 'utr5', 'utr3', 'coding_in')){
          target <- c(cds, utr5, utr3, intron)
        } else {
          target <- c(exon, intron)
        }
        
        regionbymid <- findOverlaps(target, mid_x)
        mr_name <- names(target[regionbymid %>% queryHits()])
        mmid_name <- names(mid_x[regionbymid %>% subjectHits()]) %>%
          strsplit(split='-') %>% unlist()
        mmid_name <- mmid_name[seq(2, length(mmid_name), 2)]
        regionbymid <- regionbymid[which(mr_name == mmid_name)]
        
        mr_in <- target[regionbymid %>% queryHits()]
        names(mr_in) <- names(mid_x[regionbymid %>% subjectHits])
        mr_in$target <- 1
        mm_t_in <- c(mm_t_in, mr_in)
      }
    } else {
      all_em_t <- GRanges()
      all_em_t_in <- GRanges()
      all_em_t_out <- GRanges()
      mm_t_in <- GRanges()
    }
    
    if (t_type %in% c('cds', 'utr5', 'utr3', 'pintron')){
      family <- c(cds, utr5, utr3, intron)
    } else {
      family <- c(exon, intron)
    }
    fname <- names(family)
    names(family) <- NULL
    
    mt_name <- names(sm_t)
    family <- split(family, fname)[mt_name]
    names(family) <- sm_name
    
    family <- family %>% unlist()
    names(sm_t) <- sm_name
    
    redundance <- c(sm_t, all_em_t, mm_t_in)
    
    tmp_hits <- findOverlaps(family, redundance, type='equal')
    true_hits <- names(family[queryHits(tmp_hits)]) ==
      names(redundance[subjectHits(tmp_hits)])
    
    family <- family[-queryHits(tmp_hits)[true_hits]]
    
    family <- c(family, sm_t_in, sm_t_out1, sm_t_out2,
                all_em_t_in, all_em_t_out, mm_t_in)
    family <- family %>% sort()
    family <- family %>% split(names(family))
    
    return(family)
  }
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
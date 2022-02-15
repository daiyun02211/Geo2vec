#' Resize chunkTX.
#' 
#' @description This function resizes chunkTX centered on the target site/region.
#' 
#' @param chunkTX A \code{GRanges} object whose metadata columns are generated chunkTX. 
#' @param window An integer indicates the target size which equals to 2*window + 1.
#' @param unify_strand A \code{logical} object indicating whether to unify the direction from 5'end to 3'end.
#' @param region A \code{logical} object indicating whether to resize region chunkTX.
#' @return A \code{data.frame} object contains the resized chunkTX. 
#'  
#' @examples
#' \donttest{
#' ## Resizing a large number of chunkTX is time-consuming.
#' ## Therefore, mcapply from package parallel is suggested.
#' 
#' lapply(chunkTX, resizeChunk, window=17, region=FALSE)
#' mcapply(chunkTX, resizeChunk, window=17, region=FALSE, mc.cores=8)
#' }
#' @import GenomicRanges
#' @export
resizeChunk <- function(chunkTX, window=17, unify_strand=TRUE, region=FALSE){
  meta <- as.data.frame(mcols(chunkTX))
  
  psize <- 2*window+1
  
  pad <- data.frame("target" = rep(0, psize), "exon" = rep(0, psize),
                    "intron" = rep(0, psize), "cds" = rep(0, psize),
                    "utr3" = rep(0, psize), "utr5" = rep(0, psize),
                    "length" = rep(0, psize))
  
  if (region){
    tidx <- meta$target == 1
    ltarget <- as.integer(sum(tidx) / 2) 
    rtarget <- sum(tidx) - ltarget
    ltstart <- (window + 1) - ltarget
    rtend <- window + rtarget
    
    if (unify_strand){
      if (as.vector(strand(chunkTX)@values) == '+'){
        pad[max(1, ltstart - (min(which(tidx)) - 1)):
              min(psize, rtend + nrow(meta) - max(which(tidx))),] = 
          meta[max(1, min(which(tidx)) - (ltstart - 1)):
                 min(nrow(meta), max(which(tidx)) + psize - rtend),]
      } else {
        pad[max(1, ltstart - (nrow(meta) - max(which(tidx)))):
              min(psize, rtend + min(which(tidx) - 1)),] = 
          meta[min(max(which(tidx)) + (ltstart - 1), nrow(meta)):
                 max(1, min(which(tidx)) - (psize - rtend)),]
      }
    } else {
      pad[max(1, ltstart - (min(which(tidx)) - 1)):
            min(psize, rtend + nrow(meta) - max(which(tidx))),] = 
        meta[max(1, min(which(tidx)) - (ltstart - 1)):
               min(nrow(meta), max(which(tidx)) + psize - rtend),]
    }
  }else{
    idx <- which(meta$target == 1)
    if (unify_strand){
      if (strand(chunk)@values %>% as.vector() == '+'){
        pad[max(window-idx+2, 1):min(2*window+1, window+1+nrow(meta)-idx),] = 
          meta[max(idx-window, 1):min(idx+window, nrow(meta)),]
      } else {
        pad[max(window+1-nrow(meta)+idx, 1):min(2*window+1, window+idx),] = 
          meta[min(idx+window, nrow(meta)):max(idx-window, 1),]
      }
    } else {
      pad[max(window-idx+2, 1):min(2*window+1, window+1+nrow(meta)-idx),] = 
        meta[max(idx-window, 1):min(idx+window, nrow(meta)),]
    }
  }
  return(pad)
}
#' Resize onehotTX
#' 
#' @description This function resizes onehotTX centered on the target site.
#' 
#' @param onehotTX A \code{GRanges} object whose metadata columns are generated onehotTX. 
#' @param window An integer indicates the target size which equals to 2*window + 1.
#' @param unify_strand A \code{logical} object indicating whether to unify the direction from 5'end to 3'end.
#' @return A \code{data.frame} object contains the resized onehotTX. 
#'  
#' @examples
#' \donttest{
#' ## Resizing a large number of onehotTX is time-consuming.
#' ## Therefore, mcapply from package parallel is suggested.
#' 
#' lapply(onehotTX, resizeOH, window=250)
#' mcapply(onehotTX, resizeOH, window=250, mc.cores=8)
#' }
#' @import GenomicRanges
#' @export
resizeOH <- function(onehotTX, window=250, unify_strand=T){
  meta <- onehotTX %>% mcols() %>% as.data.frame()
  num_chunk <- nrow(meta)
  mat <- matrix(0, nrow=2*window+1, ncol=5, byrow=T)
  idx <- which(meta$target == 1)
  mat[window+1,] <- meta[idx, 2:6] %>% as.numeric()
  if (idx > 1){
    for (i in c((idx-1):1)){
      num_rep <- meta[i, 7]
      gap <- sum(meta[c(idx:(i+1)), 7])
      mat[(window-gap-num_rep+2):(window-gap+1),] <- meta[i, 2:6] %>% as.numeric() %>% rep(num_rep) %>% matrix(nrow=num_rep, byrow=T)
    }
  }
  if (idx < num_chunk){
    for (i in c((idx+1):num_chunk)){
      num_rep <- meta[i, 7]
      gap <- sum(meta[c(idx:(i-1)), 7])
      mat[(window+1+gap):(window+gap+num_rep),] <- meta[i, 2:6] %>% as.numeric() %>% rep(num_rep) %>% matrix(nrow=num_rep, byrow=T)
    }
  }
  if (unify_strand & strand(onehot)@values %>% as.vector() == '-'){
    return(mat[nrow(mat):1,])
  } else {
    return(mat)
  }
}
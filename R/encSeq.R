#' @title Extract sequence from \code{GRanges} and generate encodings.
#'
#' @description A function used to get sequence features from \code{\link{GRanges}}
#' 
#' @param GRanges A \code{\link{GRanges}} object provides coordinates of interest.
#' 
#' @param BSgenome A \code{\link{BSgenome}} object containing the sequence of the genome.
#' 
#' @param window An integer value indicates the width of flanking region upstream and downstream.
#' 
#' @param type A character indicates the type of sequence encodings. Currently, only 'token' and 'onehot' are supported.
#' 
#' @return A matrix contains the sequence encodings of regions from GRanges - window to GRanges + window.
#' @examples 
#' 
#' ## Generate sequence encoding for example data.
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' input <- import.bed(system.file("extdata", "example.bed", package = "Geo2vec"))
#' encoding <- encSeq(input, BSgenome.Hsapiens.UCSC.hg38, window=250, type='token')
#'
#' @importFrom Biostrings getSeq strsplit
#' @importFrom stats model.matrix
#' @export
encSeq <- function(GRanges, BSgenome, window = 0, type='token'){
  
  stopifnot(class(GRanges) =="GRanges")
  
  seq <- getSeq(BSgenome, GRanges + window) %>% as.character()
  
  seq_df <- do.call(rbind, strsplit(seq,"")) %>% data.frame() %>% 
    lapply(factor, levels = c("A", "C", "G", "T")) %>% data.frame()
  
  options(na.action='na.pass')
  
  if (type == 'token') {
    seq_out <- data.matrix(seq_df)
    if(anyNA(seq_out)) seq_out[is.na(seq_out)] <- 0
  } else if (type == 'onehot') {
    seq_out <- model.matrix(~ -1 + ., seq_df, contrasts.arg = 
                              lapply(seq_df, contrasts, contrasts=FALSE))
    if(anyNA(seq_out)) seq_out[is.na(seq_out)] <- 0
  }
  
  return(seq_out)
}
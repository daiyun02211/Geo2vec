#' Convert GRanges to data frame.
#' 
#' @description This function converts \code{GRanges} objects containing Geo2vec encodings to data frames.
#' 
#' @param GRanges A \code{GRanges} object whose metadata columns are generated Geo2vec encodings. 
#' @param org_site_id A set of index indicating which entries in the original input Geo2vec came from.
#' @return A \code{data.frame} object contains the Geo2vec encodings. 
#' 
#' @examples
#' \donttest{
#' ## Convert chunkTX to dataframe for saving.
#' 
#' chunkTX <- site_chunkTX(input, txdb, exon_only = T, long_tx = T, mRNA = T)
#' out_df <- geo2df(chunkTX)
#' } 
#' @import magrittr
#' @import GenomicRanges
#' @export
geo2df <- function(GRanges, org_site_id=NULL){
  if (class(GRanges) == 'GRanges'){
  } else {
    GRanges <- unlist(GRanges)
  }
  out_names <- names(GRanges) %>% strsplit('[.]') %>% unlist()
  out_names <- out_names[seq(1, length(out_names), 2)]
  if (!is.null(org_site_id)){
    out_site_map <- out_names %>% strsplit('-') %>% unlist()
    out_site_id <- out_site_map[seq(1, length(out_site_map), 2)] %>% as.integer()
    reindex_site_id <- org_site_id[out_site_id]
    out_tx_id <- out_site_map[seq(2, length(out_site_map), 2)]
    out_names <- paste0(reindex_site_id, '-', out_tx_id)
  }
  range_df <- data.frame(mapping=out_names,
                         seqnames=seqnames(GRanges),
                         start=start(GRanges),
                         end=end(GRanges),
                         strand=strand(GRanges))
  out_df <- mcols(GRanges)
  out_df <- cbind(range_df, out_df)
  out_df %>% head()
  row.names(out_df) <- 1:nrow(out_df)
  return(out_df)
}

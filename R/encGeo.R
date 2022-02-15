#' Generate Geo2vec encodings.
#' 
#' @description This function accepts GRanges data and generate Geo2vec encodings for each input.
#' 
#' @param input A \code{GRanges} object contains genomic coordinates.
#' @param txdb A transcript database, currently only TxDb and EnsDb are supported.
#' @param exon_only A \code{logical} object indicating whether to consider only sites from exonic regions.
#' @param long_tx A \code{logical} object indicating whether to consider only the longest transcript for each site.
#' @param mRNA A \code{character} object indicating whether to consider only protein coding transcripts.
#' @param ngrid An integer indicating the number of grids, i.e., how many pieces will the transcript be divided into.
#' @param window An integer indicates the width of flanking region which equals to 2*window + input width.
#' @param type An \code{character} indicating the type of Geo2vec, should be one of the chunkTX, gridTX, landmarkTX, onehotTX.
#' @param unify_strand A \code{logical} object indicating whether to unify the direction from 5'end to 3'end for gridTX. For chunkTX and gridTX, unify_strand is inside the resize function.
#' @return A \code{\link{GRanges}} object whose metacolumns are the generated Geo2vec encodings.
#' @examples 
#' 
#' library(EnsDb.Hsapiens.v86)
#' input <- import.bed(system.file("extdata", "example.bed", package = "Geo2vec"))
#' 
#' ## Generate chunkTX for example data.
#' encoding <- encGeo(input, EnsDb.Hsapiens.v86, type='chunkTX', exon_only=T, long_tx=T, mRNA=T)
#' 
#' ## Generate gridTX for example data.
#' encoding <- encGeo(input, EnsDb.Hsapiens.v86, type='gridTX', ngrid=40, exon_only=T, long_tx=T, mRNA=T)
#' 
#' ## Generate landmarkTX for example data.
#' encoding <- encGeo(input, EnsDb.Hsapiens.v86, type='landmarkTX', long_tx=T, mRNA=T)
#' 
#' ## Generate onehotTX for example data.
#' encoding <- encGeo(input, EnsDb.Hsapiens.v86, type='onehotTX', window=250, exon_only=T, long_tx=T, mRNA=T)
#' 
#' @import magrittr
#' @export
encGeo <- function(input, txdb, exon_only = TRUE, long_tx = TRUE, 
                   mRNA = TRUE, ngrid = 40, window = 0,
                   type = 'chunkTX', unify_strand = TRUE){
  if (!(type %in% c('chunkTX', 'gridTX', 'landmarkTX', 'onehotTX'))){
    stop("Error: Encoding type should be one of the following: chunkTX, gridTX, landmarkTX, and onehotTX.")
  }
  
  db_type <- class(txdb)[1]
  if (!(db_type %in% c('TxDb', 'EnsDb'))){
    stop("Error: Currently only TxDb and EnsDb are supported.")
  }
  
  if (type == 'gridTX'){
    if ((ngrid %% 1 == 0) & (ngrid > 0)){} else {
      stop("Error: ngrid should be a positive integer.")
    }
  }
  
  if (type == 'onehotTX'){
    if ((window %% 1 == 0) & (window > 0)){} else {
      stop("Error: window should be a positive integer.")
    }
  }
  
  if (mean(width(input) == 1)){
    if (type == 'chunkTX'){
      encoding <- site_chunkTX(input, txdb, exon_only, long_tx, mRNA)
    } else if (type == 'gridTX'){
      encoding <- gridTX(input, txdb, ngrid, exon_only, long_tx, mRNA, unify_strand)
    } else if (type == 'landmarkTX'){
      encoding <- site_landmarkTX(input, txdb, long_tx, mRNA)
    } else if (type == 'onehotTX'){
      encoding  <- site_onehotTX(input, txdb, window, exon_only, long_tx, mRNA)
    }
  } else if (min(width(input) > 1)){
    if (type == 'chunkTX'){
      encoding <- region_chunkTX(input, txdb, exon_only, long_tx, mRNA)
    } else if (type == 'gridTX'){
      encoding <- gridTX(input, txdb, ngrid, exon_only, long_tx, mRNA, unify_strand)
    } else if (type == 'onehotTX'){
      encoding <- region_onehotTX(input, txdb, window, exon_only, long_tx, mRNA)
    } else if (type == 'landmarkTX'){
      stop("Error: Currently landmarkTX only supports single-nucleotide level data.")
    }
  } else {
    site_idx <- which(width(input) == 1)
    region_idx <- which(width(input) > 1)
    if (type == 'chunkTX'){
      site_enc <- site_chunkTX(input[site_idx], txdb, exon_only, long_tx, mRNA)
      region_enc <- region_chunkTX(input[region_idx], txdb, exon_only, long_tx, mRNA)
      encoding <- c(site_enc, region_enc)
    } else if (type == 'onehotTX'){
      site_enc <- site_onehotTX(input[site_idx], txdb, window, exon_only, long_tx, mRNA)
      region_enc <- region_onehotTX(input[site_idx], txdb, window, exon_only, long_tx, mRNA)
      encoding <- c(site_enc, region_enc)
    } else if (type == 'gridTX'){
      encoding <- gridTX(input, txdb, ngrid, exon_only, long_tx, mRNA, unify_strand)
    } else if (type == 'landmarkTX'){
      stop("Error: Currently landmarkTX only supports single-nucleotide level data.")
    }
  }
  return (encoding)
}
    
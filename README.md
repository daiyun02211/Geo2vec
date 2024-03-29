# Geo2vec
## Introduction
**Geo2vec** (**geo**graphic representation of transcript as **vec**tors) explored different strategies for encoding sub-molecular geographic information of ribonucleotides. Three encoding methods, i.e., landmarkTX, gridTX, and chunkTX, as well as the widely used one-hot method are currently supported. LandmarkTX is a lightweight encoding scheme directly capturing the position of the target ribonucleotide (or site) relative to transcript landmarks, i.e., the distances to the two edges of the exon, coding sequence (CDS), and transcript, respectively. Meanwhile, gridTX and chunkTX are designed to describe the landscape of the entire transcript through grids (of equal widths) or regions (with unequal width), respectively. 
## Installation
To install Geo2vec from Github, please use the following command in R consol.
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("daiyun02211/Geo2vec")
```
## Usage
Example m6A coordinates can be found in inst/extdata:
```
input <- import.bed(system.file("extdata", "example.bed", package = "Geo2vec"))
```
It is recommended to use the function **encGeo** to generate the encoding. Different encodings can be selected by the parameter *type*:
```
library(EnsDb.Hsapiens.v86)
encoding <- encGeo(input, EnsDb.Hsapiens.v86, type='onehotTX', window=50, exon_only=T, long_tx=T, mRNA=T)
encoding <- encGeo(input, EnsDb.Hsapiens.v86, type='landmarkTX', long_tx=T, mRNA=T)
encoding <- encGeo(input, EnsDb.Hsapiens.v86, type='gridTX', ngrid=40, exon_only=T, long_tx=T, mRNA=T)
encoding <- encGeo(input, EnsDb.Hsapiens.v86, type='chunkTX', exon_only=T, long_tx=T, mRNA=T)
```
The current version of the package supports transcription annotation packages in the form of TxDb (e.g., TxDb.Hsapiens.UCSC.hg19.knownGene) and EnsDb.

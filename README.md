# Geo2vec
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
encoding <- encGeo(input, EnsDb.Hsapiens.v86, type='onehotTX', window=50, exon_only=T, long_tx=T, mRNA=T)
encoding <- encGeo(input, EnsDb.Hsapiens.v86, type='landmarkTX', long_tx=T, mRNA=T)
encoding <- encGeo(input, EnsDb.Hsapiens.v86, type='gridTX', ngrid=40, exon_only=T, long_tx=T, mRNA=T)
encoding <- encGeo(input, EnsDb.Hsapiens.v86, type='chunkTX', exon_only=T, long_tx=T, mRNA=T)
```

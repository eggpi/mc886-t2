library("jpeg")
library("parallel")

BDR_DIR = "cropped/bdr"
BDC_DIR = "cropped/bdc"

psapply <- function(data, fn) {
    return(simplify2array(mclapply(data, fn, mc.cores = 8)))
}

euclid.dist <- function(v1, v2) {
    return(sqrt(sum((v1 - v2) ^ 2)))
}

make.feature.vectors.from.images <- function(image.files) {
    fvs <- psapply(image.files, function(imf) {
        return(c(t(readJPEG(imf))))
    })

    fvs <- t(fvs)
    rownames(fvs) <- names(image.files)
    return(fvs)
}

load.feature.vectors.from.image.dir <- function(dir) {
    image.files <- Sys.glob(file.path(dir, "*.jpg"))
    names(image.files) <- basename(image.files)
    return(make.feature.vectors.from.images(image.files))
}

bdc <- load.feature.vectors.from.image.dir(BDC_DIR)
bdr <- load.feature.vectors.from.image.dir(BDR_DIR)

library("jpeg")
library("parallel")
library("FNN")
library("MASS")

BDR_DIR = "cropped/bdr"
BDC_DIR = "cropped/bdc"

# the variance to retain when selecting principal components
PCA_CUMULATIVE_VARIANCE = 0.95

psapply <- function(data, fn) {
    return(simplify2array(mclapply(data, fn, mc.cores = 8)))
}

get.person.from.image.name <- function(image.name) {
    return(substr(image.name, 1, 5))
}

make.feature.vectors.from.images <- function(image.files) {
    fvs <- t(psapply(image.files, function(imf) {
        return(c(t(readJPEG(imf))))
    }))

    # turn fvs into matrix
    dim <- dim(fvs)
    fvs <- matrix(fvs, nrow = dim[1], ncol = dim[2])
    return(fvs)
}

load.feature.vectors.from.image.dir <- function(dir) {
    image.files <- Sys.glob(file.path(dir, "*.jpg"))
    fvs <- make.feature.vectors.from.images(image.files)
    rownames(fvs) <- psapply(image.files, basename)
    return(fvs)
}

# our custom distance
custom.dist <- function(v1, v2) {
    return(sum(abs(v1 - v2)))
}

# query a fv in bdr.
# returns the index of the row in bdr that best fits the fv.
custom.query.one.fv <- function(fv, bdr) {
    bdr.images <- rownames(bdr)
    distances <- sapply(bdr.images, function(imname) {
        return(custom.dist(fv, bdr[imname,]))
    })

    return(match(min(distances), distances))
}

# query bdr for the vectors in bdc using a custom distance function.
# returns a vector of the indices in bdr for the closest match for
# each fv in bdc.
custom.query.bdr <- function(bdc, bdr) {
    bdr.images <- rownames(bdr)
    return(psapply(rownames(bdc), function(imname) {
        return(custom.query.one.fv(bdc[imname,], bdr))
    }))
}

# query bdr for the vectors in bdc using euclidean distance
euclidean.query.bdr <- function(bdc, bdr) {
    return(get.knnx(bdr, bdc, k = 1, algorithm = "kd_tree")[[1]])
}

compute.correct.classifications <- function(query.results, bdc, bdr) {
    return(sum(psapply(1:(dim(bdc)[1]), function(bdc.idx) {
        bdr.idx <- query.results[bdc.idx]

        queried.image <- rownames(bdc)[bdc.idx]
        queried.person <- get.person.from.image.name(queried.image)

        result.image <- rownames(bdr)[bdr.idx]
        result.person <- get.person.from.image.name(result.image)
        return(result.person == queried.person)
    })))
}

apply.pca <- function(bdc, bdr) {
    # FIXME should we apply PCA on all feature vectors
    # instead of bdr only?
    pca.result <- princomp(bdr, scale = TRUE)

    variances <- pca.result$sdev ^ 2
    variance.proportions <- variances / sum(variances)
    cumulative.variances <- cumsum(variance.proportions)
    components <- which(cumulative.variances > PCA_CUMULATIVE_VARIANCE)[1]

    rotation <- loadings(pca.result)[,1:components]
    bdr <- bdr %*% rotation
    bdc <- bdc %*% rotation

    return(list(bdc, bdr))
}

apply.lda <- function(bdc, bdr) {
    # FIXME should we apply LDA on all feature vectors
    # instead of bdr only?
    grouping <- psapply(rownames(bdr), get.person.from.image.name)
    lda.results <- lda(bdr, grouping)
    scaling <- lda.results$scaling
    bdr <- bdr %*% scaling
    bdc <- bdc %*% scaling

    return(list(bdc, bdr))
}

query.bdr <- function(bdc, bdr, dist) {
    if (dist == "euclidean") {
        query.results <- euclidean.query.bdr(bdc, bdr)
    } else {
        query.results <- custom.query.bdr(bdc, bdr)
    }

    return(compute.correct.classifications(query.results, bdc, bdr))
}

results <- rep(0, 5)
names(results) <- c("euclid,plain", "euclid,pca", "euclid,lda",
                    "euclid,lda+pca", "custom,lda+pca")

bdc <- load.feature.vectors.from.image.dir(BDC_DIR)
bdr <- load.feature.vectors.from.image.dir(BDR_DIR)

print(query.bdr(bdc, bdr, "custom"))

results["euclid,plain"] <- query.bdr(bdc, bdr, "euclidean")

pca.results <- apply.pca(bdc, bdr)
bdc.pca <- pca.results[[1]]
bdr.pca <- pca.results[[2]]

results["euclid,pca"] <- query.bdr(bdc.pca, bdr.pca, "euclidean")

lda.results <- apply.lda(bdc, bdr)
bdc.lda <- lda.results[[1]]
bdr.lda <- lda.results[[2]]

results["euclid,lda"] <- query.bdr(bdc.lda, bdr.lda, "euclidean")

lda.pca.results <- apply.lda(bdc.pca, bdr.pca)
bdc.lda.pca <- lda.pca.results[[1]]
bdr.lda.pca <- lda.pca.results[[2]]

results["euclid,lda+pca"] <- query.bdr(bdc.lda.pca, bdr.lda.pca, "euclidean")
results["custom,lda+pca"] <- query.bdr(bdc.lda.pca, bdr.lda.pca, "custom")

#' Title
#'
#' @param dagPeptides 
#' @param dagBackground 
#' @param group 
#' @param bgNoise 
#'
#' @return
#' @export
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' 
#' 
testDAU <- function(dagPeptides,
                    dagBackground,
                    group = c("null", "classic", "charge", "chemistry", "hydrophobicity"),
                    bgNoise = NA) 
{
    if (missing(dagPeptides) || class(dagPeptides) != "dagPeptides") {
        stop(
            "dagPeptides should be an object of dagPeptides class.\n
            Please try ?fetchSequence to get help.",
            call. = FALSE
        )
    }
    if (missing(dagBackground) ||
        class(dagBackground) != "dagBackground") {
        stop(
            "dagBackground should be an object of dagBackground class .\n
            Please try ?buildBackgroundModel to get help.",
            call. = FALSE
        )
    }
    if (!is.na(bgNoise)) {
        if (bgNoise < 0 || bgNoise > 1)
            stop("The background noise is a number in the range of 0 to 1", 
                 call. = FALSE)
    }
    group <- match.arg(group)
    
    exp <- dagPeptides@peptides
    bg <- dagBackground@background
    
    AA <- get("namehash", envir = cacheEnv)
    classic <- get("classic", envir = cacheEnv)$group
    charge <- get("charge", envir = cacheEnv)$group
    chemistry <-get("chemistry", envir = cacheEnv)$group
    hydrophobicity <- get("hydrophobicity", envir = cacheEnv)$group
    
    coln <- if (group == "null") as.character(AA) else names(get(group))
    
    ##convert by group
    convert <- function(x, gtype) 
    {
        d <- dim(x)
        x <- as.character(x)
        for (i in 1:length(gtype)) 
        {
            id <- x %in% gtype[[i]]
            name <- names(gtype)[i]
            x[id] <- name
        }
        matrix(x, nrow = d[1], ncol = d[2])
    }
    groupAA <- function(dat, group) 
    {
        dat <- switch(
            group,
            classic = convert(dat, classic),
            charge = convert(dat, charge),
            hydrophobicity = convert(dat, hydrophobicity),
            chemistry = convert(dat, chemistry),
            null = dat
        )
        dat
    }
    bg <- lapply(bg, function(.bg)
        groupAA(.bg, group))
    exp <- groupAA(exp, group)
    if (ncol(exp) != ncol(bg[[1]]))
    {
        stop("The length of background is different from inputs",
             call. = FALSE)
    }
    counts <- function(mat, coln) 
    {
        num <- apply(mat, 2, function(.ele) {
            cnt <- table(.ele)[coln]
            names(cnt) <- coln ## just in case not all coln in the dataset
            cnt[is.na(cnt)] <- 0
            total <- sum(cnt)
            cnt / total})
    }
    bg <- lapply(bg, counts, coln)
    exp <- counts(exp, coln)
    rownames(exp) <- coln
    bg <- lapply(1:ncol(exp), function(i) {
        do.call(cbind, lapply(bg, function(.bg) {
            .bg[, i]
        }))
    })
    
    ## Add background noise to the background models
    if (!is.na(bgNoise)) 
    {
        rdirichlet <- function (n, alpha) {
            l <- length(alpha)
            x <-
                matrix(rgamma(l * n, alpha),
                       ncol = n,
                       byrow = TRUE)
            sm <- x %*% rep(1, n)
            t(x / as.vector(sm))
        }
        bg <- lapply(bg, function(col_freqs) {
            (1 - bgNoise) * col_freqs + bgNoise * rdirichlet(nrow(col_freqs), rep(1, ncol(col_freqs)))
        })
    }
    
    ##Z-score = (x-mu)/std
    ## The standard deviation is not correctly calculated?
    std <- do.call(cbind, lapply(bg, function(.bg) {
        apply(.bg, 1, sd, na.rm = TRUE)
    }))
    mu <- do.call(cbind, lapply(bg, function(.bg) {
        apply(.bg, 1, mean, na.rm = TRUE)
    }))
    
    
    ##difference
    exp[is.na(exp)] <- 0
    diff <- exp - mu
    diff[is.na(diff)] <- 0
    
    zscore <- diff / std
    
    rownames(diff) <- coln
    rownames(zscore) <- coln
    coln <- c()
    if (dagPeptides@upstreamOffset > 0) 
    {
        coln <- paste("AA", -1 * (dagPeptides@upstreamOffset:1), sep = "")
    }
    coln <- c(coln, "AA0")
    if (dagPeptides@downstreamOffset > 0) 
    {
        coln <- c(coln,
                  paste("AA", 1:dagPeptides@downstreamOffset, sep = ""))
    }
    if (length(coln) == ncol(diff)) 
    {
        colnames(diff) <- colnames(zscore) <- coln
    } 
    else
    {
        colnames(diff) <-
            colnames(zscore) <- paste("AA", 1:ncol(diff), sep = "")
    }
    pvalue <- 2 * pnorm(-abs(zscore))
    
    ## return the test results as an object of testDAUresults class
    new(
        "testDAUresults",
        group = group,
        difference = diff,
        zscore = zscore,
        pvalue = pvalue,
        background = mu,
        motif = exp,
        upstream = dagPeptides@upstreamOffset,
        downstream = dagPeptides@downstreamOffset
    )
}
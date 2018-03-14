#' Conduct differential usage test of amino acids.
#'
#' Test differential usage of amino acids with or without grouping in
#' experimental sets and background sets.
#'
#' @param dagPeptides An object of Class \code{\link{dagPeptides}}.
#' @param dagBackground An object of Class \code{\link{dagBackground}}.
#' @param testType A character vector of length 1. The available options are
#' "fisher" and "z-test", that is "Fisher's exact test" and "Z-test". When the 
#' the number of sequences in the background set is too samll to perform non-
#' replacement subsamplint, "Fisher's exact test" is suggested.
#' @param groupingScheme A character vector of length 1. Available choices are
#' "no", "classic", "charge", "chemistry",and "hydrophobicity". It is used
#' to group amino acids into groups.
#' @param bgNoise A numeric vector of length 1. It should be in the interval of
#' (0, 1).
#'
#' @return An object of Class \code{\link{testDAUresults}}.
#' @export
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' 
#' 
testDAU <- function(dagPeptides,
                    dagBackground,
                    testType = c("fisher", "z-test"),
                    groupingScheme = c("no", "classic", "charge", 
                                       "chemistry", "hydrophobicity"),
                    bgNoise = NA) 
{
    if (missing(dagPeptides) || class(dagPeptides) != "dagPeptides") 
    {
        stop(
            "dagPeptides should be an object of dagPeptides class.\n
            Please try ?fetchSequence to get help.",
            call. = FALSE)
    }
    if (missing(dagBackground) || class(dagBackground) != "dagBackground") 
    {
        stop(
            "dagBackground should be an object of dagBackground class .\n
            Please try ?buildBackgroundModel to get help.",
            call. = FALSE)
    }
    if (!is.na(bgNoise)) 
    {
        if (bgNoise < 0 || bgNoise > 1)
        {
            stop("The background noise is a number in the range of 0 to 1", 
                 call. = FALSE)
        }
    }
    testType <- match.arg(testType)
    if(!testType %in% c("fisher", "z-test"))
    {
        stop("The type of testing should be Fisher's exact test or Z-test.", 
             call. = FALSE)
    }

    exp <- dagPeptides@peptides
    bg <- dagBackground@background
    
    groupingScheme <- match.arg(groupingScheme)
    AA <- get("no", envir = cachedEnv)$symbol
    
    if(groupingScheme == "no")
    {
        coln <- as.character(AA)
    } else 
    {
        coln <- names(get(groupingScheme, envir = cachedEnv)$group)    
    } 
    
    ## helper function used to group AAs based on groupingScheme
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
            no = dat
        )
        dat
    }
    
    bg <- lapply(bg, function(.bg)
        groupAA(.bg, groupingScheme))
    
    exp <- groupAA(exp, groupingScheme)
    if (ncol(exp) != ncol(bg[[1]]))
    {
        stop("The length of background is different from inputs", call. = FALSE)
    }
    
    ## get counts for each amino acid at each position if Fisher's exact test
    ## otherwise get proportions
    counts <- function(mat, coln, testType) 
    {
        if (testType == "fisher")
        {
            num <- apply(mat, 2, function(.ele) {
                cnt <- table(.ele)[coln]
                ## just in case not all coln in the dataset
                names(cnt) <- coln 
                cnt[is.na(cnt)] <- 0
                cnt})
        } else
        {
            num <- apply(mat, 2, function(.ele) {
                cnt <- table(.ele)[coln]
                ## just in case not all coln in the dataset
                names(cnt) <- coln 
                cnt[is.na(cnt)] <- 0
                total <- sum(cnt)
                cnt / total})
        }
    }
    bg <- lapply(bg, counts, coln, testType = testType)
    exp <- counts(exp, coln, testType = testType)
    exp[is.na(exp)] <- 0
    
    rownames(exp) <- coln
    
    
    bg <- lapply(1:ncol(exp), function(i) {
        do.call(cbind, lapply(bg, function(.bg) {
            .bg[, i]
        }))
    })
    
    ## Add background noise to the background model
    if (!is.na(bgNoise) && testType == "z-test") 
    {
        rdirichlet <- function (n, alpha) {
            l <- length(alpha)
            x <-matrix(rgamma(l * n, alpha),
                       ncol = n,
                       byrow = TRUE)
            sm <- x %*% rep(1, n)
            t(x / as.vector(sm))
        }
        bg <- lapply(bg, function(col_freqs) {
            (1 - bgNoise) * col_freqs + bgNoise * 
                rdirichlet(nrow(col_freqs), rep(1, ncol(col_freqs)))
        })
    }
    
    ##Z-score = (x-mu)/std
    ## The standard deviation is not correctly calculated?
    mu <- do.call(cbind, lapply(bg, function(.bg) {
        apply(.bg, 1, mean, na.rm = TRUE)
    }))
    
    if (testType == "z-test")
    {
        ## std <- do.call(cbind, lapply(bg, function(.bg) {
        ##    apply(.bg, 1, sd, na.rm = TRUE)}))
        ## A better estimate of sigma = ssqrt(p*(1-p)/n)
        std <- mu
        n <-  dagBackground@numSubsamples
        std <- mu*(1-mu)/n
    }

    
    ##difference

    diff <- exp - mu
    diff[is.na(diff)] <- 0
    if (testType == "z-test")
    {
        zscore <- diff / std
        pvalue <- 2 * pnorm(-abs(zscore))
        oddsRatio = NULL
    } else
    {
        
        
        non_bg <- sweep(x= bg, MARGIN= 2, STATS = colSums(bg), FUN = "-") 
        non_exp <- sweep(x= exp, MARGIN= 2, STATS = colSums(exp), FUN = "-")
        
        testResults <- mapply(function(x, X, y, Y){
            fisher <- fisher.test(matrix(c(x, X, y, Y), nrow =2, byrow =TRUE),
                                  alternative = "two.sided")
            c(pvalues=fisher$p.value, statistics = fisher$estimate)
        })
        
        testOut <- do.call(rbind, testResults)
        pvalue <- matrix(testOut$pvalues, nrow = nrow(exp), byrow =FALSE)
        oddsRatio <-  matrix(testOut$pvalues, nrow = nrow(exp), byrow =FALSE)
        zscore  <- NULL
        
    }
   
    
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
    
    
    ## return the test results as an object of testDAUresults class
    new(
        "testDAUresults",
        group = groupingScheme,
        difference = diff,
        zscore = zscore,
        oddsRatio = oddsRatio,
        pvalue = pvalue,
        background = mu,
        motif = exp,
        upstream = dagPeptides@upstreamOffset,
        downstream = dagPeptides@downstreamOffset
    )
}
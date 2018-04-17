#' Add a custom grouping scheme.
#' 
#' Add a custom grouping scheme for grouping amino acids as desired.
#' 
#' @param color A vector of character, with length match the number of desired 
#' groups and names as the same as the groups' names. This vector specifies the 
#' different colors for visualizing the different groups of amino acids.
#' @param symbol A vector of character, with length match the number of desired 
#' groups and names as the same as the groups' names. This vector specifies the 
#' different symbols for visualizing the different groups of amino acids.
#' @param group A list with names as the same as the groups' names.
#'
#' @return Add the custom grouping scheme to the cached environment.
#' @export
#' 
#'
#' @examples
#' ## Add a grouping scheme based on the BLOSUM50 level 3 
#' color = c(LVIMC = "#33FF00", AGSTP = "#CCFF00",
#'          FYW = '#00FF66', EDNQKRH = "#FF0066")
#' symbol = c(LVIMC = "L", AGSTP = "A", FYW = "F", EDNQKRH = "E")
#' group = list(
#'    LVIMC = c("L", "V", "I", "M", "C"), 
#'    AGSTP = c("A", "G", "S", "T", "P"),
#'    FYW = c("F", "Y", "W"),
#'    EDNQKRH = c("E", "D", "N", "Q", "K", "R", "H"))
#' addGroupingScheme(color = color, symbol = symbol, group = group) 
#' 
 
addGroupingScheme <- function(color = vector("character"), 
                              symbol = vector("character"),
                              group=list())
{
    if (length(color) < 1 || length(symbol) < 1 || length(group) < 1)
    {
        stop("Too few groups!")
    }
    if (length(color) != length(symbol) || length(color) != length(group) ||
        length(symbol) != length(group))
    {
        stop("Wrong grouping specification:", 
             "The length of color, symbol and group should be the same!")
    }
    if(names(color) != names(symbol) || names(color) != names(group) ||
       names(symbol) != names(group))
    {
        stop("The names of color, symbol and group should be the same!") 
    }
    custom <- list(color = color, symbol = symbol, group = group)
    assign("custom", custome, envir = cachedEnv)
}



#' Differential usage test of amino acids or groups.
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
#' "no", "classic", "charge", "chemistry", "hydrophobicity", "BLOSM50_L1", 
#' "BLOSM50_L2", "BLOSM50_L3", "BLOSM50_L4", "BLOSM50_L5", "BLOSM50_L6", 
#' "BLOSM50_L7", "BLOSM50_L8", "chemistry_property_Mahler", 
#' "contact_potential_Maiorov", "protein_blocks_Rogov", 
#' "sequence_alignment_Dayhoff", "structure_alignments_Mirny",
#' and "custom". If "custom" is used, users must define a grouping scheme using a
#' a list containing sublist named as "color", "symbol" and "group" using the 
#' function addGroupingScheme. It is used to group amino acids into groups of
#' similarities.
#' 
#' @param bgNoise A numeric vector of length 1. It should be in the interval of
#' (0, 1).
#'
#' @return An object of Class \code{\link{testDAUresults}}.
#' @import stats
#' @import methods
#' @export
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' dat <- unlist(read.delim(system.file(
#'                                    "extdata", "grB.txt", package = "dagLogo"),
#'                          header = F, as.is = TRUE))
#' head(dat)
#' ##prepare proteome from a fasta file
#' proteome <- prepareProteome(fasta = system.file("extdata",
#'                                                 "HUMAN.fasta",
#'                                                 package = "dagLogo"), 
#'                             species = "Homo sapiens")
#' ##prepare object of dagPeptides
#' ##prepare an object of dagPeptides
#' seq <- formatSequence(seq = dat, proteome = proteome, upstreamOffset = 14,
#'                      downstreamOffset = 15)
#' bg_fisher <- buildBackgroundModel(seq, background = "wholeProteome", 
#'                                   proteome = proteome, testType = "fisher")
#' bg_ztest <- buildBackgroundModel(seq, background = "wholeProteome",
#'                                    proteome = proteome, testType = "ztest")
#'  
#' t0 <- testDAU(seq, dagBackground = bg_ztest)
#' ## grouded classically: nonpolar_aliphatic = c("A", "G", "L", "M", "I", "V"), 
#' ## polar_uncharged = c("C", "P", "Q", "S", "T"), aromatic = c("F", "W", "Y"), 
#' ## positively_charged = c("H", "K", "N", "R"), negatively_charged = c("D", "E").
#' t1 <- testDAU(dagPeptides = seq, dagBackground = bg_ztest, 
#'              groupingScheme = "classic")

#' ## grouded on the basis of charge: positive = c("H", "K", "R"), 
#' ## neutral = c("A", "C", "F", "G", "I", "L", "M", "N", "P", "Q", "S", "T", 
#' ## "V","W","Y"), negative = c("D", "E")
#' t2 <- testDAU(dagPeptides = seq, dagBackground = bg_ztest, 
#'               groupingScheme = "charge")

#' ## grouped on the basis of their chemical property: hydrophobic = c("A", "F",
#' ## "I", "L", "M", "P", "V", "W"), polar = c("C", "G", "S", "T", "Y"), 
#' ## basic = c("H", "K", "R"), neutral = c("N", "Q"), acidic = c("D", "E")
#' t3 <- testDAU(dagPeptides = seq, dagBackground = bg_ztest, 
#'               groupingScheme = "chemistry")

#' ## grouped on the basis of hydrophobicity: hydrophilic = c("D", "E", "K", "N", 
#' ## "Q", "R"), neutral = c("A", "G", "H", "P", "S", "T"), hydrophobic = c("C", 
#' ## "F", "I", "L", "M", "V", "W", "Y")
#' t4 <- testDAU(dagPeptides = seq, dagBackground = bg_ztest, 
#'               groupingScheme = "hydrophobicity")                                   


testDAU <- function(dagPeptides,
                    dagBackground,
                    groupingScheme = c("no", "classic", "charge", 
                                       "chemistry", "hydrophobicity", 
                                       "BLOSM50_L1", "BLOSM50_L2", 
                                       "BLOSM50_L3", "BLOSM50_L4", 
                                       "BLOSM50_L5", "BLOSM50_L6", 
                                       "BLOSM50_L7", "BLOSM50_L8", 
                                       "chemistry_property_Mahler", 
                                       "contact_potential_Maiorov", 
                                       "protein_blocks_Rogov", 
                                       "sequence_alignment_Dayhoff", 
                                       "structure_alignments_Mirny",
                                       "custom"),
                    bgNoise = NA) 
{
    if (missing(dagPeptides) || class(dagPeptides) != "dagPeptides") 
    {
        stop("dagPeptides should be an object of dagPeptides class.\n
            Please try ?fetchSequence to get help.",
            call. = FALSE)
    }
    if (missing(dagBackground) || class(dagBackground) != "dagBackground") 
    {
        stop("dagBackground should be an object of dagBackground class .\n
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

    exp <- dagPeptides@peptides
    bg <- dagBackground@background
    
    groupingScheme <- match.arg(groupingScheme)
    AA <- get("no", envir = cachedEnv)$symbol
    
    if (groupingScheme == "no")
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
            classic = convert(dat, get("classic", envir = cachedEnv)$group),
            charge = convert(dat, get("charge", envir = cachedEnv)$group),
            hydrophobicity = convert(dat, get("hydrophobicity", envir = cachedEnv)$group),
            chemistry = convert(dat, get("chemistry", envir = cachedEnv)$group),
            BLOSM50_L1 = get("BLOSM50_L1", envir = cachedEnv)$group,
            BLOSM50_L2 = get("BLOSM50_L2", envir = cachedEnv)$group,
            BLOSM50_L3 = get("BLOSM50_L3", envir = cachedEnv)$group,
            BLOSM50_L4 = get("BLOSM50_L4", envir = cachedEnv)$group,
            BLOSM50_L5 = get("BLOSM50_L5", envir = cachedEnv)$group,
            BLOSM50_L6 = get("BLOSM50_L6", envir = cachedEnv)$group,
            BLOSM50_L7 = get("BLOSM50_L7", envir = cachedEnv)$group,
            BLOSM50_L8 = get("BLOSM50_L8", envir = cachedEnv)$group,
            chemistry_property_Mahler = get("chemistry_property_Mahler", 
                                            envir = cachedEnv)$group, 
            contact_potential_Maiorov = get("contact_potential_Maiorov",
                                            envir = cachedEnv)$group, 
            protein_blocks_Rogov = get("protein_blocks_Rogov", 
                                       envir = cachedEnv)$group,
            sequence_alignment_Dayhoff = get("sequence_alignment_Dayhoff",
                                             envir = cachedEnv)$group, 
            structure_alignments_Mirny = get("structure_alignments_Mirny", 
                                             envir = cachedEnv)$group,
            custom = get("custom", envir = cachedEnv)$group,
            no = dat
        )
        dat
    }
    
    bg <- lapply(bg, function(.bg) groupAA(.bg, groupingScheme))
    
    exp <- groupAA(exp, groupingScheme)
    if (ncol(exp) != ncol(bg[[1]]))
    {
        stop("The length of background is different from inputs", call. = FALSE)
    }
    
    ## get counts for each amino acid at each position if Fisher's exact test
    ## otherwise get proportions
    testType <- dagBackground@testType
    counts <- function(mat, coln, testType) 
    {
        num <- apply(mat, 2, function(.ele) {
            cnt <- table(.ele)[coln]
            ## just in case not all coln in the dataset
            names(cnt) <- coln 
            cnt[is.na(cnt)] <- 0
            total <- sum(cnt)
            list(freq = cnt, percent = cnt / total)})
    }
    ## a list of list
    bg <- lapply(bg, counts, coln, testType = testType)
 
    ## a list
    exp <- counts(exp, coln, testType = testType)
    ## get frequency of each amino acid at each position
    exp_freq <- do.call(cbind, lapply(exp, function(.ele){
       .ele$freq
    }))
    exp_freq[is.na(exp_freq)] <- 0
    ## get percentage of each amino acid at each position
    exp_percent <- do.call(cbind, lapply(exp, function(.ele){
        .ele$percent
    }))
    exp_percent[is.na(exp_percent)] <- 0
    rownames(exp_freq) <- rownames(exp_percent ) <- coln
    
    ## Add background noise to the background model
    if(testType == "z-test")
    {
        per_sample_bg_percent <- lapply(bg, function(.bg) {
                do.call(cbind, lapply(.bg, function(.ele){
                    .ele$percent
                }))
            })
        per_position_bg_percent <- lapply(1:ncol(exp_percent),function(i) {
            do.call(cbind, lapply(per_sample_bg_percent, function(.ele){
            .ele[, i]
            }))
        })
        
        if (!is.na(bgNoise)) 
        {
            rdirichlet <- function (n, alpha) {
                l <- length(alpha)
                x <-matrix(rgamma(l * n, alpha),
                           ncol = n,
                           byrow = TRUE)
                sm <- x %*% rep(1, n)
                t(x / as.vector(sm))
            }
            bg_percent <- lapply(per_position_bg_percent, function(col_pct) {
                (1 - bgNoise) * col_pct + bgNoise * 
                    rdirichlet(nrow(col_pct), rep(1, ncol(col_pct)))
            })
        }
        
        ## get mean proportion of each amino acid at each position
        mu_percent <- do.call(cbind, lapply(bg_percent, function(.bg) {
            apply(.bg, 1, mean, na.rm = TRUE)
        }))
        
        diff_percent <- exp_percent - mu_percent
        diff_percent[is.na(diff_percent)] <- 0
        n <-  dagBackground@numSubsamples
        std_percent <- mu_percent*(1-mu_percent)/n
        statistics <- diff_percent / std_percent
        pvalue <- 2 * pnorm(-abs(statistics))
    } else 
    {
        bg <- bg[[1]]
        bg_freq <- do.call(cbind, lapply(bg, function(.ele){
            .ele$freq
        }))
        mu_percent <- do.call(cbind, lapply(bg, function(.ele){
            .ele$percent
        }))

        diff_percent <- exp_percent - mu_percent
        
        ## Fisher exact test
        non_bg_freq <- -sweep(x= bg_freq, MARGIN= 2, STATS = colSums(bg_freq), FUN = "-") 
        non_exp_freq <- -sweep(x= exp_freq, MARGIN= 2, STATS = colSums(exp_freq), FUN = "-")
        
        testResults <- mapply(function(x, X, y, Y){
            fisher <- fisher.test(matrix(c(x, X, y, Y), nrow =2, byrow =TRUE),
                                  alternative = "two.sided")
            c(pvalues=fisher$p.value, statistics = fisher$estimate)
        }, exp_freq, non_exp_freq, bg_freq, non_bg_freq, SIMPLIFY = FALSE)
        
        testOut <- as.data.frame(do.call(rbind, testResults))
        colnames(testOut) <- c("pvalues", "statistics")
        pvalue <- matrix(testOut$pvalues, nrow = nrow(exp_percent), byrow =FALSE)
        statistics <- matrix(testOut$statistics, nrow = nrow(exp_percent), byrow =FALSE)
    }
    rownames(diff_percent) <- rownames(pvalue) <- rownames(statistics) <- coln
   
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
    if (length(coln) == ncol(diff_percent)) 
    {
        colnames(diff_percent) <- colnames(statistics) <- colnames(pvalue) <- coln
    } else
    {
        coln <- paste("AA", 1:ncol(diff_percent), sep = "")
        colnames(diff_percent) <- colnames(statistics) <- colnames(pvalue) <- coln
    }

    
    ## return the test results as an object of testDAUresults class
    new(
        "testDAUresults",
        group = groupingScheme,
        difference = diff_percent,
        statistics = statistics,
        pvalue = pvalue,
        background = mu_percent,
        motif = exp_percent,
        testType = dagBackground@testType,
        upstreamOffset = dagPeptides@upstreamOffset,
        downstreamOffset = dagPeptides@downstreamOffset
    )
}
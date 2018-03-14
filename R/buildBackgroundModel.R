#' Create a new \code{backgroundModel} object for \code{testDAU}.
#'
#' @param background A list each element of which is a vector of aligned peptide
#' sequences of the same length..
#' @param numSubsamples An integer, the number of random sampling to get 
#' background sequence sets.
#'
#' @return An object of \code{\link{dagBackground}} class.
#' @export
#' @keywords internal
#'
#' @examples
initiateBackgroundModel <- function(background, numSubsamples = 1L, testType)
{
    new("dagBackground",
        background = background,
        numSubsamples = numSubsamples,
        testType = testType)
}


#' Build a background model for Z-test.
#'
#' @param dagPeptides An object of \code{\link{dagPeptides}} class containing 
#' peptide sequences as the input set.
#' @param matches A character vector with the matched subsequences.
#' @param numSubsamples An integer, the number of random sampling.
#' @param rand.seed An integer, the seed used to perform random sampling
#' @param uniqueSeq A logical vector indicating whether only unique peptide sequences
#' are included in a background model.
#' @param replacement A logical vector of length 1, indicating whether replacement 
#' is allowed for random sampling.
#'
#' @return An object of \code{\link{dagBackground}} class.
#' @export
#' @keywords internal
#'
#' @examples
#' 
buildZTestBackgroundModel <- function(dagPeptides,
                                      matches,
                                      numSubsamples = 30L,
                                      rand.seed = 1,
                                      testType = testType,
                                      uniqueSeq = TRUE,
                                      replacement = FALSE)
{
    numSubsamples <- as.integer(numSubsamples)
    if (numSubsamples < 2)
    {
        stop("numSubsamples should be greater than 1", call. = FALSE)
    }
    
    #### sampling the same number of sequences from background as the experiment set
    set.seed(rand.seed)
    n <- nrow(dagPeptides@data)
    
    if (length(matches) < n)
    {
        stop("Too few matches in the background. Please try different parameters.",
             call. = FALSE)
    }
    
    ## subsampleing: samplesize = n; number of subsamples = numSubsamples
    background <- lapply(seq_len(numSubsamples), function(p) {
        s <- sample(matches, n, replace = replacement, prob = NULL)
        do.call(rbind, strsplit(s, "", fixed = TRUE))
    })
    
    initiateBackgroundModel(background = background, 
                            numSubsamples = numSubsamples, testType = testType)
}

#' Build background models for DAU tests
#' 
#' A method used to build background models for testing differential amino acid usage
#'
#' @param dagPeptides An object of \code{\link{dagPeptides}} class containing 
#' peptide sequences as the input set.
#' @param background A character vector with defaults: "wholeProteome" and "inputSet", 
#' "nonInputSet", indicating what set of peptide sequences should be considered to 
#' generate a background model.
#' @param model A character vector with defaults: "any" and "anchored", indicating 
#' whether an anchoring position should be applied to generate a background model.
#' @param targetPosition A character vector with defaults: "any", "Nterminus" and 
#' "Cterminus", indicating whether which part of protein sequences of choice should
#' be used to generate a background model.
#' @param uniqueSeq A logical vector indicating whether only unique peptide sequences
#' are included in a background model.
#' @param numSubsamples An integer, the number of random sampling.
#' @param rand.seed An integer, the seed used to perform random sampling
#' @param replacement A logical vector of length 1, indicating whether replacement 
#' is allowed for random sampling.
#' @param testType A character vector of length 1. Available options are "ztest" 
#' and "fisher".

#' @import 
#' @importFrom
#'
#' @return An object of \code{\link{dagBackground}} class.
#' @export
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' dat <- unlist(read.delim(system.file("extdata", "grB.txt", package="dagLogo"), 
#'                                      header=F, as.is=TRUE))
#' ##prepare proteome from a fasta file
#' proteome <- prepareProteome(fasta=system.file("extdata", 
#'                                              "HUMAN.fasta",
#'                                              package="dagLogo"))
#' ##prepare object of dagPeptides
#' seq <- formatSequence(seq=dat, proteome=proteome, 
#'                      upstreamOffset=14, downstreamOffset=15)
#' background <- buildBackgroundModel(seq, background="wholeProteome", proteome=proteome)

buildBackgroundModel <- function(dagPeptides,
                                 background = c("wholeProteome", "inputSet", "nonInputSet"),
                                 model = c("any", "anchored"),
                                 targetPosition = c("any", "Nterminus", "Cterminus"),
                                 uniqueSeq = TRUE,
                                 numSubsamples = 30L,
                                 rand.seed = 1,
                                 replacement = FALSE,
                                 testType = c("ztest", "fisher"),
                                 proteome) 
{
    if (missing(dagPeptides) || class(dagPeptides) != "dagPeptides") 
    {
        stop("dagPeptides should be an object of dagPeptides.", call. = FALSE)
    }
    background <- match.arg(background)
    targetPosition <- match.arg(targetPosition)
    model <- match.arg(model)
    testType <- match.arg(testType)
    numSubsamples <- as.integer(numSubsamples)
    if (numSubsamples < 2)
    {
        stop("numSubsamples should be greater than 1")
    }
    
    length <- dagPeptides@upstreamOffset + dagPeptides@downstreamOffset + 1
    
    ## Decide what set of peptide sequences should be used to build the background model
    if (background != "inputSet") 
    {
        if (missing(proteome) || class(proteome) != "Proteome") 
        {
            stop("proteome should be an object of the Proteome Class. \n
                Try ?prepareProteome to get help", call. = FALSE)
        }
        
        ## Whole proteome as background
        if (background == "wholeProteome") 
        {
            SequenceStr <- proteome@proteome$SEQUENCE
        } else
        {
            proteome.s <- proteome@proteome[!proteome@proteome$SEQUENCE %in% 
                                                dagPeptides@data$peptide, ]
            SequenceStr <- proteome.s$SEQUENCE
        }
    } else
    {
        SequenceStr <- dagPeptides@data$peptide
    }
    
    ## Apply the anchoring position restriction or not for extract subsequences for
    ## building background model. Different patterns are used indifferent scenario.
    if (model == "anchored") 
    {
        anchorAA <- table(dagPeptides@data$anchorAA)
        anchorAA <- anchorAA[order(anchorAA, decreasing = TRUE)]
        
        if (length(anchorAA) > 1) 
        {
            model <- "any"
            warning("Anchor amino acid is not unique. Model is set to 'any'.")
            anchorAA <- paste("[", paste(names(anchorAA), collapse = ""), "]",
                              sep = "")
        } else
        {
            anchorAA <- names(anchorAA)[1]
            pattern <-paste("([A-Z]{",dagPeptides@upstreamOffset,"}",
                            anchorAA,
                            "[A-Z]{",dagPeptides@downstreamOffset,"})",
                            sep = "")
        }
    }
    
    if (model == "any") 
    {
        length <- dagPeptides@upstreamOffset + dagPeptides@downstreamOffset + 1
        pattern <- paste("([A-Z]{", length , "})", sep = "")
    }
    
    ## Decide how to extract subsequences based on defined pattern: C-terminus only,
    ## N-terminus only or anywhere of peptide sequences. If C-terminus only or 
    ## N-terminus only is set, very restrictive extracting will be adopted. Isn't 
    ## these two method too restrictive? 
    if (targetPosition == "Cterminus")
    {
        pattern <- paste(pattern, "$", sep = "")
    } else if (targetPosition == "Nterminus")
    {
        pattern <- paste("^", pattern, sep = "")
    }
    matches <- gregexpr(pattern, SequenceStr)
    matches <- unlist(regmatches(SequenceStr, matches))
    if(uniqueSeq)
    {
        matches <- unique(matches)
    }
    
    ## build background model based on the type of hypothesis test
    if (testType == "fisher")
    {
        numSubsamples <- 1L
        background <- do.call(rbind, strsplit(matches, "", fixed = TRUE))
        backgroundModel <-
            initiateBackgroundModel(background = background,
                                    numSubsamples = numSubsamples,
                                    testType = testType)
        
    } else
    {
        backgroundModel <-
            buildZTestBackgroundModel(
                dagPeptides = dagPeptides,
                matches = matches,
                testType = testType,
                numSubsamples = as.integer(numSubsamples),
                rand.seed = rand.seed,
                uniqueSeq = uniqueSeq,
                replacement = replacement
            )
    }
    backgroundModel
}
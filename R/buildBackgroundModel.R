#' Build background models for DAU tests
#' 
#' A method used to build background models for testing differential amino acid usage
#'
#' @param dagPeptides An object of \link[dagPeptides]{dagPeptides} class containing 
#' peptide sequences as the input set.
#' @param bg A character vector with defaults: "wholeProteome" and "inputSet", 
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
#' @param proteome An object of \link[Proteome]{Proteome} class containing a whole
#' set of peptide sequences used for building background model using restriction of
#' choice.
#' @import 
#' @importFrom
#'
#' @return An object of \link[dagBackground]{dagBackground} class.
#' @export
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' 

buildBackgroundModel <- function(dagPeptides,
                                 bg = c("wholeProteome", "inputSet", "nonInputSet"),
                                 model = c("any", "anchored"),
                                 targetPosition = c("any", "Nterminus", "Cterminus"),
                                 uniqueSeq = TRUE,
                                 numSubsamples = 30L,
                                 rand.seed = 1,
                                 replacement = FALSE,
                                 proteome) 
{
    if (missing(dagPeptides) || class(dagPeptides) != "dagPeptides") 
    {
        stop("dagPeptides should be an object of dagPeptides.", call. = FALSE)
    }
    bg <- match.arg(bg)
    targetPosition <- match.arg(targetPosition)
    model <- match.arg(model)
    numSubsamples <- as.integer(numSubsamples)
    if (numSubsamples < 2)
    {
        stop("numSubsamples should be greater than 1")
    }

    ## Decide what set of peptide sequences should be used to build the background model
    if (bg != "inputSet") 
    {
        if (missing(proteome) || class(proteome) != "Proteome") 
        {
            stop("proteome should be an object of the Proteome Class. \n
                Try ?prepareProteome to get help", call. = FALSE)
        }
        
        if (bg == "wholeProteome") 
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
            anchorAA <- paste("[", paste(names(anchorAA), collapse = ""), "]", sep = "")
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
    set.seed(rand.seed)
    

    if (length(matches) <= n)
    {
        stop("too less matches in background. Please try different parameters.",
             call. = FALSE)
    }
    
    ##  subsampling size
    n <- nrow(dagPeptides@data)
    background <- lapply(seq_len(numSubsamples), function(p) {
        s <- sample(matches, n, replace = replacement, prob = NULL)
        if (uniqueSeq) { s <- unique(s)}
        do.call(rbind, strsplit(s, "", fixed = TRUE))
    })
    
    new("dagBackground",
        background = background,
        numSubsamples = numSubsamples)
    }
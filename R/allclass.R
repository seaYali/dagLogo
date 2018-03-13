#' Class dagPeptides.
#' An S4 class to represent formatted, aligned peptides for DAGlogo analysis.
#' 
#' @slot data A data frame with columns: IDs, anchorAA, anchorPos, peptide and anchor.
#' @slot peptides A matrix of character, each element is a single-character symbol
#' for a amino acid.
#' @slot upstreamOffset An integer in the interval (0, 20): the upstream offset
#' relative to the anchoring position.
#' @slot downstreamOffset An integer in the interval (0, 20): the downstream 
#' offset relative to the anchoring position.
#' @slot type A character vector of length 1. Available options :"Uniprot", and 
#' "fasta".
#' 
#' @name dagPeptides-class
#' @rdname dagPeptides-class
#' @exportClass dagPeptides
#' @import methods
#' 
#' @author Jianhong Ou


setClass(
    "dagPeptides",
    representation(
        data = "data.frame",
        peptides = "matrix",
        upstreamOffset = "numeric",
        downstreamOffset = "numeric",
        type = "character"
    ),
    validity = function(object) {
        re <- TRUE
        if (object@upstreamOffset < 0 ||
            object@downstreamOffset < 0)
            re <-
                "upstreamOffset and downstreamOffset should be a integer greater than 0"
        peptides <- as.character(object@peptides)
        peptides <-
            peptides[(!is.na(peptides)) & (peptides != "NA")]
        if (!all(1 == nchar(peptides)))
            re <-
            "peptides must be a matrix with one amino acide in each position"
        re
    }
)

#' Class Proteome.
#' 
#' An S4 class to represent a whole proteome for DAGlogo analysis.
#' 
#' @slot proteome A data frame.
#' @slot type A character vector of length 1. Available options :"Uniprot", and 
#' "fasta".
#' @slot species A character vector of length 1, a conventional Latin name for 
#' a species.
#'
#' @name Proteome-class
#' @rdname Proteome-class
#' @exportClass Proteome
#' 
#' @import methods
#' 
#' @author Jianhong Ou

setClass(
    "Proteome",
    representation(
        proteome = "data.frame",
        type = "character",
        species = "character"
    ),
    validity = function(object) {
        re <- TRUE
        if (!object@type %in% c("fasta", "UniProt"))
            re <- "type must be fasta or UniProt"
        if (object@type == "UniProt" &&
            is.null(object@proteome$ENTREZ_GENE))
            re <-
                "when type equals to UniProt, ENTREZ_GENE column is required for proteome"
        if (is.null(object@proteome$SEQUENCE))
            re <- "proteome sequence is required"
        re
    }
)

#' Class dagBackground.
#' 
#' An S4 class to represent a background of a formatted, aligned peptide for 
#' dagoLogo analysis.
#' 
#' @slot background A list of data frame. Within each n-by-1 dataframe is a the 
#' aligned pedtide of same length.
#' @slot numSubsamples An integer. That is the length of the \code{background} list
#'
#' @name dagBackground-class
#' @rdname dagBackground-class
#' @exportClass dagBackground
#' @import methods
#' 
#' @author Jianhong Ou 


setClass("dagBackground",
         representation(background = "list",
                        numSubsamples = "integer"))


#' Class testDAUresults.
#' 
#' An S4 class to represent a DAU statistical test result from DAGlogo analysis.
#' 
#' @slot group A character vector of length 1, the type of method for grouping 
#' amino acid.
#' @slot difference A numeric matrix consisting of differences of amino acid. 
#' proportions between the test set and the background set of aligned, formatted 
#' peptides at each position.
#' @slot zscore A numeric matrix consisting of Z-scores. This slot is \code{NULL}
#' for an Fisher's exact test result.
#' @slot pvalue A numeric matrix consisting of p-values.
#' @slot background A numeric matrix consisting mino acid proportions in the 
#' background set of aligned, formatted peptides at each position.
#' @slot motif A numeric matrix consisting of proportions for DAGLogo. 
#' @slot upstreamOffset A positive integer between (0, 20): the upstream offset
#' relative to the anchoring position.
#' @slot downstreamOffset A positive integer between (0, 20): the upstream offset
#' relative to the anchoring position.
#'
#' @name testDAUresults-class
#' @rdname testDAUresults-class
#' @exportClass testDAUresults
#' @import methods
#' 
#' @author Jianhong Ou, Hiabo Liu

setClass(
    "testDAUresults",
    representation(
        group = "character",
        difference = "matrix",
        test = "character",
        zscore = "matrix",
        pvalue = "matrix",
        background = "matrix",
        motif = "matrix",
        upstreamOffset = "numeric",
        downstreamOffset = "numeric"
    ),
    validity = function(object) {
        re <- TRUE
        if (object@upstreamOffset < 0 || object@downstreamOffset < 0)
        {
            re <-
                "upstreamOffset and downstreamOffset should be a integer greater than 0"
        }

        if(!test %in% c("fisher", "ztest"))
        {
            re <- 'The test type must be "fisher", "ztest"'
        } else if (test == "fisher")
        {
            if (!is.null(object@zscore) || ncol(object@difference) == 0 || 
                ncol(object@pvalue) == 0)
            {
                re <-
                    "slots zscore, difference and pvalue could not be empty" 
            }
            if (any(dim(object@pvalue) != dim(object@difference)))
            {
                re <-
                    "dim of slots zscore, difference and pvalue should be identical"
            }
        } else
        {
            
            if (ncol(object@zscore) == 0 ||
                ncol(object@difference) == 0 || ncol(object@pvalue) == 0)
            {
                re <-
                    "slots zscore, difference and pvalue could not be empty"
            }
            if (any(dim(object@pvalue) != dim(object@difference)) || 
                any(dim(object@pvalue) != dim(object@zscore)))
            {
                re <-
                    "dim of slots zscore, difference and pvalue should be identical"
            }
        }
        re
    }
)

#' Format already aligned peptide sequences.
#' 
#' Convert already aligned peptide sequences into an object of \code{dagPeptides} 
#' class.
#' 
#' @param seq A vector of aligned peptide sequences of the same length
#' @param proteome An object of \code{Proteome} Class.
#' @param upstreamOffset An integer in the interval (0, 20): the upstream offset
#' relative to the anchoring position.
#' @param downstreamOffset An integer in the interval (0, 20): the downstream 
#' offset relative to the anchoring position.
#'
#' @return An object of dagPeptides Class
#' @export
#' @author Jianhong Ou
#' @examples
#' dat <- unlist(read.delim(system.file("extdata", "grB.txt", package="dagLogo"), 
#'                                      header=F, as.is=TRUE))
#' head(dat)
#' ##prepare proteome from a fasta file
#' proteome <- prepareProteome(fasta=system.file("extdata", "HUMAN.fasta", 
#'                                                package="dagLogo"))
#' ##prepare object of dagPeptides
#' seq <- formatSequence(seq=dat, proteome=proteome, upstreamOffset=14, 
#'                      downstreamOffset=15)
#' 
#' 

formatSequence <-function(seq,
                          proteome,
                          upstreamOffset,
                          downstreamOffset) 
{
    if (missing(proteome) || class(proteome) != "Proteome") 
    {
        stop("proteome should be an object of Proteome Class. \n
             Try ?prepareProteome to get help",
             call. = FALSE)
    }
    if (missing(seq)) 
    {
        stop("seq is a required parameter.", call. = FALSE)
    }
    if (!is.character(seq)) 
    {
        seq <- as.character(seq)
    }
    
    ## the number of characters in each seq should be the same
    width <- unique(unlist(lapply(seq, nchar)))
    
    if (length(width) > 1) 
    {
        stop("seq must be characters of the same length", call. = FALSE)
    }
    
    ## determine upstreamOffset and downstreamOffset
    if (missing(upstreamOffset) && missing(downstreamOffset)) 
    {
        upstreamOffset <- floor(width / 2)
        downstreamOffset <- width - upstreamOffset - 1
    } else
    {
        if (missing(downstreamOffset)) 
        {
            downstreamOffset <- width - upstreamOffset - 1
        } else
        {
            upstreamOffset <- width - downstreamOffset - 1
        }
    }
    
    ## check validity of upstreamOffset and downstreamOffset
    if (upstreamOffset <0 || upstreamOffset > 20 || 
        downstreamOffset < 0 || downstreamOffset >20)
    {
        stop("upstreamOffset and downstreamOffset must be in the 
             interval (0, 20).", call. = FALSE)
    }
    
    ## retrieve anchorAA and anchorPos
    center <- upstreamOffset + 1
    anchorAA <- unlist(lapply(seq, function(.ele)
                                   substr(.ele, center, center)))
    
    ## find the IDs and the starting positions of peptide sequences in Proteome 
    ## which match the sequences in seq
    m <- do.call(rbind, lapply(seq, function(.seq) {
        ## a numeric vector: starting positions of the matches, -1 if no match
        .m <- regexpr(.seq, proteome@proteome$SEQUENCE) 
        ## index of sequence with true match
        .id <- which(.m != -1)
        ## index of first sequence with true match
        .id <- .id[1]
        ## starting position of the first match
        .pos <- .m[.id]
        ## return the index and starting position of the first match
        c(.id, .pos)
    }))
    
    ## starting positions of the first match in sequences in Proteome for each 
    ## aligned sequence from seq
    anchorPos <- m[, 2]
    
    ## IDs of sequences in Proteome, matching each aligned sequence in seq
    IDs <- proteome@proteome[m[, 1], "ID"]
    
    ## full sequences in Proteome, matching each aligned sequence in seq
    peptide <- proteome@proteome[m[, 1], "SEQUENCE"]
    dat <- data.frame(IDs = IDs,
                      anchorAA = anchorAA,
                      anchorPos = anchorPos,
                      peptide = peptide,
                      anchor = anchorAA)
    
    ## upstream and downstrean sequences and characters
    dat$upstream <- substr(seq, 1, upstreamOffset)
    dat$downstream <-
        substr(seq, upstreamOffset + 2, upstreamOffset + downstreamOffset + 1)
    seqchar.upstream <-
        do.call(rbind, lapply(dat$upstream, function(.seq) {
            .seq <- c(rep("NA", upstreamOffset),
                      unlist(lapply(1:nchar(.seq), function(i)
                          substr(.seq, i, i))))
            .seq <- .seq[(length(.seq) - upstreamOffset + 1):length(.seq)]
            .seq
        }))
    seqchar.downstream <-
        do.call(rbind, lapply(dat$downstream, function(.seq) {
            .seq <-
                c(unlist(lapply(1:nchar(.seq), function(i)
                    substr(.seq, i, i))),
                  rep("NA", downstreamOffset))
            .seq <- .seq[1:downstreamOffset]
            .seq
        }))
    seqchar <-cbind(seqchar.upstream,
                    as.character(dat$anchorAA),
                    seqchar.downstream)
    new(
        "dagPeptides",
        data = dat,
        peptides = seqchar,
        upstreamOffset = upstreamOffset,
        downstreamOffset = downstreamOffset,
        type = proteome@type
    )
}
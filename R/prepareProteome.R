#' Create an object of Proteome class.
#'
#' Create an object of Proteome class by downloading a whole proteome data from 
#' Uniprot for a given organism of an NCBI taxonomy ID or by using peptide sequences
#' in a fasta file.
#'
#' @param source A character vector of length 1 or NULL. A database source from 
#' where the proteome sequences are to downloaded. By default, currently it is
#' "Uniprot". If it is NULL, then \code{fastaFile} has to be specified. The 
#' priority of \code{source} is higher than \code{fastaFile}.
#' @param destDir A character vector of length 1. A destination directory with 
#' writing permission for saving downloaded sequences.
#' @param fastaFile A character vector of length 1. A fasta file name from which
#' to read in protein sequences.
#' @param ... other parameters passing to \link{download.file}.
#' @param species A character vector of length 1. The Latin name of a species 
#' confirming to the Linnaean taxonomy nomenclature system.
#'
#' @importFrom UniProt.ws availableUniprotSpecies
#' @import methods
#' @importFrom Biostrings readAAStringSet
#' @importFrom utils download.file
#' @return An object of Proteome
#' @export
#'
#' @examples 
#' ## Prepare an objecto of Proteome from the UniProt database
#' if(interactive()) 
#' {
#'     proteome <- prepareProteome(source = "UniProt", species = "Homo sapiens")
#' }
#' 
#' ## Prepare n objecto of Proteome from a fasta file
#' fasta <- system.file("extdata", "Drosophila.melanogaster.fa", package="dagLogo")
#' proteome <- prepareProteome(source = NULL, species = "Drosophila melanogaster", fastaFile=fasta)

prepareProteome <- function(source = NULL, species = NULL, 
                            destDir=getwd(), fastaFile, ...) 
{
    if (is.null(species))
    {
        stop("Parameter species must be specified!", call. = FALSE)
    }
    if (!is.null(source))
    {
        if (source != "UniProt")
        {
            stop("source must be UniProt.", call. = FALSE)
        }
        tempFile <- tempfile(pattern = gsub(" ", ".", paste0(source, species, sep =".")), 
                             tmpdir = destDir, fileext = ".fasta")
        
        speciesInfo <- availableUniprotSpecies(pattern = paste0("^", species, "$"))
        if (nrow(speciesInfo) != 1) 
        {
            stop("Parameter species is wrong. Please specify a correct Latin name 
                 for the species!", call. = FALSE)
        } else 
        {
            taxonID <- speciesInfo[1, 1]
        }
        
        ## UniProt url for downloading the proteome of a given species
        url <- paste0("http://www.uniprot.org/uniprot/?query=organism:", 
                      taxonID,"&format=fasta")
        
        download.file(url = url, destfile = tempFile, ...)
        fasta <- readAAStringSet(tempFile)
        proteome <- data.frame(
            SEQUENCE = as.character(fasta),
            ID = gsub("^>.+?\\|(.+?)\\|.+", ">\\1", names(fasta), perl = TRUE),
            stringsAsFactors = FALSE)
        return(new(
            "Proteome",
            proteome = proteome,
            type = "UniProt",
            species = species))
    } else if (!missing(fastaFile)) ## prepare Proteome from a fasta file
    {
        if (length(fastaFile) == 1 && class(fastaFile) == "character") 
        {
            fasta <- readAAStringSet(fastaFile)
        }
        if (class(fasta) != "AAStringSet") 
        {
            stop("fasta should be an object AAStringSet",
                 call. = FALSE)
        }
        proteome <- data.frame(
            SEQUENCE = as.character(fasta),
            
            ## ID may be simplified by using gsub
            ID = names(fasta),
            stringsAsFactors = FALSE)
        return(new(
            "Proteome",
            proteome = proteome,
            type = "fasta",
            species = species))
    } else
    {
        stop("At least one of arguments source and fastaFile should be provided!",
             call. = FALSE)
    }
}

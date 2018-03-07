
#' Create an object of Proteome class
#'
#' Create an object of Proteome class by downloading a whole proteome data from 
#' Uniprot for a given organism of an NCBI taxonomy ID or by using peptide sequences
#' in a fasta file.
#'
#' @param source A database source from where the proteome sequences are to 
#' downloaded. By default, currently it is "Uniprot". If it is NULL, then 
#' fastaFile has to be specified. The priority of source is higher than dataFile.
#' @param destDir A destination directory with writing permission for saving 
#' downloaded sequences.
#' @param fastaFile A fasta file name from which to read in protein sequences.
#' @param ... other parameters passing to \link{download.file}.
#' @param species The Latin name of a species confirming to the Linnaean taxonomy 
#' nomenclature system.
#'
#' @improt UniProt.ws
#' @importFrom utils download.file
#' @return An object of Proteome
#' @export
#'
#' @examples 
#' proteome <- prepareProteome(source = "Uniprot", species = "Drosophila melanogaster")
#' 
#' fasta <- system.file("extdata", "Drosophila.melanogaster.fa", package="dagLogo")
#' proteome <- prepareProteome(source = NULL, species = "Drosophila melanogaster", fastaFile=fasta)

prepareProteome <- function(source = "Uniprot", species = "unknown", destDir=getwd(), fastaFile, ...) 
{
    if (missing(species))
    {
        stop("Argument species is necessary!", call. = FALSE)
    }
    if (!is.null(source))
    {
        tempFile <- tempfile(pattern = gsub(" ", ".", paste0(source, species, sep =".")), 
                             tmpdir = destDir, fileext = ".fasta")
        
        speciesInfo <- availableUniprotSpecies(pattern = paste0("^", species, "$"))
        if (nrow(speciesInfo) != 1) 
        {
            stop("Argument species is wrong. Please specify a correct Latin name for the species!")
        } else 
        {
            taxonID <- speciesInfo[1, 1]
        }
        
        url <- paste0("http://www.uniprot.org/uniprot/?query=organism:", taxonID,"&format=fasta")
        download.file(url = url, destfile = tempFile, ...)
        fasta <- readAAStringSet(tempFile)
    } else if (!missing(fastaFile)) 
    {
        if (length(fastaFile) == 1 && class(fastaFile) == "character") 
        {
            fasta <- readAAStringSet(fastaFile)
        }
    } else
    {
        stop("At least one of arguments source and fastaFile should be provided!", call. = FALSE)
    }
    
    if (class(fasta) != "AAStringSet") 
    {
        stop("fasta should be an object AAStringSet",
             call. = FALSE)
    }
    proteome <- data.frame(
        SEQUENCE = as.character(fasta),
        ID = names(fasta),
        stringsAsFactors = FALSE)
    return(new(
        "Proteome",
        proteome = proteome,
        type = "fasta",
        species = species))
}

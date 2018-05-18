#' Create an object of \code{\link{Proteome}} Class.
#'
#' Create an object of \code{\link{Proteome}} class by downloading a whole 
#' proteome data from UniProt for a given organism of an NCBI taxonomy ID or 
#' species' scientific name, or by using peptide sequences in a fasta file.
#'
#' @param source A character vector of length 1 or NULL. A database source from 
#' which the proteome sequences are to be downloaded. By default, currently it is
#' "UniProt". If it is NULL, then \code{fastaFile} has to be specified. The 
#' priority of \code{source} is higher than \code{fastaFile}.
#' @param taxonID Taxonomy ID for a species of interest. Check the NCBI
#' taxonomy database: \link{https://www.ncbi.nlm.nih.gov/taxonomy} or 
#' the UniProt database \link{http://www.uniprot.org/taxonomy/}. At least one of 
#' the two parameters, \code{taxonID} and \code{species}, should be specified. 
#' If both are specified, \code{taxonID} will be used preferentially.
#' @param species A character vector of length 1. The Latin name of a species 
#' confirming to the Linnaean taxonomy nomenclature system. CAUTION: for species 
#' with different strains, attention should be paid. You can interactively choose
#' the right \code{taxonID} from an output list.
#' @param destDir A character vector of length 1. A destination directory with 
#' writing permission for saving downloaded sequences. Default is the current
#' working directory.
#' @param fastaFile A character vector of length 1. A fasta file name from which
#' to read in protein sequences.
#' @param ... other parameters passing to the function \code{\link{download.file}}.
#'
#' @importFrom UniProt.ws availableUniprotSpecies
#' @import methods
#' @importFrom Biostrings readAAStringSet
#' @importFrom utils download.file
#' @return An object of Proteome
#' @export
#'
#' @examples 
#' ## Prepare an objecto of Proteome Class from the UniProt database
#' proteome <- prepareProteome(source = "UniProt", species = "Homo sapiens")
#' 
#' ## Prepare an objecto of Proteome from a fasta file
#' fasta <- system.file("extdata", "Drosophila.melanogaster.fa", package="dagLogo")
#' proteome <- prepareProteome(source = NULL, species = "Drosophila melanogaster", 
#' fastaFile=fasta)

prepareProteome <- function(source = NULL, taxonID = NULL, species = NULL, 
                            destDir=getwd(), fastaFile, ...) 
{
    if (is.null(taxonID) && is.null(species))
    {
        stop("Parameter species or taxonID must be specified!", call. = FALSE)
    }
    if (!is.null(source))
    {
        if (!source %in% c("UniProt", "Uniprot", "uniprot"))
        {
            stop("Source must be UniProt.", call. = FALSE)
        } else
        {
            source <- "UniProt"
        }
        tempFile <- tempfile(pattern = gsub(" ", ".", paste0(source, species, sep =".")), 
                             tmpdir = destDir, fileext = ".fasta")
        if (!is.null(species) && is.null(taxonID))
        {
            speciesInfo <- availableUniprotSpecies(pattern = paste0("^", species))
            numberOfEntries <- nrow(speciesInfo)
            if ( numberOfEntries != 1) 
            {
                if(numberOfEntries == 0)
                {
                    stop("No matching species in the UniProt database. ",
                         "Parameter species might be wrong. ", 
                         "Please specify a correct Latin name for the species!",
                         call. = FALSE)
                } else 
                {
                    if (interactive())
                    {
                        print(speciesInfo) 
                        message("There is multiple matching entries.")
                        taxonID <- readline(prompt="Please select the right taxon ID based on the printout above: taxon ID: ")
                    } else 
                    {
                        stop("There are multiple entries in the Uniprot database", 
                             " matching your input species. ", 
                             "Please specify the exact one by referring to the ", 
                             "NCBI taxonomy database: https://www.ncbi.nlm.nih.gov/taxonomy ", 
                             "or the UniProt database http://www.uniprot.org/taxonomy.")
                    }
                }
            } else 
            {
                taxonID <- speciesInfo[1, 1]
            }
        } 
        ## check validity of taxon ID 
        species <- lookupUniprotSpeciesFromTaxId(taxonID)
        message(paste0("Downloading data for species: ", species))
        
        ## UniProt url for downloading the proteome of a given species
        url <- paste0("http://www.uniprot.org/uniprot/?query=organism:", 
                      taxonID,"&format=fasta")
        
        download.file(url = url, destfile = tempFile, ...)
        fasta <- readAAStringSet(tempFile)
        proteome <- data.frame(
            SEQUENCE = as.character(fasta),
            ID = gsub("^.+?\\|(.+?)\\|.+", "\\1", names(fasta), perl = TRUE),
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
        return(new("Proteome",
                   proteome = proteome,
                   type = "fasta",
                   species = species))
    } else
    {
        stop("At least one of arguments source and fastaFile should be provided!",
             call. = FALSE)
    }
}

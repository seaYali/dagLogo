#' An object of \code{Proteome} representing the E. coli proteome.
#'
#' A dataset containing the E. coli proteome. The variables are as follows:
#'
#' \itemize{
#'   \item proteome. A data frame containg ENTREZ_GENE IDs, protein SEQUENCES, 
#'                   and UniProtKB IDs and length of the sequences
#'   \item type. The data source: Uniprot
#'   \item species. Escherichia coli
#' }
#'
#' @format A object of Proteome class
#' @source \url{http://www.uniprot.org/}
"ecoli.proteome"

#' An object of \code{Proteome} representing the E. coli proteome.
#'
#' A dataset containing the E. coli proteome. The variables are as follows:
#'
#' \itemize{
#'   \item proteome. A data frame containg ENTREZ_GENE IDs, protein SEQUENCES, 
#'                   and UniProtKB IDs and length of the sequences
#'   \item type. The data source: Uniprot
#'   \item species. Drosophila melanogaster
#' }
#'
#' @format A object of Proteome class
#' @source \url{http://www.uniprot.org/}
"proteome.example"

#' An object of \code{dagPeptides} representing acetylated Lysine-containg peptides.
#'
#' A dataset containing the acetylated Lysine-containg peptides from Drosophila 
#' melanogaster. The variables are as follows:
#'
#' \itemize{
#'   \item data. A data frame containg entrezgene IDs, anchorAA, anchorPos,
#'               peptide, anchor, upstreamOffset and downstreamOffset.
#'   \item peptides. A matrix of individual symbols of aligned amino acids.
#'   \item upstreamOffset. 
#'   \item downstreamOffset.
#'   \item type. The ID type
#' }
#'
#' @format A object of \code{dagPeptides} class
"seq.example"
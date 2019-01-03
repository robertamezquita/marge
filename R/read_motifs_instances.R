#' Read Instances of Motifs Across Regions
#'
#' Companion function to \code{find_motifs_instances()} to read in and tidy
#' the results.
#' 
#' @param path File path pointing to results from \code{find_motifs_instances()}.
#'
#' @return A tibble with the following columns:
#' \itemize{
#'   \item \code{region_id} The identifier given to the regions provided
#'     for searching for a given motif. This region id is typically the
#'     fourth column (following the positional coordinates) in a bed file.
#'   \item \code{offset} The location of the motif with respect to the
#'     center of the region.
#'   \item \code{sequence} The genomic sequence at the specified location
#'     within the region.
#'   \item \code{motif_name} The readable name of the sequence.
#'   \item \code{strand} The strand on which the motif was found.
#'   \item \code{motif_score} The similarity score between the motif and
#'     sequence. A log odds score of the motif matrix, where higher scores
#'     are better matches.
#' }
#'
#' @export

read_motifs_instances <- function(path) {
    dat <- suppressMessages(readr::read_tsv(path))
    colnames(dat) <- c("region_id", "offset", "sequence",
                       "motif_name", "strand", "motif_score")
    return(dat)
}

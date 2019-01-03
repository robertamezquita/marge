#' HOMER Database of Motifs
#'
#' Dataset containing all known motifs included with HOMER by default.
#'
#' @format A tibble enumerating 3229 motifs from various motif sets
#'   (organisms), as described by:
#' \describe{
#'   \item{mset}{Motif set (in most cases, organism) from which motif
#'     was derived.}
#'   \item{consensus}{Consensus nucleotide sequence.}
#'   \item{motif_name}{Human readable motif name.}
#'   \item{motif_family}{Family which motif belongs to.}
#'   \item{experiment}{Experiment that identified the motif.}
#'   \item{accession}{Public identifier for raw/processed data from experiment.}
#'   \item{database}{Database where motif was defined from.}
#'   \item{log_odds_detection}{Log odds threshold to set for defining when a motif
#'     is enriched in a given sequence.}
#'   \item{motif_pwm}{List vector with motif position weight matrices (PWMs).}
#'   \item{log_p_value_detection}{Negative log p-value with which motif was
#'     detected (not used, but included for completeness).}
#' }
#'
#' @source \url{http://homer.ucsd.edu/homer/}
"HOMER_motifs"

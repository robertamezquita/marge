#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

. <- NULL

`%||%` <- function(a, b) if (is.null(a)) b else a

.calc_free_mem <- function() {
    ## Currently broken
    ##    system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE) %>%
    ##        as.numeric / 1000
    return(4000)
}

.parse_pcts <- function(x) {
    stringr::str_replace(x, "%", "") %>%
        as.numeric() * 0.01
}

## Add on motif_pwm to known_results data using HOMER_motifs.rda
.append_known_pwm <- function(known_results) {
    data(HOMER_motifs, envir = environment())
    hm <- HOMER_motifs %>%
##        dplyr::filter(rlang::UQ(rlang::sym('log_odds_detection')) > 0) %>%
        dplyr::select("motif_name", "motif_family", "experiment", "accession",
                      "motif_pwm", "log_odds_detection")
    dplyr::inner_join(known_results, hm, 
                      by = c("motif_name", "motif_family", "experiment", "accession"))
}    

## Calculate consensus sequence from a PWM tibble
.calc_consensus <- function(motif_pwm) {
    wm <- apply(motif_pwm, 1, which.max)
    consensus <- colnames(motif_pwm)[wm] %>%
        paste0(., collapse = '')
    return(consensus)
}

.write_bed <- function(x, path) {
    readr::write_tsv(x, path, col_names = FALSE)
}

.read_bed <- function(file) {
    readr::read_tsv(file, col_names = FALSE)
}

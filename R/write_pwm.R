#' Write Motif PWM in HOMER format
#'
#' Create a minimal HOMER formatted motif file based on PWM and required
#' motif info.
#'
#' @param motif_pwm tibble with alphabet columns
#' @param motif_name name of the sequence
#' @param consensus the accepted consensus motif sequence
#'   [optional: calculated from pwm if not provided]
#' @param log_odds_detection the threshold for detecting given motif
#' @param file path to new motif file; suffix should be '.motif'
#' @param append whether to overwrite existing file [default: TRUE, e.g.
#'   add on to existing file]
#'
#' @return nothing; called for its side-effect of producing HOMER motif file
#'
#' @export
#'
#' @examples
#' \dontrun{
#' knw_pwm <- read_known_results("inst/extdata/knownResults.txt", homer_dir = FALSE)
#'
#' ## Write multiple motifs
#' pwalk(list(motif_pwm = knw_pwm$motif_pwm,
#'            motif_name = knw_pwm$motif_name,
#'            log_odds_detection = knw_pwm$log_odds_detection),
#'       write_homer_motif, file = "test.motif", append = TRUE)
#' }

write_homer_motif <- function(motif_pwm, motif_name, log_odds_detection,
                            consensus = NULL,
                            file, append = TRUE) {
    ## Calculate consensus sequence if not provided
    if (is.null(consensus)) {
        consensus <- .calc_consensus(motif_pwm)
    }

    ## Create HOMER style header
    header <- paste(
        paste0('>', consensus),
        motif_name,
        log_odds_detection,
        sep = '\t'
    )

    ## Process PWM to printable form
    motif_rows <- apply(motif_pwm, 1, function(x) {
        paste0(paste(x, collapse = '\t'), '\n')
    })

    ## Write to file
    cat(header, file = file, append = append)
    cat('\n', file = file, append = TRUE)
    cat(motif_rows, file = file, append = TRUE)
}


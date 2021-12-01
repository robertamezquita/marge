#' Read Known Enriched Motifs HOMER output
#'
#' Reads in results from a known motif enrichment analysis created by
#' HOMER.
#'
#' Following an analysis using \code{\link{find_motifs_genome}}, a HOMER directory
#' is created which analyses the enrichment of known motifs. This
#' function reads and parses in that file into a tidy format.
#'
#' @param path path to the HOMER directory where all outputs are
#'   stored
#' @param homer_dir does the path point to a HOMER directory;
#'   if \code{FALSE}, path must point to the file directly
#'   [default: TRUE]
#'
#' @return a tibble with the following columns:
#' \itemize{
#'   \item \code{motif_name} the readable name of the sequence
#'   \item \code{motif_family} of transcription factors the motif belongs to
#'   \item \code{experiment} from where the motif was identified
#'   \item \code{accession} source publication or online repository ID
#'   \item \code{database} the curator of the motif
#'   \item \code{consensus} the accepted consensus motif sequence
#'   \item \code{log_p_value} -log10(p-value) significance of enrichment
#'   \item \code{fdr} Benjamini-Hochberg corrected p-value
#'   \item \code{tgt_num} number of times motif appears in target sequences
#'   \item \code{tgt_pct} percent of times motif appears in target sequences
#'   \item \code{bgd_num} number of times motif appears in background sequences
#'   \item \code{bgd_pct} percent of times motif appears in background sequences
#' }
#' 
#' @importFrom stringr str_replace
#' @importFrom lazyeval interp
#' @importFrom dplyr vars contains
#' 
#' @export
#' @seealso \code{\link{read_denovo_results}}, \code{\link{find_motifs_genome}}
#' 

read_known_results <- function(path, homer_dir = TRUE) {
    if (homer_dir) {
        path <- paste0(path, "/knownResults.txt")
    } 
    if (!file.exists(path)) {
        warning(paste("File", path, "does not exist"))
        return(NULL)
    }

    ## Read in raw file
    col_spec <- readr::cols('c', 'c', 'd', '-', 'd', 'd', 'c', 'd', 'c')
    raw <- readr::read_tsv(path, col_types = col_spec)
    colnames(raw) <- c('motif_name', 'consensus',
                       'log_p_value', 'fdr',
                       'tgt_num', 'tgt_pct', 'bgd_num', 'bgd_pct')

    ## Parse down all the combined columns
    tmp <- tidyr::separate(data = raw, 
                           col = motif_name, 
                           into = c('motif_name', 'experiment', 'database'),
                           sep = '/', 
                           extra = 'drop') %>%
        dplyr::mutate(log_p_value = "-log10(log_p_value)")
    parsed <- .parse_homer_subfields(tmp) %>%
        dplyr::mutate(dplyr::across(contains('pct'), .parse_pcts))

    ## Add on motif_pwm from the HOMERdb (see utils.R)
    known <- .append_known_pwm(parsed)
    names(known$motif_pwm) <- known$motif_name
    return(known)
}


#' Read Denovo Enriched Motifs HOMER Output
#'
#' Reads in results from a denovo motif enrichment analysis created by
#' HOMER.
#' 
#' Following an analysis using \code{\link{find_motifs_genome}}, a HOMER directory
#' is created which analyses the enrichment of denovo motifs. This
#' function reads and parses in that file into a tidy format.
#'
#' @param path path to the HOMER directory where all outputs are
#'   stored
#' @param homer_dir does the path point to a HOMER directory;
#'   if \code{FALSE}, path must point to the file directly
#'   [default: TRUE]
#'
#' @return a tibble with the following columns:
#' \itemize{
#'   \item \code{consensus} the consensus sequence of the denovo motif
#'   \item \code{motif_name} name of the motif
#'   \item \code{log_odds_detection} threshold used to determine bound vs. unbound sites
#'   \item \code{motif_pwm} a list column with PWMs for each motif
#'   \item \code{log_p_value_detection} from the original experiment used to ID motif
#'   \item \code{tgt_num} number of times motif appears in target sequences
#'   \item \code{tgt_pct} percent of times motif appears in target sequences
#'   \item \code{bgd_num} number of times motif appears in background sequences
#'   \item \code{bgd_pct} percent of times motif appears in background sequences
#'   \item \code{log_p_value} final enrichment from experiment -log10(p-value)
#'   \item \code{tgt_pos} average position of motif in target sequences, where
#'     0 = start of sequences
#'   \item \code{tgt_std} standard deviation of position in target sequences
#'   \item \code{bgd_pos} average position of motif in background sequences,
#'     where 0 = start of sequences
#'   \item \code{bgd_std} standard deviation of position in background sequences
#'   \item \code{strand_bias} log ratio of + strand occurrences to - strand occurrences
#'   \item \code{multiplicity} average number of occurrences per sequence in
#'     sequences with 1 or more binding sites
#' }
#'
#' @export
#' @seealso \code{\link{read_known_results}}, \code{\link{find_motifs_genome}}

read_denovo_results <- function(path, homer_dir = TRUE) {
    if (homer_dir) {
        path <- paste0(path, "/homerMotifs.all.motifs")
    }
    if (!file.exists(path)) {
        warning(paste("File", path, "does not exist"))
        return(NULL)
    }
    
    ## Suppress too few values warning when reading FDR column with NA values
    motifs <- suppressWarnings(read_motif(path))
    names(motifs$motif_pwm) <- motifs$motif_name
    return(motifs)
}


#' Read HOMER Motif Files
#'
#' The core function to read in motif files, whether from the HOMER database,
#' from HOMER denovo motif enrichment results, or even custom motifs. In all
#' cases, these files must be in the HOMER-format. See below for more details.
#'
#' To read-in a HOMER-formatted motif, at a minimum, the first three
#' fields are required to properly ID the motif:
#' 
#' \itemize{
#'   \item \code{">" + Consensus sequence} The dominant or likeliest sequence
#'   \item \code{Motif name} Should be unique
#'   \item \code{Log odds detection threshold} determines bound vs. unbound sites
#' }
#'
#' The remaining extra fields of HOMER-formatted motifs are described at the
#' URL below, and primarily meant for interpreting motifs from HOMER's
#' own database. To read more about the HOMER format, see:
#' \url{http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html}
#'
#'
#' Note that HOMER also has additional information in the motif name
#' regarding its origin and identity. See the internal function \code{.parse_homer_subfields}
#' for more info and to break this field up.
#'
#' Subsequent lines (after the ">") describe the position weight matrix (PWM),
#' with columns in order of A, C, G, T describing the probabilities of
#' per position of each nucleotide.
#'
#' Note that it is possible to combine complete information (HOMER-formatted) motifs
#' with minimal motifs. Simply use \code{dplyr::bind_rows} for easy concatenation
#' despite column spec differences.
#'
#' @param path location of motif file
#' 
#' @return at minimum, a tibble with the following columns:
#' \itemize{
#'   \item \code{consensus} the consensus sequence of the denovo motif
#'   \item \code{motif_name} name of the motif
#'   \item \code{log_odds_detection} threshold used to determine bound vs. unbound sites
#'   \item \code{motif_pwm} a list column with PWMs for each motif
#' }
#'
#' The following columns are presented when available from complete \code{*.motif*} files
#' or from HOMER results directories:
#' 
#' \itemize{
#'   \item \code{log_p_value_detection} from the original experiment used to ID motif
#'   \item \code{tgt_num} number of times motif appears in target sequences
#'   \item \code{tgt_pct} percent of times motif appears in target sequences
#'   \item \code{bgd_num} number of times motif appears in background sequences
#'   \item \code{bgd_pct} percent of times motif appears in background sequences
#'   \item \code{log_p_value} final enrichment from experiment -log10(p-value)
#'   \item \code{tgt_pos} average position of motif in target sequences, where
#'     0 = start of sequences
#'   \item \code{tgt_std} standard deviation of position in target sequences
#'   \item \code{bgd_pos} average position of motif in background sequences,
#'     where 0 = start of sequences
#'   \item \code{bgd_std} standard deviation of position in background sequences
#'   \item \code{strand_bias} log ratio of + strand occurrences to - strand occurrences
#'   \item \code{multiplicity} average number of occurrences per sequence in
#'     sequences with 1 or more binding sites
#' }
#' 
#' @importFrom stringr str_replace str_split str_detect
#' @importFrom lazyeval interp
#' @importFrom dplyr vars contains select_ rename_ group_by group_by_ bind_cols
#' @importFrom dplyr mutate_ mutate_at mutate_all mutate 
#' @importFrom tidyr separate separate_ nest
#' @importFrom purrr map
#' @importFrom readr read_tsv
#' @export

## Note: Field 'X6' is not read-in from using denovo method on
## Homer database files; 6th field is not described in HOMER docs

read_motif <- function(path) {
    ## Read in raw file
    all <- suppressMessages(suppressWarnings(readr::read_tsv(path, col_names = FALSE)))

    ## Calculate separation of each motif using '>' for extraction
    gt <- stringr::str_detect(all$X1, '>')

    ## --------------------------------------------------------
    ## Construct base motif info
    motif_info <- all[gt, 1:3] %>%
      dplyr::rename(consensus = 'X1', motif_name = 'X2', log_odds_detection = 'X3') %>%
      dplyr::mutate(consensus = gsub("^>", "", consensus))
    
    ## Munge motif PWMs
    ## Get motifs - where no '>' is detected
    ## Each motif instance is separated by a '>'
    ## Use cumsum to add up where '>' occur consecutively
    motif_pwm <- all[!gt, c('X1', 'X2', 'X3', 'X4')] %>%
      dplyr::mutate(motif_id = cumsum(gt)[!gt]) %>%
      dplyr::rename(A = 'X1', C = 'X2', G = 'X3', T = 'X4') %>%
      dplyr::group_by(motif_id) %>%
      dplyr::mutate(across(all_of(c('A', 'C', 'G', 'T')), as.numeric)) %>%
      tidyr::nest(motif_pwm = c('A', 'C', 'G', 'T'))

    ## Combine PWM + info
    motif_info <- dplyr::bind_cols(motif_info, motif_pwm)

    ## Try to parse subfields based on detecting splits
    ## Else it returns the same table
    motif_info <- .parse_homer_subfields(motif_info)
    
    ## For motif files with extra info
    if (ncol(all[gt, ]) > 3) {
        motif_info <- dplyr::mutate(.data = motif_info, log_p_value_detection = all[gt, ]$X4)
    }
    ## Return early if incomplete (old?) HOMER motif
    ## or custom with less than 7 columns
    if (ncol(all[gt, ]) < 7) {
        return(motif_info)
    } 

    ## --------------------------------------------------------
    ## Munge occurence and stats for complete info HOMER motifs
    motif_info <- dplyr::mutate(.data = motif_info, 
                                occurrence = all[gt, ]$X6, 
                                stats = all[gt, ]$X7)
    
    ## Helper formatting functions
    .drop_prior <- function(x, pattern = ":") {
      stringr::str_split(x, pattern) %>%
            purrr::map(function(x) { x[2] }) %>% unlist
    }
    .format_pct <- function(x) {
      stringr::str_replace(x, '%\\)', '') %>%
            as.numeric * 0.01
    }

    ## Munge motif position statistics
    stats <- motif_info[, 'stats'] %>%
        tidyr::separate(col = stats,
                 into = c('tgt_pos', 'tgt_std', 'bgd_pos', 
                          'bgd_std', 'strand_bias', 'multiplicity'), 
                 sep = ',') %>%
        dplyr::mutate(across(everything(), .drop_prior)) %>%
        dplyr::mutate(across(everything(), as.numeric))

    ## Munge motif occurrence metrics
    occurrence <- motif_info[, 'occurrence'] %>%
      tidyr::separate(col = occurrence, 
                      into = c('tgt_num', 'bgd_num', 'log_p_value', 'fdr'), 
                      sep = ',', 
                      fill = "right") %>%
      dplyr::mutate(across(everything(), .drop_prior)) %>% 
      tidyr::separate(col = tgt_num, into = c('tgt_num', 'tgt_pct'), sep = '\\(') %>%
      tidyr::separate(col = bgd_num, into = c('bgd_num', 'bgd_pct'), sep = '\\(') %>%
      dplyr::mutate(across(contains('_pct'), .format_pct)) %>%
      dplyr::mutate(across(everything(), as.numeric)) %>%
      dplyr::rename("-log10(log_p_value)" = "log_p_value")

    ## Put it all together
    motif_final <- motif_info %>%
      dplyr::select(-occurrence, -stats) %>%
      bind_cols(occurrence, stats)

    return(motif_final)
}


#' Parse HOMER subfields
#'
#' An internal function that parses the \code{motif_name} and
#' \code{experiment} subfields present in HOMER-formatted motifs.
#'
#' HOMER-formatted motifs contain three key fields:
#'
#' \itemize{
#'   \item \code{motif_name} human-readable name
#'   \item \code{experiment} how the motif was identified
#'   \item \code{database} what database the motif is from
#' }
#'
#' Furthermore, each motif may contain additional information enclosed
#' in parentheses. This includes:
#'
#' \itemize{
#'   \item \code{motif_family} included in the \code{motif_name} field,
#'     describes the transcription factor family.
#'   \item \code{accession} included in the \code{experiment} field,
#'     provides the specific ID or publication of the raw/processed data.
#' }
#'
#' @param motif_tbl a motif_tbl constructed while running \code{read_*_results}
#'   and consequently \code{read_motif}, or a custom tibble that
#'   has the above field separators as described with columns
#'   \code{motif_name} and \code{experiment}
#' 
#' @return Adds the above additional columns to an existing motif
#'   tibble returned by \code{read_motif} or \code{read_denovo_results}
#'
#' @seealso \code{\link{read_motif}}

.parse_homer_subfields <- function(motif_tbl) {
    cond <- stringr::str_detect(motif_tbl$motif_name, "/") %>%
        sum(., na.rm = TRUE) > 0
    if (cond) {
        motif_tbl <- tidyr::separate(data = motif_tbl, col = motif_name,
                                     into = c('motif_name', 'experiment', 'database'), 
                                     sep = '/', extra = "drop", fill = "right")
    }

    ## Detect if parentheses are present in motif_name
    ## to break apart into motif_name vs. motif_family
    cond <- stringr::str_detect(motif_tbl$motif_name, '\\(') %>%
        sum(., na.rm = TRUE) > 0
    if (cond) {
        motif_tbl <- tidyr::separate(data = motif_tbl, 
                                     col = motif_name, 
                                     into = c('motif_name', 'motif_family'), 
                                     sep = '\\(', extra = "drop", fill = "right")
        motif_tbl$motif_family <- stringr::str_replace(motif_tbl$motif_family, '\\)', '')
    }

    ## Detect If parentheses are present in experiment
    ## to break apart into experiment vs. accession
    if ("experiment" %in% colnames(motif_tbl)) {
        cond <- stringr::str_detect(motif_tbl$experiment, '\\(') %>%
            sum(., na.rm = TRUE) > 0
        if (cond) {
            motif_tbl <- tidyr::separate(data = motif_tbl, 
                                         col = experiment, 
                                         into = c('experiment', 'accession'), 
                                         sep = '\\(', extra = "drop", fill = "right")
            motif_tbl$accession <- stringr::str_replace(motif_tbl$accession, '\\)', '')
        }
    }
    return(motif_tbl)
}
#' Read de novo enriched motifs from HOMER .html output
#'
#' Reads in results from a .html file produced by HOMER de novo motif enrichment analysis.
#'
#' Following an analysis using \code{\link{find_motifs_genome}}, a HOMER directory
#' is created which analyses the enrichment of known motifs. This
#' function reads and parses in the .html file into a tidy format.
#' 
#' It might be more useful to read in the .html file instead of the .motifs file because
#' the .html contains information on known similar motifs.
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

read_denovo_html_results = function(path, homer_dir = TRUE) {
  if (homer_dir == TRUE) {
    path = paste0(path, "/homerResults.html")
  } 
  if (!file.exists(path)) {
    warning(paste("File", path, "does not exist"))
    return(NULL)
  }
  
  ## Read in html file
  html = readLines(path)
  
  testReturn = html
  return(testReturn)
}
#' Read de novo enriched motifs from HOMER .info.html output files
#'
#' Reads in results from all individual info.html files produced by HOMER de novo motif enrichment analysis.
#'
#' Following an analysis using \code{\link{find_motifs_genome}}, a HOMER directory
#' is created which analyses the enrichment of known motifs. This
#' function reads and parses in the .html file into a tidy format.
#'
#' It might be more useful to read in the .html files instead of the .motifs file because
#' the .html contain information on known similar motifs.
#'
#' @param path path to the HOMER directory where all outputs are
#'   stored
#' @param homer_dir does the path point to a HOMER directory;
#'   if \code{FALSE}, path must point to the file directly
#'   [default: TRUE]
#'
#' @return a list with the following structure:
#' \itemize{
#'   \item \code{Motif_information} a dataframe with most of the info returned by \code{\link{read_denovo_results}}
#'   \item \code{Matches_to_known_motifs} Dataframe with information on related known motifs matched to the de novo motif
#' }
#' The \code{Motif_information} dataframe has the following columns:
#' \itemize{
#'   \item \code{motif_name} name of the motif
#'   \item \code{consensus} the consensus sequence of the denovo motif
#'   \item \code{p_value} final enrichment from experiment (p-value)
#'   \item \code{log_p_value} final enrichment from experiment -log10(p-value)
#'   \item \code{info_content} information Content per bp
#'   \item \code{tgt_num} number of times motif appears in target sequences
#'   \item \code{tgt_pct} percent of times motif appears in target sequences
#'   \item \code{bgd_num} number of times motif appears in background sequences
#'   \item \code{bgd_pct} percent of times motif appears in background sequences
#'   \item \code{tgt_pos} average position of motif in target sequences, where 0 = start of sequences
#'   \item \code{bgd_pos} average position of motif in background sequences
#'   \item \code{strand_bias} strand Bias (log2 ratio + to - strand density)
#'   \item \code{multiplicity} multiplicity (# of sites on avg that occur together)
#' }
#' #' The \code{Matches_to_known_motifs} dataframe has the following columns:
#' \itemize{
#'   \item \code{motif_name} name of the motif
#'   \item \code{motif_family} family of the motif
#'   \item \code{ID} ID of the motif or the experiment
#'   \item \code{database} database of the motif
#'   \item \code{rank} match rank
#'   \item \code{score} log odds score of the motif matrix, higher scores are better matches
#'   \item \code{offset} from the center of the region
#'   \item \code{orientation} forward or reverse (on average?)
#'   \item \code{original_alignment} consensus sequence for comparison
#'   \item \code{matched_alignment} matched sequence with mismatches in "-" from the consensus. Can be offset.
#' }
#'#' @importFrom stringr str_replace str_split str_detect
#' @importFrom tidyr separate separate_
#' @export
# 
# library(stringr)
# library(tidyr)
# ### test params
# ###
read_denovo_html_results = function(path, homer_dir = TRUE) {
  if (homer_dir) {
    path = paste0(path, "/homerResults")
  }
  if (!file.exists(path)) {
    warning(paste("No files found"))
    return(NULL)
  }
  
  # Match correct html files
  filenames = list.files(path, pattern = "*.info.html")
  df_list = list()
  
  for (f in filenames) {
    print(f)
    ## Read in  html file
    html = readLines(paste(path, f, sep = "/"))
    
    # get number of motif from file
    mypattern = "motif([^<]*).info.html"
    n = gsub(mypattern, '\\1', grep(mypattern, f, value = TRUE))
    
    # Create dataframe
    df = data.frame(matrix(ncol = 13, nrow = 1))
    colnames(df) =   c(
      'motif_name',
      'consensus',
      'p_value',
      'log_p_value',
      "info_content",
      'tgt_num',
      'tgt_pct',
      'bgd_num',
      'bgd_pct',
      'tgt_pos',
      'bgd_pos',
      'strand_bias',
      'multiplicity'
    )
    # Main header
    mypattern = '<H2>Information for ([^<]*)</H2>'
    df$motif_name = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = paste('.*-([^<]*) \\(Motif ', n, '\\)</H2>', sep = "")
    df$consensus = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    # Other characteristics: p-value, log-pvalue etc as in main function
    mypattern = '<TR><TD>p-value:</TD><TD>([^<]*)</TD></TR>'
    df$p_value = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>log p-value:</TD><TD>([^<]*)</TD></TR>'
    df$log_p_value = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Information Content per bp:</TD><TD>([^<]*)</TD></TR>'
    df$info_content = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Number of Target Sequences with motif</TD><TD>([^<]*)</TD></TR>'
    df$tgt_num = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Percentage of Target Sequences with motif</TD><TD>([^<]*)</TD></TR>'
    df$tgt_pct = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Number of Background Sequences with motif</TD><TD>([^<]*)</TD></TR>'
    df$bgd_num = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Percentage of Background Sequences with motif</TD><TD>([^<]*)</TD></TR>'
    df$bgd_pct = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Average Position of motif in Targets</TD><TD>([^<]*)</TD></TR>'
    df$tgt_pos = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Average Position of motif in Background</TD><TD>([^<]*)</TD></TR>'
    df$bgd_pos = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Strand Bias \\(log2 ratio \\+ to \\- strand density\\)</TD><TD>([^<]*)</TD></TR>'
    df$strand_bias = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    mypattern = '<TR><TD>Multiplicity \\(# of sites on avg that occur together\\)</TD><TD>([^<]*)</TD></TR>'
    df$multiplicity = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
  
    df_list[[n]][["Motif_information"]] = df
    
    ########### new information
    mypattern = '<H4>([^<]*)</H4>'
    length_df = length(gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE)))

    df = data.frame(matrix(ncol = 9, nrow = length_df))
    colnames(df) =   c(
      'motif_name',
      'ID',
      'database',
      "rank",
      'score',
      'offset',
      'orientation',
      "original_alignment",
      "matched_alignment"
    )
    
    # Known motif matches
    mypattern = '<H4>([^<]*)</H4>'
    df$motif_name = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    df <- tidyr::separate(data = df, 
                          col = motif_name, 
                          into = c('motif_name', 'ID', 'database'), 
                          sep = '/',
                          extra = "drop",
                          fill = "right")
    
    ## Detect if parentheses are present in motif_name
    ## to break apart into motif_name vs. motif_family
    cond <- stringr::str_detect(df$motif_name, '\\(') %>%
      sum(., na.rm = TRUE) > 0
    if (cond) {
      df <- tidyr::separate(data = df, 
                            col = motif_name, 
                            into = c('motif_name', 'motif_family'),
                            sep = '\\(', 
                            extra = "drop", 
                            fill = "right")
      df$motif_family <- stringr::str_replace(df$motif_family, '\\)', '')
    }
    
    # Ranks
    mypattern = '<TR><TD>Match Rank:</TD><TD>([^<]*)</TD></TR>'
    df$rank = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    # Scores
    mypattern = '<TR><TD>Score:</TD><TD>([^<]*)</TD</TR>'
    df$score = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    # Offset
    mypattern = '<TR><TD>Offset:</TD><TD>([^<]*)</TD</TR>'
    df$offset = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    # Orientation
    mypattern = '<TR><TD>Orientation:</TD><TD>([^<]*)</TD></TR>'
    df$orientation = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    # Alignment
    # original
    mypattern = paste('.*-([^<]*) \\(Motif ', n, '\\)</H2>', sep = "")
    df$original_alignment = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    # matched
    mypattern = '.+>([^<]+)</FONT></TD></TR></TABLE>'
    df$matched_alignment = gsub(mypattern, '\\1', grep(mypattern, html, value = TRUE))
    
    df_list[[n]][["Matches_to_known_motifs"]] = df
  }
  return(df_list)
}
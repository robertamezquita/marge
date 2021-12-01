#' Find Motifs over Regions
#'
#' Calls \code{findMotifsGenome.pl} to run motif analysis over a given
#' set of regions.
#'
#' \code{find_motifs_genome} runs the core HOMER motif enrichment function
#' from the R system, and in the process generates (as a side-effect)
#' a HOMER results directory.
#'
#' This results directory is inspectable via a file browser, and contains
#' a summary of the results as HTML files as well as text files.
#'
#' For our purposes, within the directory two key files exist:
#'
#' \itemize{
#'   \item \code{knownResults.txt} known motifs that are enriched
#'   \item \code{homerResults.all.motif} denovo motifs that are enriched
#' }
#'
#' These two text files are the core results (all else can be discarded
#' by setting \code{keep_minimal} to \code{TRUE}, and are parsed downstream
#' by \code{\link{read_known_results}} and \code{\link{read_denovo_results}}.
#'
#' @param x \code{data.frame} with the first three columns being
#'   chromosome, start, and end coordinates, with a fourth column corresponding
#'   to a region identifier; extra columns may be kept; x may alternately be a
#'   path to an existing bed file of this format
#' @param path where to write HOMER results
#' @param genome ID of installed genome; check installed genomes using
#'   \code{list_homer_packages()}; examples include "hg38" and "mm10";
#'   add an 'r' at the end to mask repeats, e.g. "mm10r"
#' @param motif_length vector of motif lengths to consider [default is
#'   \code{c(8, 10, 12)}]
#' @param scan_size size of sequence to scan; this can be a numeric to
#'   specify the number of bases to scan centered on the region, or
#'   alternately can be set to "given" to scan the entire region;
#'   if using "given", will use the "-chopify" option to cut large
#'   background sequences to average of target sequence size
#'   [default: \code{100}]
#' @param optimize_count number of motifs to optimize [default: 8]
#' @param background \code{data.frame} containing coordinates of desired
#'   regions to use as the background; alternately may be a path to an
#'   existing bedfile; the default, "automatic", creates a background based
#'   on the GC content and (scan) size of the target sequences
#' @param local_background a numeric scalar specifying number of equal size
#'   regions around peaks to use as the local background
#'   [by default this is not used, e.g. default: \code{FALSE}]
#' @param only_known turns off searching for denovo motifs
#'   [default: \code{FALSE}]
#' @param only_denovo turns off searching for known motif enrichment
#'   [default: \code{FALSE}]
#' @param fdr_num number of randomizations to perform to calculate FDR 
#'   [default: 0]
#' @param cores number of cores to use 
#'   [default: \code{max(1, parallel::detectCores() - 2)}]
#' @param cache number in MB to use as cache to store sequences in memory
#'   [default: 8000Mb (8Gb)]
#' @param overwrite overwrite an existing HOMER results directory
#'   [default: \code{FALSE}]
#' @param keep_minimal remove all extra clutter from results, keep only
#'   the essentials (\code{knownResults.txt} and \code{homerMotifs.all.motifs}
#'   [default: \code{FALSE}]
#' @param scale_logos whether to scale sequence logos by information content 
#'   [default: \code{FALSE}]
#' @param showHomerStdErr Whether to print HOMER stderr to the R console 
#'   [default: \code{FALSE}]
#'   
#' @return Nothing; called for its side-effect of producing HOMER results
#' 
#' @export
#' @seealso \code{\link{read_known_results}}, \code{\link{read_denovo_results}}

find_motifs_genome <- function(x, path, genome,
                               motif_length = c(8, 10, 12),
                               scan_size = 100,
                               optimize_count = 8,
                               background = "automatic",
                               local_background = FALSE,
                               only_known = FALSE,
                               only_denovo = FALSE,
                               fdr_num = 0,
                               cores = unlist(options("mc.cores")) %||% max(1, parallel::detectCores() - 2),
                               cache = 8000, # syze in MB
                               overwrite = FALSE,
                               keep_minimal = FALSE,
                               scale_logos = FALSE, 
                               showHomerStdErr = FALSE) {

    ## Error checking -----------------------------------------------------
    if (overwrite == FALSE && dir.exists(path)) {
        stop("Output directory exists (set `overwrite = TRUE` to bypass)")
    }
    if (background != "automatic" && local_background) {
        stop("`background` and `local_background` are mutually exclusive; use only one")
    }
    if (only_known && only_denovo) {
        stop("Both `only_known` and `only_denovo` set to `TRUE`; pick one")
    }
        
    ## File format for `x` and `background` (bed vs data.frame)
    ## Write bed files if data.frame; assign if a character
    if ("data.frame" %in% class(x)) {
        target_bed <- tempfile("target_")
        .write_bed(x, path = target_bed)
    } else {
        if (!file.exists(x)) {
            stop("Check that your bed file for `x` exists")
        }        
        target_bed <- x
    }
    if (!("automatic" %in% background)) {
        if ("data.frame" %in% class(background)) {
            background_bed <- tempfile("background_")
            .write_bed(background, path = background_bed)
        } else {
            if (!file.exists(background)) {
                stop("Check that your bed file for `background` exists")
            }        
            background_bed <- background
        }
    }

    ## Run findMotifsGenome.pl ---------------------------------------------
    ## Make HOMER results output dir
    system(paste("mkdir -p", path))
    
    ## Construct and run command
    homer_base <- get_homer_bin()
    
    cmd <- paste(
        paste0(homer_base, "findMotifsGenome.pl"),
        target_bed, genome, path,
        "-len", paste0(motif_length, collapse = ","),
        "-size", scan_size,
        "-S", optimize_count,
        "-p", cores,
        "-cache", cache,
        "-fdr", fdr_num
    )
    if (!("automatic" %in% background)) {
        cmd <- paste(cmd, "-bg", background_bed)
    }
    if (local_background) {
        cmd <- paste(cmd, "-local", local_background)
    }
    if (only_known) {
        cmd <- paste(cmd, "-nomotif")
    }
    if (only_denovo) {
        cmd <- paste(cmd, "-noknown")
    }
    if (scan_size == "given") {
        cmd <- paste(cmd, "-chopify")
    }
    if (scale_logos) {
        cmd <- paste(cmd, "-bits")
    }
    system(cmd, ignore.stderr = !showHomerStdErr)

    ## Remove extraneous files if desired
    if (keep_minimal) {
        extra_files <- c("homerResults.html",
                         "knownResults.html",
                         "homerMotifs.motifs*",
                         "motifFindingParameters.txt",
                         "seq.autonorm.tsv",
                         "*tmp*")
        extra_dirs <- c("homerResults",
                        "knownResults",
                        "randomizations")
        remove_extra <- paste(c(paste0("rm -f ", path, "/", extra_files),
                                paste0("rm -Rf ", path, "/", extra_dirs)),
                              collapse = "; ")
        system(remove_extra)
    }

    ## Cleanup tmp files that are created by HOMER in cwd
    ## Should make sure these were created by HOMER..
    system("rm -f *.tmp")
    
}

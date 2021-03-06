% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_motifs_genome.R
\name{find_motifs_genome}
\alias{find_motifs_genome}
\title{Find Motifs over Regions}
\usage{
find_motifs_genome(
  x,
  path,
  genome,
  motif_length = c(8, 10, 12),
  scan_size = 100,
  optimize_count = 8,
  background = "automatic",
  local_background = FALSE,
  only_known = FALSE,
  only_denovo = FALSE,
  fdr_num = 0,
  cores = parallel::detectCores(),
  cache = .calc_free_mem()/4,
  overwrite = FALSE,
  keep_minimal = FALSE,
  scale_logos = FALSE
)
}
\arguments{
\item{x}{\code{data.frame} with the first three columns being
chromosome, start, and end coordinates, with a fourth column corresponding
to a region identifier; extra columns may be kept; x may alternately be a
path to an existing bed file of this format}

\item{path}{where to write HOMER results}

\item{genome}{ID of installed genome; check installed genomes using
\code{list_homer_packages()}; examples include "hg38" and "mm10";
add an 'r' at the end to mask repeats, e.g. "mm10r"}

\item{motif_length}{vector of motif lengths to consider [default is
\code{c(8, 10, 12)}]}

\item{scan_size}{size of sequence to scan; this can be a numeric to
specify the number of bases to scan centered on the region, or
alternately can be set to "given" to scan the entire region;
if using "given", will use the "-chopify" option to cut large
background sequences to average of target sequence size
[default: \code{100}]}

\item{optimize_count}{number of motifs to optimize [default: 8]}

\item{background}{\code{data.frame} containing coordinates of desired
regions to use as the background; alternately may be a path to an
existing bedfile; the default, "automatic", creates a background based
on the GC content and (scan) size of the target sequences}

\item{local_background}{a numeric scalar specifying number of equal size
regions around peaks to use as the local background
[by default this is not used, e.g. default: \code{FALSE}]}

\item{only_known}{turns off searching for denovo motifs}

\item{only_denovo}{turns off searching for known motif enrichment}

\item{fdr_num}{number of randomizations to perform to calculate FDR [default: 0]}

\item{cores}{number of cores to use [default: all cores available]}

\item{cache}{number in MB to use as cache to store sequences in memory
[default: calculates free memory and divides by 4]}

\item{overwrite}{overwrite an existing HOMER results directory
[default: \code{FALSE}]}

\item{keep_minimal}{remove all extra clutter from results, keep only
the essentials (\code{knownResults.txt} and \code{homerMotifs.all.motifs}
[default: \code{FALSE}]}

\item{scale_logos}{whether to scale sequence logos by information content 
[default: \code{FALSE}]}
}
\value{
Nothing; called for its side-effect of producing HOMER results
}
\description{
Calls \code{findMotifsGenome.pl} to run motif analysis over a given
set of regions.
}
\details{
\code{find_motifs_genome} runs the core HOMER motif enrichment function
from the R system, and in the process generates (as a side-effect)
a HOMER results directory.

This results directory is inspectable via a file browser, and contains
a summary of the results as HTML files as well as text files.

For our purposes, within the directory two key files exist:

\itemize{
  \item \code{knownResults.txt} known motifs that are enriched
  \item \code{homerResults.all.motif} denovo motifs that are enriched
}

These two text files are the core results (all else can be discarded
by setting \code{keep_minimal} to \code{TRUE}, and are parsed downstream
by \code{\link{read_known_results}} and \code{\link{read_denovo_results}}.
}
\seealso{
\code{\link{read_known_results}}, \code{\link{read_denovo_results}}
}

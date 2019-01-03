# marge <img src="man/figures/Marge_Simpson.png" align="right" height=205 width=132 />

## Overview 

The aim of `marge` is to provide an R API to `HOMER` for the analysis of genomic data, utilizing a tidy framework to accelerate organization and visualization of analyses.

## Installing `HOMER`

First, running `marge` requires having a working installation of `HOMER` on your computer. Please see the [HOMER website](http://homer.ucsd.edu/homer/) for more information on installing and configuring `HOMER` and to learn more about the methodology. In particular, note that you should install your desired genomes in addition to installing `HOMER` using the `./configureHomer.pl` script.

Note that working with a conda installation of `HOMER` is not well tested at this time. A potential workaround is below. I recommend installing directly from source.

## Installing `marge`

To install the latest development version of `marge`, simply do:


```r
devtools::install_github('robertamezquita/marge', ref = 'master')
```

## Before Running `marge`

While `marge` will do its best to find `HOMER`, there are certain environments where it will not be able to do so, specifically, with regards to RStudio and conda installs of `HOMER`. In these cases, a custom path to `HOMER` can be provided if it is not found by the package's utilities by setting `options('homer_path' = "/path/to/homer-4.10")`. This can be set in your `~/.Rprofile` so it loads the correct path automagically each time. 

## Usage

`marge` is currently semi-stable. The package currently includes the ability to:

* Run a motif analysis: `find_motifs_genome()` - runs the `HOMER` script `findMotifsGenome.pl` via R, and outputs a results directory in the default `HOMER` style
* Read in results: `read_*_results()` - read in either `denovo` or `known` enriched motifs with the `read_denovo_results()` or `read_known_results()` functions, pointing to the `HOMER` directory that was created in the previous step. The `read_*` functions produce tibbles summarizing the motif enrichment results into a tidy format for easier visualization and analysis. See the reference pages of each for more details.
* Write motifs in HOMER compatible format with `write_homer_motif()`
* Find specific motif instances across regions using supplied PWMs with `find_motifs_instances()` and read in the results with `read_motifs_instances()`
* Access the HOMER database of known motifs by inspecting the `HOMER_motifs` object

Further details can be found in the associated vignette, describing installation and typical workflows encompassing basic/advanced usage schemas.

## Compared to HOMER alone

Like the actual Homer Simpson, `HOMER` is made better with the addition of `marge`. With the continually increasing throughput in conducting sequencing analysis, `marge` provides a native R framework to work from end to end with motif analyses - from processing to storing to visualizing these results, all using modern tidy conventions.

---

#### Contact

* Robert A. Amezquita (aut, cre, cph): [robert.amezquita@fredhutch.org](robert.amezquita@fredhutch.org)

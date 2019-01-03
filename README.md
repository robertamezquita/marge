# marge <img src="man/figures/Marge_Simpson.png" align="right" height=205 width=132 />

## Overview 

The aim of [`marge`](https://marge.aerobatic.io) is to provide an R API to `HOMER` for the analysis of genomic data, utilizing a tidy framework to accelerate organization and visualization of analyses.

If here courtesy of Bitbucket, check out the docs at: [`marge`](https://marge.aerobatic.io)

## Installing `HOMER`

First, running `marge` requires having a working installation of `HOMER` on your computer. Please see the [HOMER website](http://homer.ucsd.edu/homer/) for more information on installing and configuring `HOMER` and to learn more about the methodology. In particular, note that you should install your desired genomes in addition to installing `HOMER` using the `./configureHomer.pl` script.

## Installing `marge`

To install the latest development version of `marge`, navigate to the `marge` bitbucket [downloads page](https://bitbucket.org/robert_amezquita/marge/downloads/) to download and build, or simply do:


```r
devtools::install_bitbucket('robert_amezquita/marge', ref = 'master')
```

To install a stable version, simply navigate to the [downloads page](https://bitbucket.org/robert_amezquita/marge/downloads/), navigate to the tab called "Tags", and change the `ref` argument from `master` to your desired release (for example, `v0.0.3_carl`). 

## Usage

`marge` is currently in semi-active development, the package currently includes the ability to:

* Run a motif analysis: `find_motifs_genome()` - runs the `HOMER` script `findMotifsGenome.pl` via R, and outputs a results directory in the default `HOMER` style
* Read in results: `read_*_results()` - read in either `denovo` or `known` enriched motifs with the `read_denovo_results()` or `read_known_results()` functions, pointing to the `HOMER` directory that was created in the previous step. The `read_*` functions produce tibbles summarizing the motif enrichment results into a tidy format for easier visualization and analysis. See the reference pages of each for more details.
* Write motifs in HOMER compatible format with `write_homer_motif()`
* Find specific motif instances across regions using supplied PWMs with `find_motifs_instances()` and read in the results with `read_motifs_instances()`
* Access the HOMER database of known motifs by inspecting the `HOMER_motifs` object

Further details can be found in the associated [vignette](https://marge.aerobatic.io/articles/marge-workflow), describing installation and typical workflows encompassing basic/advanced usage schemas.

## Compared to HOMER alone

Like the actual Homer Simpson, `HOMER` is made better with the addition of `marge`. With the continually increasing throughput in conducting sequencing analysis, `marge` provides a native R framework to work from end to end with motif analyses - from processing to storing to visualizing these results, all using modern tidy conventions.

---

#### Contact

* Robert A. Amezquita (aut, cre, cph): [robert.amezquita@fredhutch.org](robert.amezquita@fredhutch.org)

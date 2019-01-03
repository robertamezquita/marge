# marge 0.0.4

* added more information regarding installation hiccups and finding HOMER
* migrating to Github for long-term, removing references to Bitbucket

# marge 0.0.3

* completely reworked how `marge` checks for HOMER installation - 
  * Uses `options('homer_path')` to hold the installation base path
  * Tries the default (e.g. if HOMER is a part of `$PATH`) first, susses out HOMER base
  * Failing that, asks for user to manually set it, and provides code to set in `.Rprofile`
  * Instructions encoded into messages following success/failure
* added an `onAttach.R` to check for proper HOMER install on load with link to online docs


# marge 0.0.2

* lots of little fixes for various functions, code cleanup


# marge 0.0.1.9999

* bug fix for `find_motifs_genome` when using a custom background - added an `fdr_num` argument to add number of iterations for FDR calculation (background file was being fed to `-fdr` as the `<num>` argument); set to 0 by default since runtime takes forever when set when doing denovo enrichment
* reads in new `fdr` field for denovo results, as the `fdr` arg for `find_motifs_genome` has now been enabled by default always
* bug fix for `find_motifs_genome` - `motif_length` arg was only using the first number
* added `find_motifs_instances` function to get where specific motifs (based on PWM) occur in peaks
* added `write_pwm` function to support `find_motifs_instances`
* created new datafile `HOMER_motifs` which includes all motifs used by HOMER - includes their PWMs, organism/motif set definition info
* removed placeholder functions `[read|write]_bed` that were initially included for testing
* added a new vignette


# marge 0.0.1

* initial release of `marge` which includes the most basic utilities needed for analysis (see below)
* `find_motifs_genome` - runs motif analysis using `HOMER`; creates a `HOMER` results directory
* `read_denovo_results` - reads the `homerMotifs.all.motifs` file created in the `HOMER` results directory
* `read_known_results` - reads the `knownResults.txt` file created in the `HOMER` results directory

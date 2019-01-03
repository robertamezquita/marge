## Read in the motifs from HOMER DB
library(tidyverse)
devtools::load_all()

## Construct paths to HOMER db motif sets (all)
base <- marge::get_homer_base() %>% dirname()
mset <- c('vertebrates', 'insects', 'plants', 'worms', 'yeast')
paths <- paste0(base, '/data/knownTFs/', mset, '/all.motifs')

## Function to read in each motif and add 'mset' (motif set) var
.read_org_motifs <- function(path, mset) {
    read_denovo_results(
        path = path,
        homer_dir = FALSE
    ) %>%
        data.frame(mset = mset, .) %>%
        as_tibble
}

org_motifs <- map2(paths, mset, .read_org_motifs) %>%
    bind_rows()

## Read in additional motifs from all without an mset
all <- read_denovo_results(paste0(base, '/data/knownTFs/all.motifs'), homer_dir = FALSE)
amiss <- all[!(all$consensus %in% org_motifs$consensus), ]
amiss$mset <- 'various'

HOMER_motifs <- bind_rows(org_motifs, amiss)


## Save data
usethis::use_data(HOMER_motifs, overwrite = TRUE, compress = 'xz')



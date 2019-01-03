context("Test running HOMER")

test_that("Check test data", {
    expect_s3_class(test_data, "data.frame")
})

test_that("Running HOMER findMotifsGenome.pl", {
    skip("Need to limit time for `system` call to `findMotifsGenome.pl`")
    
    ## Check if HOMER installed
    skip_if_not(check_homer() == 0, "HOMER not installed")
    
    ## Check if genomes are installed
    packages <- list_homer_packages()
    genomes <- packages %>%
        dplyr::filter_("package %in% c('hg38', 'mm10')")
    skip_if_not(nrow(genomes) > 0, "Skipping tests..genomes (hg38/mm10) not installed")
    genome <- genomes$package[1]
    
    ## Set temp path for outputs
    path <- tempfile("test_homer_")

    find_motifs_genome(test_data,
                       path = path,
                       genome = genome,
                       motif_length = 5,
                       scan_size = 10,
                       optimize_count = 1,
                       background = "automatic",
                       local_background = FALSE,
                       only_known = FALSE,
                       only_denovo = FALSE,
                       cores = 1,
                       cache = 500,
                       overwrite = TRUE)
  
    expect_true(file.exists(paste0(path, "/homerMotifs.all.motifs")))
    expect_true(file.exists(paste0(path, "/knownResults")))
})

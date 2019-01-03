context("Test reading HOMER files")

test_that("Check read_known_results", {
    file <- system.file("extdata", "", package = "marge")
    dat <- read_known_results(file)
    expect_s3_class(dat, "data.frame")
})

test_that("Check read_denovo_results", {
    file <- system.file("extdata", "", package = "marge")
    dat <- read_denovo_results(file)
    expect_s3_class(dat, "data.frame")
})

test_that("Check reading complete/incomplete motif files", {
    ## Incomplete motif
    ap1 <- system.file("extdata/ap1.motif", package = "marge")
    expect_s3_class(read_motif(ap1), "data.frame")
    ## Complete motif
    znf7 <- system.file("extdata/znf7.motif", package = "marge")
    expect_s3_class(read_motif(znf7), "data.frame")
    ## Check sub function
    expect_s3_class(.parse_homer_subfields(read_motif(ap1)), "data.frame")
    expect_s3_class(.parse_homer_subfields(read_motif(znf7)), "data.frame")
})

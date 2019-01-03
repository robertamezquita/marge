context("HOMER install")

test_that("HOMER is present on system", {
    skip_if_not(check_homer() == 0, "Skipping tests..HOMER not installed")
    expect_message(check_homer())
    expect_length(get_homer_base(), 1)
    expect_s3_class(list_homer_packages(), "tbl")
})




## Create example data for testing purposes
test_data <- tibble::tibble(
    chr   = c("chr6", "chr7", "chr7", "chr12"),
    start = c(102207441, 44688352, 92093830, 103669461),
    end   = c(102207830, 44688664, 92094134, 103669829),
    row   = c("region_1", "region_2", "region_3", "region_4"),
    lfc   = c(1, 2, 3, 4),
    gene  = c("A", "B", "C", "D")
)



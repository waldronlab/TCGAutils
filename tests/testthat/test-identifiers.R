context("Identifier tests")

.sectionNums <- function(bcode) {
    filler <- .uniqueDelim(bcode)
    unique(lengths(strsplit(bcode, filler)))
}

test_that("TCGAbarcode works", {
    example("TCGAbarcode")
    expect_identical(.sectionNums(TCGAbarcode(barcodes)), 3L)

    expect_identical(.sectionNums(TCGAbarcode(barcodes, sample = TRUE)), 4L)

    expect_identical(
        .sectionNums(
            TCGAbarcode(barcodes, sample = TRUE, portion = TRUE)), 5L)

    expect_identical(
        .sectionNums(
            TCGAbarcode(barcodes, sample = TRUE, portion = TRUE, plate = TRUE)),
        6L)
    expect_identical(
        .sectionNums(
            TCGAbarcode(barcodes, sample = TRUE, portion = TRUE,
                plate = TRUE, center = TRUE)),
        7L)
})

test_that("TCGAbiospec works", {
    bc0 <- TCGAbarcode(barcodes)
    expect_error(TCGAbiospec(bc0))
    bc1 <- TCGAbarcode(barcodes, sample = TRUE)
    expect_identical(dim(TCGAbiospec(bc1)), c(length(bc1), .sectionNums(bc1)))
    bc2 <- TCGAbarcode(barcodes, sample = TRUE, portion = TRUE)
    expect_identical(dim(TCGAbiospec(bc2)), c(length(bc2),
        .sectionNums(bc2)+1L))
    bc3 <- TCGAbarcode(barcodes, sample = TRUE, portion = TRUE, plate = TRUE)
    expect_identical(dim(TCGAbiospec(bc3)), c(length(bc3),
        .sectionNums(bc3)+1L))
    bc4 <- TCGAbarcode(barcodes, sample = TRUE, portion = TRUE,
        plate = TRUE, center = TRUE)
    expect_identical(dim(TCGAbiospec(bc4)), c(length(bc4),
        .sectionNums(bc4)+1L))
    expect_identical(names(TCGAbiospec(barcodes)), c("submitter_id",
        "sample_definition", "sample", "vial", "portion", "analyte", "plate",
        "center"))
})

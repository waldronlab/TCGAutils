context("ID translation testing")

test_that("barcodeToUUID translates correctly", {
    pt1 <- c("TCGA-06-6391", "TCGA-06-6700")
    case_id <- barcodeToUUID(pt1)
    expect_true("case_id" %in% names(case_id))
    expect_equal(
        case_id[["case_id"]],
        c("7a4c0a14-ac97-4c2b-a9cc-68cb561b2494",
        "3dddfc44-7bb1-4974-8a65-a84fd4bac484")
    )
    analytes <- c("TCGA-AA-A00L-10A-01X", "TCGA-AA-A00L-10A-01D")
    analyte_ids <- barcodeToUUID(analytes)
    expect_true("analyte_ids" %in% names(analyte_ids))
    expect_identical(
        analyte_ids[["analyte_ids"]],
        c("2f3031e8-7dac-4444-9cf7-a5dff4751b9e",
        "1c429d23-89eb-4c35-bef3-9eff2508d9d5")
    )
    portions <- c("TCGA-AA-A00L-10A-01", "TCGA-AA-A00L-01A-31")
    portion_ids <- barcodeToUUID(portions)
    expect_true("portion_ids" %in% names(portion_ids))
    expect_identical(
        portion_ids[["portion_ids"]],
        c("c72ff462-a355-49fa-8275-c34ef5dd91c9",
        "7d25aecc-9068-463b-adc5-71ec2f4ba7aa")
    )
})

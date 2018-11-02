context("ID translation testing")

test_that("barcodeToUUID translates correctly", {
    pt1 <- c("TCGA-06-6391", "TCGA-06-6700")
    caseids <- barcodeToUUID(pt1)[["case_id"]]
    expect_identical(
        caseids,
        c("7a4c0a14-ac97-4c2b-a9cc-68cb561b2494",
        "3dddfc44-7bb1-4974-8a65-a84fd4bac484")
    )
})

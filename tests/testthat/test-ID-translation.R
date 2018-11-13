context("ID translation testing")

test_that("barcodeToUUID translates correctly", {
    pts <- c("TCGA-06-6391", "TCGA-06-6700")
    case_id <- barcodeToUUID(pts)
    expect_true("case_id" %in% names(case_id))
    expect_equal(
        case_id[["case_id"]],
        c("7a4c0a14-ac97-4c2b-a9cc-68cb561b2494",
        "3dddfc44-7bb1-4974-8a65-a84fd4bac484")
    )
    samps <- c("TCGA-06-6700-01A", "TCGA-AD-6888-01A")
    samp_id <- barcodeToUUID(samps)
    expect_true("sample_ids" %in% names(samp_id))
    expect_equal(
        samp_id[["sample_ids"]],
        c("8d35786c-5edb-4a84-b3e5-c401b8c73bd6",
        "ecf0f65b-bf3c-4d0e-899a-f209247cbe97")
    )
    analytes <- c("TCGA-AA-A00L-10A-01X", "TCGA-AA-A00L-10A-01D",
        "TCGA-12-0653-10A-01D")
    analyte_ids <- barcodeToUUID(analytes)
    expect_true("analyte_ids" %in% names(analyte_ids))
    expect_equal(
        analyte_ids[["analyte_ids"]],
        c("4b6a77dc-7a2a-459e-a7a0-253f950f1c8c",
        "1c429d23-89eb-4c35-bef3-9eff2508d9d5",
        "63645523-bb46-40b3-899b-c3fa5fefd121")
    )
    portions <- c("TCGA-AA-A00L-10A-01", "TCGA-AA-A00L-01A-31",
        "TCGA-12-0653-10A-01")
    portion_ids <- barcodeToUUID(portions)
    expect_true("portion_ids" %in% names(portion_ids))
    expect_equal(
        portion_ids[["portion_ids"]],
        c("c72ff462-a355-49fa-8275-c34ef5dd91c9",
        "7d25aecc-9068-463b-adc5-71ec2f4ba7aa",
        "03209a36-67a0-48df-a9f7-a0cedd0db82f")
    )
    aliquots <- c("TCGA-12-0653-10A-01D-0333-01",
        "TCGA-12-0653-10A-01D-0334-04", "TCGA-AA-3556-01A-01D-1953-10")
    aliquot_ids <- barcodeToUUID(aliquots)
    expect_true("aliquot_ids" %in% names(aliquot_ids))
    expect_equal(aliquot_ids[["aliquot_ids"]],
        c("51ddbc44-1cae-454f-bc67-5c5cc3d9e853",
        "2f0fe3f0-6a24-47ee-acba-df9c04d89532",
        "2303247f-9691-4b38-bac2-8a30d6e08cc9")
    )
})

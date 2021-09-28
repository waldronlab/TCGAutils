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


test_that("UUIDtoBarcode translates correctly", {
    exuuid <- "6b7d7a7f-f16d-472d-9b7b-3482c434cc99"
    entities <- c(
        "TCGA-ND-A4W6-01A-11D-A28Q-01", "TCGA-N9-A4Q8-01A-31D-A28Q-01",
        "TCGA-NG-A4VW-01A-11D-A28Q-01", "TCGA-N9-A4Q4-01A-11D-A28Q-01",
        "TCGA-NA-A4QW-01A-11D-A28Q-01", "TCGA-N6-A4VF-01A-31D-A28Q-01",
        "TCGA-N5-A4RD-01A-11D-A28Q-01", "TCGA-N8-A56S-01A-11D-A28Q-01",
        "TCGA-N5-A4RS-01A-11D-A28Q-01", "TCGA-NA-A5I1-01A-21D-A28Q-01",
        "TCGA-ND-A4WF-01A-11D-A28Q-01", "TCGA-N5-A4RU-01A-31D-A28Q-01",
        "TCGA-N5-A4RT-01A-11D-A28Q-01", "TCGA-N5-A4RA-01A-11D-A28Q-01",
        "TCGA-ND-A4WA-01A-12D-A28Q-01", "TCGA-NA-A4R0-01A-11D-A28Q-01",
        "TCGA-N8-A4PN-01A-11D-A28Q-01", "TCGA-N9-A4PZ-01A-22D-A28Q-01",
        "TCGA-N5-A4RF-01A-11D-A28Q-01", "TCGA-NG-A4VU-01A-11D-A28Q-01",
        "TCGA-QN-A5NN-01A-11D-A28Q-01", "TCGA-N9-A4Q1-01A-11D-A28Q-01",
        "TCGA-N6-A4VC-01A-11D-A28Q-01", "TCGA-N7-A4Y0-01A-12D-A28Q-01",
        "TCGA-N8-A4PQ-01A-11D-A28Q-01", "TCGA-N5-A4RJ-01A-11D-A28Q-01",
        "TCGA-N9-A4Q3-01A-11D-A28Q-01", "TCGA-N5-A4RM-01A-11D-A28Q-01",
        "TCGA-N8-A4PI-01A-21D-A28Q-01", "TCGA-N6-A4VG-01A-31D-A28Q-01",
        "TCGA-N6-A4V9-01A-11D-A28Q-01", "TCGA-N5-A59F-01A-11D-A28Q-01",
        "TCGA-N7-A4Y5-01A-12D-A28Q-01", "TCGA-NA-A4QX-01A-11D-A28Q-01",
        "TCGA-NA-A4QV-01A-11D-A28Q-01", "TCGA-N6-A4VE-01A-11D-A28Q-01",
        "TCGA-N6-A4VD-01A-11D-A28Q-01", "TCGA-N8-A4PM-01A-11D-A28Q-01",
        "TCGA-NF-A4X2-01A-11D-A28Q-01", "TCGA-N8-A4PL-01A-11D-A28Q-01",
        "TCGA-N5-A59E-01A-11D-A28Q-01", "TCGA-N8-A4PO-01A-11D-A28Q-01",
        "TCGA-ND-A4WC-01A-21D-A28Q-01", "TCGA-N7-A4Y8-01A-11D-A28Q-01",
        "TCGA-NF-A4WU-01A-11D-A28Q-01", "TCGA-N7-A59B-01A-11D-A28Q-01",
        "TCGA-QM-A5NM-01A-11D-A28Q-01", "TCGA-N5-A4RV-01A-21D-A28Q-01",
        "TCGA-NA-A4QY-01A-11D-A28Q-01", "TCGA-N9-A4Q7-01A-11D-A28Q-01",
        "TCGA-NF-A4WX-01A-11D-A28Q-01", "TCGA-NF-A5CP-01A-12D-A28Q-01",
        "TCGA-NA-A4R1-01A-11D-A28Q-01", "TCGA-N5-A4RO-01A-11D-A28Q-01",
        "TCGA-N8-A4PP-01A-11D-A28Q-01", "TCGA-N5-A4RN-01A-12D-A28Q-01"
    )
    resframe <- UUIDtoBarcode(exuuid, from_type = "file_id")
    expect_true(
        all(entities %in% resframe[["associated_entities.entity_submitter_id"]])
    )

    aliquot_uuids <- c(
        "f8c7d038-1182-42d0-8787-b84b5ca57eaf",
        "b37ea112-340e-4613-8514-d8a8bd47410f",
        "4a9967bf-444c-4573-a082-121a30be7f3b"
    )
    aliquot_submitters <- c(
        "TCGA-UF-A71A-06A-11D-A390-01",
        "TCGA-BB-4224-01A-01D-1432-01",
        "TCGA-CN-4735-01A-01D-1432-01"
    )
    resframe <- UUIDtoBarcode(aliquot_uuids, from_type = "aliquot_ids")
    expect_true(
        all(
            aliquot_uuids %in% resframe[[1L]]
        )
    )
    expect_true(
        all(
            aliquot_submitters %in% resframe[[2L]]
        )
    )
})


test_that("UUIDtoBarcode shows multiple entries per file_id", {

    file_ids <- c(
        "01ef8a08-1de5-4ceb-be51-979418465f1a",
        "81179af5-7b83-46ec-814d-877e97c884a2",
        "edce3d0c-efad-4b00-9a53-f26d714b755d"
    )

    restabs <- lapply(file_ids, UUIDtoBarcode, "file_id")
    results <- do.call(rbind, restabs)

    expect_identical(results, UUIDtoBarcode(file_ids, "file_id"))

    file_ids <- c(
        "01ef8a08-1de5-4ceb-be51-979418465f1a",
        "82279af5-7b83-46ec-814d-877e97c884a2",
        "edce3d0c-efad-4b00-9a53-f26d714b755d"
    )

    expect_warning(UUIDtoBarcode(file_ids, "file_id"))
})

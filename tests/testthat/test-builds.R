context("Build information testing")

test_that("translateBuild works correctly", {
    buildDF <- human_builds()
    ncbinos <- as.character(34:38)
    resbuilds <- vapply(ncbinos, translateBuild, character(1L))

    expect_identical(unname(resbuilds), buildDF[["UCSC"]])

    ucscnos <- paste0("hg", c(16:19, 38))
    resbuilds <- vapply(ucscnos, translateBuild, character(1L), "NCBI")
    expect_identical(unname(resbuilds), buildDF[["NCBI"]])

    ## UCSC (default 'to')
    expect_identical(translateBuild("Grch37"), "hg19")
    expect_identical(translateBuild("GrCh37"), "hg19")
    expect_identical(translateBuild("grch37"), "hg19")

    expect_identical(
        translateBuild("hg19", to = "NCBI"),
        "GRCh37"
    )
    expect_identical(
        translateBuild("HG19", to = "NCBI"),
        "GRCh37"
    )
    expect_identical(
        translateBuild("hG19", to = "NCBI"),
        "GRCh37"
    )
    expect_true(
            is.na(translateBuild(NA_character_))
    )
    expect_true(
            is.na(translateBuild("33"))
    )
})


test_that("uniformBuilds is returning the appropriate output", {
    build <- rep(c("GRCh37", "hg19"), times = c(5, 1))
    rebuild <- uniformBuilds(build)
    expect_identical(1L, length(unique(rebuild)))

    ## NA imputed to rest of builds
    build <- c(rep(c("GRCh37", "hg19"), times = c(5, 1)), "NA")
    rebuild <- uniformBuilds(build)
    expect_identical(1L, length(unique(rebuild)))

    build <- c(rep(c("GRCh37", "hg19"), times = c(2, 1)), "NA")
    expect_error(uniformBuilds(build, cutoff = 0.2))

    # NA prop > 0.2
    build <- c(rep(c("GRCh37", "hg19"), times = c(7, 1)), "NA", "NA")
    expect_error(uniformBuilds(build, cutoff = 0.2))

    # NA converted to main build annotation
    build <- c(rep(c("GRCh37", "hg19"), times = c(7, 2)), NA_character_)
    rebuild <- uniformBuilds(build, cutoff = 0.2)
    expect_identical(1L , length(unique(rebuild)))

    # if build numbers identical then replace with high prop
    build <- rep(c("GRCh37", "37"), times = c(7, 2))
    rebuild <- uniformBuilds(build, cutoff = 0.2)
    expect_identical(rebuild, rep("GRCh37", length(rebuild)))

    build <- c(rep(c("GRCh37", "37"), times = c(7, 2)), NA_character_)
    rebuild <- uniformBuilds(build, cutoff = 0.2)
    expect_identical(rebuild, rep("GRCh37", length(rebuild)))

    build <- c(rep(c("GRCh37", "37"), times = c(7, 2)), rep(NA_character_, 3))
    expect_error(uniformBuilds(build, cutoff = 0.2))
})


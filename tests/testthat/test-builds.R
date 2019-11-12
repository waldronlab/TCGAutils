context("Build information testing")

test_that("translateBuild works correctly", {
    ncbinos <- as.character(34:38)
    ucscnos <- paste0("hg", c(16:19, 38))
    resbuilds <- vapply(ncbinos, translateBuild, character(1L))

    Map(function(x, y) {

        expect_identical(x, y)

    }, ucscnos, unname(resbuilds))

    expect_identical(translateBuild("Grch37"), "hg19")
    expect_identical(translateBuild("GrCh37"), "hg19")
    expect_identical(
        translateBuild("hg19", to = "NCBI"),
        "GRCh37"
    )
    expect_true(is.na(translateBuild(NA_character_)))
    expect_warning(
        translateBuild(NA_character_),
        "build could not be matched"
    )
    expect_true(is.na(translateBuild("33")))
    expect_warning(
        translateBuild("33"),
        "build could not be matched"
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


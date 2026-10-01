## A few small helpers live in R/shared-helpers.R, the same file byte for
## byte in illume and illumex (dev/check_shared_helpers.R compares the two
## repositories), rather than one package reaching into the other's
## internals. This checks the installed illumex agrees, so a change made to
## one copy and not the other fails here as well. The code is compared, not
## the function objects, whose environments are different namespaces by
## construction.

test_that("the helpers illume shares with illumex are the same code", {
  shared <- c("%||%", "ilm_bq", "ilm_wrap", "ilm_and", "ilm_progress",
              "ilm_pch", "ilm_pch_spec", "ilm_pch_table", "ilm_pch_names",
              "ilm_rng_restore", "ILM_POW10", "ilm_disp_scale",
              "ilm_disp_prod_err", "ilm_disp_mantissa", "ilm_disp_digits",
              "ilm_disp_value", "ilm_disp_round", "ilm_disp_signif_digits",
              "ilm_disp_signif", "ilm_disp_text", "ilm_disp")
  for (f in shared) {
    mine <- get(f, envir = asNamespace("illume"))
    theirs <- get(f, envir = asNamespace("illumex"))
    expect_identical(deparse(mine), deparse(theirs),
                     label = paste0(f, " in illume"),
                     expected.label = paste0(f, " in illumex"))
  }
})

test_that("a list of nothing is said as nothing", {
  and <- get("ilm_and", envir = asNamespace("illume"))
  expect_identical(and(character(0)), "")
  expect_identical(and(""), "")
  expect_identical(and("a"), "a")
  expect_identical(and(c("a", "b", "c")), "a, b and c")
})

test_that("no name the shared file defines is defined anywhere else in R/", {
  ## A name defined both in shared-helpers.R and in another file is silently
  ## overwritten by whichever file collates later: illume's ilm_disp_scale(),
  ## the dispersion on its natural scale, was replaced by the shared file's
  ## rounding helper of that name, and the marginal effects' standard errors
  ## came out NA. The namespace keeps only the winner, so the sources are
  ## parsed: beside the tests in a source tree, and in 00_pkg_src under
  ## R CMD check.
  cand <- c(test_path("..", "..", "R"),
            test_path("..", "..", "00_pkg_src", "illume", "R"))
  rdir <- cand[dir.exists(cand)][1]
  skip_if(is.na(rdir), "the package's R/ sources are not beside the tests")
  top <- function(f) unlist(lapply(parse(f, keep.source = FALSE), function(e)
    if (is.call(e) && is.name(e[[1]]) &&
        as.character(e[[1]]) %in% c("<-", "=", "<<-") &&
        (is.name(e[[2]]) || is.character(e[[2]]))) as.character(e[[2]])))
  fs <- list.files(rdir, "[.][Rr]$", full.names = TRUE)
  sh <- fs[basename(fs) == "shared-helpers.R"]
  expect_length(sh, 1L)
  others <- unlist(lapply(setdiff(fs, sh), top))
  expect_identical(intersect(top(sh), others), character(0))
})

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
              "ilm_rng_restore")
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

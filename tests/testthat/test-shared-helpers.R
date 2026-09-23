## illume keeps its own copies of a few small helpers that illumex also has,
## rather than reaching into illumex's internals. A change made to one copy
## and not the other would leave the two packages quietly disagreeing -- about
## a plotting character, a backticked name, a progress bar -- so this fails the
## moment they differ. The code is compared, not the function objects, whose
## environments are different namespaces by construction.

test_that("the helpers illume shares with illumex are the same code", {
  shared <- c("%||%", "ilm_bq", "ilm_wrap", "ilm_progress",
              "ilm_pch", "ilm_pch_spec", "ilm_pch_table", "ilm_pch_names")
  for (f in shared) {
    mine <- get(f, envir = asNamespace("illume"))
    theirs <- get(f, envir = asNamespace("illumex"))
    expect_identical(deparse(mine), deparse(theirs),
                     label = paste0(f, " in illume"),
                     expected.label = paste0(f, " in illumex"))
  }
})

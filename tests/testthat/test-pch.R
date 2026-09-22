## Plotting characters by name
##
## The value of naming them is that a name can be CHECKED and a number cannot:
## pch = 26 is accepted by graphics and draws nothing, and 16 against 19 is a
## difference no reader can interpret. So these tests care about two things --
## that every name lands on the code base R actually uses, and that a wrong one
## stops rather than drawing something arbitrary.

test_that("the names land on the codes base R uses", {
  expect_identical(ilm_pch("filled circle"), 16L)
  expect_identical(ilm_pch("open circle"), 1L)
  expect_identical(ilm_pch("open square"), 0L)
  expect_identical(ilm_pch("filled square"), 15L)
  expect_identical(ilm_pch("filled triangle"), 17L)
  expect_identical(ilm_pch("filled diamond"), 18L)
  expect_identical(ilm_pch("plus"), 3L)
  expect_identical(ilm_pch("cross"), 4L)
  expect_identical(ilm_pch("star"), 8L)
  expect_identical(ilm_pch("circle fill"), 21L)
})

test_that("spelling it any reasonable way works", {
  for (v in c("filled circle", "Filled Circle", "FILLED CIRCLE",
              "filled_circle", "filled-circle", "filled.circle",
              "filledcircle", "  filled circle  "))
    expect_identical(ilm_pch(v), 16L, info = v)
})

test_that("every code from 0 to 25 has a name, and every alias reaches it", {
  ## a table that covers only the codes someone happened to need is a table
  ## people stop trusting
  spec <- ilm_pch_spec()
  codes <- vapply(spec, function(s) as.integer(s[1]), 0L)
  expect_identical(sort(codes), 0:25)
  for (s in spec)
    expect_identical(ilm_pch(s[-1]), rep(as.integer(s[1]), length(s) - 1L),
                     info = s[2])
})

test_that("numbers pass through, and impossible ones do not", {
  expect_identical(ilm_pch(16), 16)
  expect_identical(ilm_pch(c(16L, 17L)), c(16L, 17L))
  expect_identical(ilm_pch(0), 0)
  expect_identical(ilm_pch(25), 25)
  ## 26 draws nothing at all in base graphics, silently
  expect_error(ilm_pch(26), "0 to 25")
  expect_error(ilm_pch(-1), "0 to 25")
  expect_error(ilm_pch(3.5), "0 to 25")
})

test_that("a single character stays a glyph", {
  ## base R draws pch = "x" as the letter x. Translating it to 4 would be a
  ## silent change of plot, and "x" is in the table as an alias for cross.
  expect_identical(ilm_pch("x"), "x")
  expect_identical(ilm_pch("*"), "*")
  expect_identical(ilm_pch("o"), "o")
  expect_identical(ilm_pch(c("filled circle", "x")), c("16", "x"))
})

test_that("an unknown name stops and lists the ones that exist", {
  expect_error(ilm_pch("blob"), "not a plotting character")
  expect_error(ilm_pch("blob"), "filled circle")
})

test_that("NULL, NA and empty input are left alone", {
  expect_null(ilm_pch(NULL))
  expect_identical(ilm_pch(character(0)), character(0))
  expect_identical(ilm_pch(c("filled circle", NA)), c(16L, NA_integer_))
})

test_that("a name reaches the plotting functions themselves", {
  ## the helper being right is not the same as it being wired in
  skip_if_not_installed("tinyplot")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(ilm_plot_scatter(mtcars, "mpg", "wt", pch = "filled circle"))
  expect_silent(ilm_plot_line(mtcars, "mpg", "wt", pch = "open square"))
  expect_silent(ilm_plot_stat_error(mtcars, "mpg", "cyl",
                                    pch = "filled diamond"))
  expect_silent(ilm_plot_box(mtcars, "mpg", x = "cyl", pch = "star"))
  expect_error(ilm_plot_scatter(mtcars, "mpg", "wt", pch = "blob"),
               "not a plotting character")
  expect_silent(ilm_plot_scatter(mtcars, "mpg", "wt", pch = 17))
  expect_silent(ilm_plot_scatter(mtcars, "mpg", "wt", pch = "x"))
})

test_that("every function that documents `pch` actually takes one", {
  ## The first attempt routed all thirteen call sites through a do.call()
  ## funnel so the translation lived in one place. do.call() puts the
  ## EVALUATED arguments into the call, tinyplot deparses its arguments to
  ## title the legend, and a numeric `by` column came back as a deparsed
  ## 32-element vector that blew the legend width past what the device could
  ## take -- "invalid graphics state", from an example that had been in the
  ## package for months. It passed every test here and failed R CMD check.
  ##
  ## So the promise is kept per function instead, and this checks that the
  ## documentation and the signature agree: anything whose help says `pch`
  ## takes a name has to have a `pch` argument that reaches ilm_pch().
  fns <- c("ilm_plot", "ilm_plot_scatter", "ilm_plot_line",
           "ilm_plot_stat_error", "ilm_plot_box", "ilm_plot_violin")
  for (fn in fns) {
    f <- get(fn, envir = asNamespace("illume"))
    expect_true("pch" %in% names(formals(f)), info = fn)
    ## after `...`, so no existing positional call can be captured by it
    nm <- names(formals(f))
    expect_gt(match("pch", nm), match("...", nm))
    expect_match(paste(deparse(body(f)), collapse = " "), "ilm_pch", info = fn)
  }
})

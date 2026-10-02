## An AR(1) correlation is fitted as atanh(rho), which has no end. Past about
## 11 the stationary terms lose their digits and the objective falls away
## spuriously; a 12-point series from a simulation study went there, to
## atanh(rho) = -19, and was graded FAIL with a gradient of 11, against a
## true minimum at rho = -0.55. Short of the floor the objective flattens
## towards the edge, and the same data rounded to four places stopped there,
## at rho = -1 and an objective of 0.666, graded well. A fit near the edge or
## past the floor (ilm_rho_floor) is refitted inside it (ilm_rho_refit).

ar12 <- function() data.frame(g = factor("1"), t = 1:12,
  y = c(4.9582, 5.2958, 5.4162, 5.4551, 5.7439, 4.3653, 4.8777, 6.2922,
        4.8298, 5.9805, 5.5264, 6.7742))

test_that("a correlation left at the edge is brought back to the true minimum", {
  ## by maximum likelihood, where the study's series and its true minimum
  ## were found (a gaussian model with a correlation over time is fitted by
  ## REML by default, Craig's item 249)
  f <- suppressMessages(suppressWarnings(ilm_model(y ~ 1, data = ar12(), family = "gaussian",
                                                   ar = ilm_ar1(~ t | g), reml = FALSE,
                                                   verbose = FALSE)))
  pe <- f$opt$par
  expect_lte(abs(pe[["rho_raw"]]), illume:::ilm_rho_floor)
  expect_equal(f$rho, -0.545, tolerance = 0.01)
  expect_lt(max(abs(f$obj$gr(pe))), 1e-3)
  expect_true(f$ok)
  expect_identical(f$checks$status[f$checks$check == "gradient"], "OK")
})

test_that("the refit leaves a fit away from the edge as it was, and keeps a tie", {
  obj <- list(fn = function(p) sum((p - c(1, 2))^2), gr = function(p) 2 * (p - c(1, 2)),
              env = new.env())
  inside <- list(par = c(a = 1, rho_raw = 2), objective = 0)
  expect_identical(illume:::ilm_rho_refit(obj, inside, 2L, integer(0), numeric(0), list()), inside)
  ## an objective that does not depend on rho: the refit stays at the floor
  flat <- list(fn = function(p) (p[1] - 1)^2, gr = function(p) c(2 * (p[1] - 1), 0),
               env = new.env())
  past <- list(par = c(a = 1, rho_raw = 12), objective = 0)
  o <- illume:::ilm_rho_refit(flat, past, 2L, integer(0), numeric(0), list())
  expect_equal(o$par[[2]], illume:::ilm_rho_floor)
  ## near the edge, inside the floor, with nothing better within: unchanged
  near <- list(par = c(a = 1, rho_raw = 5), objective = 0)
  expect_identical(illume:::ilm_rho_refit(flat, near, 2L, integer(0), numeric(0), list()), near)
})

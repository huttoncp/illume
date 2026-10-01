## ---------------------------------------------------------------------------
## The designs a fitted model needs for new rows, and where each new row falls
## among a correlation over time's cells.
##
## Built by the same code predict() uses -- the fixed design with a smooth's
## null space, a smooth's penalised basis, a random slope's variables, the
## zero part's and the dispersion model's designs -- so a prediction assembled
## from these matrices is predict()'s prediction.
## ---------------------------------------------------------------------------

#' The design matrices of a fitted model, for new rows
#'
#' Everything needed to assemble a prediction for new rows from the fit's
#' parameters or from draws of them ([ilm_draws()]): the fixed design, each
#' random term's design and the group each row belongs to, each smooth's
#' penalised basis, the zero part's and the dispersion model's designs, and
#' each row's place among the cells of a correlation over time.
#'
#' @details
#' **Groups.** A row's group is matched to the fitted levels by label; a group
#' the fit has not seen has no code and is marked `new_group`. Without the
#' grouping column the rows belong to no group the fit knows -- a prediction
#' for the typical group or averaged over the population needs none -- so
#' every row is `new_group`, with `level` and `group` missing. Only a
#' prediction for a particular group needs its column.
#'
#' **Over time.** `ar` places each row among its group's cells, as numbered by
#' [ilm_cells()]: `cell` if the row's time is one of them, and the nearest
#' cells before and after, `prev_cell` and `next_cell`, with the time from
#' each, `dt_prev` and `dt_next`, in the units of the fit's time. A time past
#' a group's last cell has a `prev_cell` and no `next_cell` -- the forecast
#' case. An AR(1) term's grid is fixed, so a time that falls between its steps
#' is an error rather than rounded to one. The time and the group are read
#' from the columns the term was built from, `ilm_rw1(~ time | group)` and
#' the like; for a term built from vectors, pass them as `time` and `group`.
#' For a term built by name, rows without its columns are placed nowhere:
#' every cell and gap is missing and every row `new_group`. A term built from
#' vectors has no columns to fall back on, so there, rows without `time` and
#' `group` are an error -- a forgotten argument would otherwise predict each
#' row as a new series.
#'
#' @param object A fitted `"ilm_model"`.
#' @param newdata A data frame of new rows, with the columns the model uses.
#' @param time,group For a correlation over time built from vectors, the new
#'   rows' times and groups. Not needed when it was built by name.
#' @return A list with
#'   \describe{
#'     \item{`X`}{the fixed design, columns as `colnames(object$X)`.}
#'     \item{`re`}{one element per grouping term: `Z`, its design (the
#'       intercept and slope columns, named by dimension as [ilm_ranef()]
#'       and [ilm_draws()]' map name them, `"(Intercept)"` for a random
#'       intercept alone); `level`, each row's group label;
#'       `group`, its position among the fitted levels, `NA` for a new group;
#'       `new_group`; and `factor`, the grouping variable.}
#'     \item{`smooth`}{one penalised basis per smooth term.}
#'     \item{`zi`, `disp`}{the zero part's and the dispersion model's designs,
#'       when the model has them. Their columns are named as the coefficients
#'       they multiply are, in `coef(object, full = TRUE)` and in
#'       [ilm_draws()]' map: `"zi:(Intercept)"`, `"disp:x"`.}
#'     \item{`ar`}{with a correlation over time, a data frame with a row per
#'       new row: `group`, `time`, `cell`, `prev_cell`, `next_cell`,
#'       `dt_prev`, `dt_next` and `new_group`.}
#'     \item{`offset`}{with an offset in the formula, its value at each new
#'       row, made from `newdata`'s own columns; `NULL` otherwise. It adds to
#'       `X %*% beta` with a coefficient of one.}
#'   }
#' @seealso [ilm_draws()], [ilm_cells()], [ilm_ranef()].
#' @examples
#' set.seed(1)
#' d <- data.frame(id = factor(rep(c("a", "b", "c"), each = 10)), t = rep(1:10, 3))
#' d$y <- rnorm(30) + rep(cumsum(rnorm(10, 0, 0.5)), 3)
#' fit <- ilm_model(y ~ 1, data = d, family = "gaussian",
#'                  ar = ilm_rw1(~ t | id), verbose = FALSE)
#' nd <- data.frame(t = c(5, 12, 3), id = c("a", "a", "z"))
#' ilm_matrices(fit, nd)$ar
#' @export
ilm_matrices <- function(object, newdata, time = NULL, group = NULL) {
  if (!inherits(object, "ilm_model"))
    stop("`object` must be a fitted ilm_model, not ", class(object)[1],
         call. = FALSE)
  if (!is.data.frame(newdata))
    stop("`newdata` must be a data frame", call. = FALSE)
  if (is.null(object$terms))
    stop("ilm_matrices() needs a model fitted with a formula, whose terms say ",
         "how to build the new rows' designs", call. = FALSE)
  nd <- ilm_newX(object, newdata)
  out <- list(X = nd$X, re = list(), smooth = list(), zi = NULL, disp = NULL,
              ar = NULL, offset = nd$offset)
  gk <- which(vapply(object$re, function(e) !identical(e$kind, "basis"), TRUE))
  env <- environment(object$formula)
  if (is.null(env)) env <- parent.frame()
  for (k in seq_along(object$re)) {
    e <- object$re[[k]]; nm <- names(object$re)[k]
    if (identical(e$kind, "basis")) {
      out$smooth[[nm]] <- ilm_basis_new(object, nd, k)
      next
    }
    Z <- ilm_re_design(object, k, newdata, e$d)
    if (is.null(Z))
      stop("the random term '", nm, "' varies with a column `newdata` does not ",
           "have", call. = FALSE)
    ## columns named as ilm_ranef() and ilm_draws()' map name the dimensions,
    ## a random intercept alone included: an unnamed column there was read by
    ## code built on it as no column at all
    colnames(Z) <- if (!is.null(colnames(e$Z))) colnames(e$Z)
      else if (e$d == 1L) "(Intercept)" else paste0("z", seq_len(e$d))
    b <- object$bars[[match(k, gk)]]
    g <- tryCatch(eval(b[[3L]], newdata, env), error = function(err) NULL)
    ## no grouping column: the rows belong to no group the fit knows -- the
    ## typical group's or the population's prediction, which needs no unit --
    ## so every row is a new group, without a level or a code
    if (is.null(g)) g <- rep(NA_character_, nrow(newdata))
    if (length(g) != nrow(newdata))
      stop("the grouping variable of the random term '", nm, "' has ",
           length(g), " values for ", nrow(newdata), " rows of `newdata`",
           call. = FALSE)
    lev <- as.character(g)
    code <- if (length(e$levels) == e$nl) match(lev, e$levels)
            else rep(NA_integer_, length(lev))
    out$re[[nm]] <- list(Z = Z, level = lev, group = code,
                         new_group = is.na(code),
                         factor = if (!is.null(e$factor)) e$factor else nm)
  }
  ## named as their coefficients are, in coef(object, full = TRUE) and in
  ## ilm_draws()' map, so a column finds its coefficient by name, as X's do
  if (!is.null(object$Zzi)) {
    out$zi <- ilm_zi_design(object$zi_formula, newdata, colnames(object$Zzi))
    colnames(out$zi) <- paste0("zi:", colnames(object$Zzi))
  }
  if (!is.null(object$disp_formula)) {
    out$disp <- ilm_disp_design(object, newdata)
    if (!is.null(out$disp))
      colnames(out$disp) <- paste0("disp:", colnames(object$Zd))
  }
  if (!is.null(object$ar)) {
    v <- object$ar$vars
    if (is.null(time) || is.null(group)) {
      ## a term built from vectors has no columns to read, so rows without
      ## times and groups are almost always a forgotten argument -- and placed
      ## nowhere, a prediction for a known series would come out as a new one
      if (is.null(v))
        stop("the correlation over time was built from vectors, so the new ",
             "rows' times and groups have to be given: pass `time` and ",
             "`group`, or build the term by name, ilm_rw1(~ time | group)",
             call. = FALSE)
      if (all(v %in% names(newdata))) {
        time <- newdata[[v[["time"]]]]; group <- newdata[[v[["group"]]]]
      }
    }
    ## a term built by name, and new rows without its columns: no place among
    ## the cells -- a prediction for no particular group, which needs none
    out$ar <- if (is.null(time) || is.null(group)) ilm_ar_unplaced(nrow(newdata))
      else {
        if (length(time) != nrow(newdata) || length(group) != nrow(newdata))
          stop("`time` and `group` need one value per row of `newdata`",
               call. = FALSE)
        ilm_ar_place(object, time, group)
      }
  }
  out
}

## Rows placed nowhere among the cells: no group, no time, every one new.
#' @keywords internal
#' @noRd
ilm_ar_unplaced <- function(n)
  data.frame(group = rep(NA_character_, n), time = rep(NA_real_, n),
             cell = rep(NA_integer_, n), prev_cell = rep(NA_integer_, n),
             next_cell = rep(NA_integer_, n), dt_prev = rep(NA_real_, n),
             dt_next = rep(NA_real_, n), new_group = rep(TRUE, n),
             stringsAsFactors = FALSE)

## Where each new row falls among its group's cells.
#' @keywords internal
#' @noRd
ilm_ar_place <- function(object, time, group) {
  ar <- object$ar; cl <- ilm_cells(object)
  tn <- suppressWarnings(as.numeric(time))
  if (anyNA(tn))
    stop("the new rows' times must be numbers, dates or date-times, with ",
         "none missing", call. = FALSE)
  if (identical(ar$type, "ar1")) {
    k <- (tn - ar$origin) / ar$step
    off <- abs(k - round(k)) > 1e-6
    if (any(off))
      stop("time ", format(time[which(off)[1L]]), " falls between the AR(1) ",
           "grid's steps, which are ", signif(ar$step, 6), " apart from ",
           format(ilm_cor_as_time(ar$origin, ar)), ". An AR(1) latent exists ",
           "only on the grid; ilm_car1() takes any time.", call. = FALSE)
  }
  gi <- match(as.character(group), ar$glev)
  ct <- as.numeric(cl$time); cg <- as.integer(cl$group)
  n <- length(tn)
  cell <- prv <- nxt <- rep(NA_integer_, n)
  dtp <- dtn <- rep(NA_real_, n)
  for (g in unique(gi[!is.na(gi)])) {
    idx <- which(cg == g); tg <- ct[idx]
    tol <- ilm_time_tol(tg)
    for (i in which(gi == g)) {
      at <- which(abs(tg - tn[i]) <= tol)
      if (length(at)) cell[i] <- idx[at[1L]]
      b <- which(tg < tn[i] - tol); a <- which(tg > tn[i] + tol)
      if (length(b)) { j <- max(b); prv[i] <- idx[j]; dtp[i] <- tn[i] - tg[j] }
      if (length(a)) { j <- min(a); nxt[i] <- idx[j]; dtn[i] <- tg[j] - tn[i] }
    }
  }
  out <- data.frame(group = as.character(group), cell = cell, prev_cell = prv,
                    next_cell = nxt, dt_prev = dtp, dt_next = dtn,
                    new_group = is.na(gi), stringsAsFactors = FALSE)
  out$time <- time
  out[, c("group", "time", "cell", "prev_cell", "next_cell", "dt_prev",
          "dt_next", "new_group")]
}

## How near a new row's time must be to one of its group's fitted times to
## be placed at it: a millionth of the smallest gap between the group's
## distinct fitted times, so that matching never merges two of them, however
## irregular the series. A tolerance scaled by the times' size did: date-times
## are seconds since 1970, about 1.7e9, so 1e-9 of that is 1.7 seconds, and a
## row one second after a fitted time was placed at it. With a single fitted
## time there is no gap to scale by, and the tolerance is an absolute 1e-9.
#' @keywords internal
#' @noRd
ilm_time_tol <- function(tg) {
  u <- sort(unique(tg))
  if (length(u) < 2L) 1e-9 else 1e-6 * min(diff(u))
}

## The dispersion model's design for new rows. `mu` is reserved -- a power of
## the fitted mean, not a column -- so it is left to whoever forms the mean.
#' @keywords internal
#' @noRd
ilm_disp_design <- function(object, newdata) {
  f <- object$disp_formula
  tl <- setdiff(attr(stats::terms(f), "term.labels"), "mu")
  f2 <- if (length(tl)) stats::reformulate(tl, env = environment(f))
        else stats::as.formula("~ 1", environment(f))
  if (is.null(object$Zd)) return(NULL)
  ilm_zi_design(f2, newdata, colnames(object$Zd))
}

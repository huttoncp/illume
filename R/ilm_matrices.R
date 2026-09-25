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
#' the fit has not seen has no code and is marked `new_group`.
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
#'
#' @param object A fitted `"ilm_model"`.
#' @param newdata A data frame of new rows, with the columns the model uses.
#' @param time,group For a correlation over time built from vectors, the new
#'   rows' times and groups. Not needed when it was built by name.
#' @return A list with
#'   \describe{
#'     \item{`X`}{the fixed design, columns as `colnames(object$X)`.}
#'     \item{`re`}{one element per grouping term: `Z`, its design (the
#'       intercept and slope columns); `level`, each row's group label;
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
              ar = NULL)
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
    b <- object$bars[[match(k, gk)]]
    g <- tryCatch(eval(b[[3L]], newdata, env), error = function(err) NULL)
    if (is.null(g) || length(g) != nrow(newdata))
      stop("`newdata` does not have the grouping variable of the random term '",
           nm, "'", call. = FALSE)
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
      if (is.null(v))
        stop("the correlation over time was built from vectors, so the new ",
             "rows' times and groups have to be given: pass `time` and ",
             "`group`, or build the term by name, ilm_rw1(~ time | group)",
             call. = FALSE)
      miss <- setdiff(v, names(newdata))
      if (length(miss))
        stop("`newdata` needs ", paste(sQuote(miss), collapse = " and "),
             " to place its rows among the cells", call. = FALSE)
      time <- newdata[[v[["time"]]]]; group <- newdata[[v[["group"]]]]
    }
    if (length(time) != nrow(newdata) || length(group) != nrow(newdata))
      stop("`time` and `group` need one value per row of `newdata`",
           call. = FALSE)
    out$ar <- ilm_ar_place(object, time, group)
  }
  out
}

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
    tol <- 1e-9 * max(1, abs(tg))
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

## ---------------------------------------------------------------------------
## Helpers the family's packages share.
##
## This file is the same, byte for byte, in each of them; illumex's
## dev/check_shared_helpers.R compares it with illume's main. A change is
## made in all of them at once, in the same bytes.
## ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

## Prose wrapped to the console, with an indent.
#' @keywords internal
#' @noRd
ilm_wrap <- function(x, width = 76L, indent = "") {
  paste0(indent, strwrap(x, width = width - nchar(indent)), collapse = "\n")
}

## Names as they must be written in R code: backticked when not syntactic,
## untouched when they are. deparse() of a symbol knows the rules, reserved
## words included, which a comparison against make.names() gets subtly wrong.
#' @keywords internal
#' @noRd
ilm_bq <- function(x) {
  if (!length(x)) return(character(0))
  vapply(as.character(x), function(v) deparse(as.name(v), backtick = TRUE), "",
         USE.NAMES = FALSE)
}

## "a", "a and b", "a, b and c"
#' @keywords internal
#' @noRd
ilm_and <- function(x) {
  if (length(x) <= 1L) return(paste(x, collapse = ""))
  paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
}

## ---------------------------------------------------------------------------
## Progress reporting.
##
## Measured before being added: against a realistic loop -- 2000 bootstrap
## replicates on 20,000 rows -- a bar updated every 1%, or even every single
## iteration, made no difference beyond the run-to-run noise. The overhead
## only shows against a body so fast that a bar would be pointless anyway.
##
## So the bar is free wherever it is worth having, and the 1% step is there to
## avoid writing to a slow console rather than to save time.
##
## The default is interactive(): visible when a person is watching, silent in
## scripts, tests and knitr. A bar written into a vignette or a test log is
## noise, and would break every expect_silent() in the suite.
## ---------------------------------------------------------------------------

## A bar, or a silent stand-in with the same shape so callers need no branch.
#' @keywords internal
#' @noRd
ilm_progress <- function(n, progress = NULL, label = NULL) {
  on <- isTRUE(progress) ||
    (is.null(progress) && interactive() && n > 1L)
  if (!on || !is.finite(n) || n < 1L)
    return(list(tick = function(i) invisible(NULL),
                done = function() invisible(NULL)))
  if (!is.null(label)) message(label)
  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  step <- max(1L, as.integer(n) %/% 100L)
  list(
    tick = function(i) {
      if (i %% step == 0L || i == n) utils::setTxtProgressBar(pb, i)
      invisible(NULL)
    },
    done = function() { utils::setTxtProgressBar(pb, n); close(pb); invisible(NULL) })
}

## ---------------------------------------------------------------------------
## The random-number stream is the user's.
##
## A function given a seed sets it for its own draws -- so its result is the
## same every time -- and used to leave the stream there, so that a user's
## next rnorm() after ilm_sim(), ilm_cluster() or ilm_anomaly() came out the
## same whatever came before it. That is a side effect nobody asked for, and
## one that makes a simulation of the user's own quietly less random.
##
## So a seeded function puts the stream back as it leaves: as it was, or
## absent if it was absent (a fresh session has no .Random.seed until the
## first draw, and should not be given one). The function's own draws are
## untouched -- set.seed() stays where it was -- so every seeded result is
## what it was. With seed = NULL nothing is set and nothing is restored: the
## function draws from the user's stream, as any R code does.
## ---------------------------------------------------------------------------

## Register, in the CALLER's frame, a restore of the global random stream to
## run when the caller exits -- normally or by an error. The deferred call is
## base R's on.exit(), reached from the caller's frame the way withr::defer()
## reaches it, with add = TRUE so a handler the caller registers later cannot
## replace it, as long as that one adds too.
#' @keywords internal
#' @noRd
ilm_rng_restore <- function(seed, envir = parent.frame()) {
  if (is.null(seed)) return(invisible(FALSE))
  genv <- globalenv()
  had <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = genv, inherits = FALSE)
  restore <- function() {
    if (had) assign(".Random.seed", old, envir = genv)
    else if (exists(".Random.seed", envir = genv, inherits = FALSE))
      rm(".Random.seed", envir = genv)
  }
  thunk <- as.call(list(function() restore()))
  do.call(base::on.exit, list(thunk, TRUE, FALSE), envir = envir)
  invisible(TRUE)
}

## ---- plotting characters by name -------------------------------------------
##
## Base R's `pch` is 26 integers nobody remembers. Nothing in the argument says
## that 16 is a filled circle, and the difference between 16, 19, 20 and 21 is
## not guessable from the numbers -- so a plot gets whichever code the author
## happened to recall, and a reader comparing two plots cannot tell whether a
## difference in the markers was meant.
##
## Names are also checkable, which numbers are not: a pch outside 0:25 is
## accepted by graphics and quietly draws nothing, while a name that is not in
## the table stops here and lists the ones that are.
##
## A SINGLE character is left alone, because base R draws it literally --
## pch = "x" means the letter x, and translating it would silently turn the
## plot into crosses.

#' Plotting characters, by name
#'
#' The lookup behind the `pch` argument of the plotting functions. Numbers pass
#' through after a range check, a single character passes through (base R draws
#' it literally), and a name becomes the number base R wants.
#'
#' @param x Numeric `pch` codes, single characters, or names such as
#'   `"filled circle"`. Case, spaces, underscores, hyphens and dots are all
#'   ignored, so `"filled_circle"` and `"Filled Circle"` are the same thing.
#' @return An integer `pch` vector, or the input unchanged where it was already
#'   numeric or a single character.
#' @keywords internal
#' @noRd
ilm_pch <- function(x) {
  if (is.null(x) || !length(x)) return(x)
  if (is.numeric(x)) {
    bad <- x[is.finite(x) & (x < 0 | x > 25 | x != as.integer(x))]
    if (length(bad))
      stop("`pch` codes run from 0 to 25; got ",
           paste(unique(bad), collapse = ", "),
           ". Names work too, such as \"filled circle\".", call. = FALSE)
    return(x)
  }
  if (!is.character(x)) return(x)
  tab <- ilm_pch_table()
  key <- tolower(gsub("[ _.-]+", "", trimws(x)))
  out <- vector("list", length(x))
  for (i in seq_along(x)) {
    if (is.na(x[i])) { out[[i]] <- NA_integer_; next }
    if (nchar(x[i]) == 1L) { out[[i]] <- x[i]; next }   # a literal glyph
    j <- match(key[i], names(tab))
    if (is.na(j))
      stop("'", x[i], "' is not a plotting character name. ",
           "The names are: ", paste(ilm_pch_names(), collapse = ", "), ".",
           call. = FALSE)
    out[[i]] <- tab[[j]]
  }
  ## a mix of glyphs and codes has to stay character, since that is the only
  ## vector that can carry both
  if (any(vapply(out, is.character, TRUE))) as.character(unlist(out))
  else as.integer(unlist(out))
}

## The names, in code order, primary name first. Aliases follow it, because a
## person reaching for "solid circle" should not have to find out that the package
## calls it something else.
#' @keywords internal
#' @noRd
ilm_pch_spec <- function() {
  list(
    c("0",  "open square", "square", "hollow square", "empty square"),
    c("1",  "open circle", "circle", "hollow circle", "empty circle"),
    c("2",  "open triangle", "triangle", "triangle up", "hollow triangle"),
    c("3",  "plus"),
    c("4",  "cross", "times"),
    c("5",  "open diamond", "diamond", "hollow diamond"),
    c("6",  "open triangle down", "triangle down", "down triangle"),
    c("7",  "square cross", "crossed square"),
    c("8",  "star", "asterisk"),
    c("9",  "diamond plus"),
    c("10", "circle plus"),
    c("11", "star of david", "double triangle"),
    c("12", "square plus"),
    c("13", "circle cross"),
    c("14", "square triangle"),
    c("15", "filled square", "solid square"),
    c("16", "filled circle", "solid circle", "point", "dot"),
    c("17", "filled triangle", "solid triangle"),
    c("18", "filled diamond", "solid diamond"),
    c("19", "large filled circle", "bold circle"),
    c("20", "small filled circle", "bullet", "small dot"),
    c("21", "circle fill", "filled circle outline", "bg circle"),
    c("22", "square fill", "filled square outline", "bg square"),
    c("23", "diamond fill", "filled diamond outline", "bg diamond"),
    c("24", "triangle fill", "filled triangle outline", "bg triangle"),
    c("25", "triangle down fill", "filled triangle down"))
}

#' @keywords internal
#' @noRd
ilm_pch_table <- function() {
  out <- list()
  for (s in ilm_pch_spec()) {
    code <- as.integer(s[1])
    for (nm in s[-1]) out[[gsub("[ _.-]+", "", nm)]] <- code
  }
  out
}

## The primary name of each code, spelled the way a person would write it --
## the lookup keys have had their spaces stripped and would read badly here.
#' @keywords internal
#' @noRd
ilm_pch_names <- function() {
  vapply(ilm_pch_spec(), function(s) sprintf("%s (%s)", s[2], s[1]), "")
}

## ---------------------------------------------------------------------------
## The display rules every printed number follows, so that each package
## shows the same digits for the same value:
##   1. rounding, half away from zero on the value's 15-figure decimal form;
##   2. a display-rule formatter returning {text, rule, digits}.
## Names carry an `ilm_disp_` prefix to keep them apart from each package's
## older helpers of the same idea (ilm_fmt_pct(): in illume a phrase with
## "under 1%"/"over 99%" ends, in illumex a plain percentage). The prefix
## alone does not prevent a clash: illume's dispersion helpers are also
## ilm_disp_*, so no shared name may be ilm_disp_design, _flat, _limit,
## _rows, _vec or _words, and a name defined here must be defined nowhere
## else in a package's R/ -- this file collates late and would silently
## replace it. Each package's tests check that (test-shared-helpers.R).
## ---------------------------------------------------------------------------

## Rounding half away from zero, on the value's 15-significant-figure decimal
## form, in whole-number arithmetic. A double holds 15 significant figures
## faithfully, so its 15-figure form is the decimal it stands for: 2.675 is
## 2.67500000000000 there, though stored as 2.67499999999999982, and 100 x
## 0.2725 is 27.2500000000000 though computed as 27.250000000000004. R's
## round() and sprintf() follow the binary value instead, and round an exact
## tie to even. No noise allowance is added: any allowance either misses
## noise or reaches a digit a long number does have (95074673.3685499 to
## four decimals is 95074673.3685, a genuine 15-figure value 0.001 units
## below the tie). The cases at both ends are in
## tests/testthat/fixtures/format_cases.csv and go with any change.

## powers of ten as literals, each exact in a double (1e22 is the largest
## that is), so scaling by one is a single correctly rounded step
ILM_POW10 <- c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11,
               1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22)

## a x 10^k for whole k, by the table: one step within 10^22, and further
## steps of at most 10^22 beyond it (values beyond 1e30 or under 1e-30,
## where a step's rounding is far below any figure shown)
#' @keywords internal
#' @noRd
ilm_disp_scale <- function(a, k) {
  s <- a
  repeat {
    st <- pmax(pmin(k, 22), -22)
    st[is.na(st)] <- 0
    if (all(st == 0)) break
    s <- ifelse(st >= 0, s * ILM_POW10[abs(st) + 1], s / ILM_POW10[abs(st) + 1])
    k <- k - st
  }
  s
}

## The error of the product p = u x v, exactly, by Dekker's split: two
## doubles of 26 bits each multiply without rounding, so no fused
## multiply-add is needed and every platform gets the same answer
#' @keywords internal
#' @noRd
ilm_disp_prod_err <- function(u, v, p) {
  cu <- 134217729 * u; uh <- cu - (cu - u); ul <- u - uh
  cv <- 134217729 * v; vh <- cv - (cv - v); vl <- v - vh
  ((uh * vh - p) + uh * vl + ul * vh) + ul * vl
}

## a x 10^k rounded half up to a whole number, deciding by the exact value.
## The scaling is carried as a double and its exact error: each step by the
## table (at most 10^22, exact in a double) adds its own error, the
## product's by Dekker's split or the quotient's exact residual, so after
## any number of steps the whole number is decided by the true value. A
## product can round onto a half from below (8465491152.780715 x 10^5 is
## 846549115278071.49... but rounds to ...071.5), and steps beyond 10^22 can
## each carry the last figure off by one; neither can now.
#' @keywords internal
#' @noRd
ilm_disp_mantissa <- function(a, k) {
  hi <- a
  lo <- rep(0, length(a))
  repeat {
    st <- pmax(pmin(k, 22), -22)
    st[is.na(st)] <- 0
    if (all(st == 0)) break
    P <- ILM_POW10[abs(st) + 1]
    up <- which(st > 0)
    dn <- which(st < 0)
    if (length(up)) {
      h <- hi[up] * P[up]
      lo[up] <- ilm_disp_prod_err(hi[up], P[up], h) + lo[up] * P[up]
      hi[up] <- h
    }
    if (length(dn)) {
      h <- hi[dn] / P[dn]
      ## beside the largest double the split would overflow, so the residual
      ## is worked at a power of two below, which scales exactly
      sc <- ifelse(abs(hi[dn]) > 1e290, 2^-60, 1)
      hs <- h * sc
      pr <- hs * P[dn]
      rho <- ((hi[dn] * sc - pr) - ilm_disp_prod_err(hs, P[dn], pr)) / sc
      lo[dn] <- (rho + lo[dn]) / P[dn]
      hi[dn] <- h
    }
    k <- k - st
  }
  ## hi + 0.5 is exact below 2^52; the error then moves the whole number
  ## only across the boundary it straddles
  M <- floor(hi + 0.5)
  f <- (hi + 0.5) - M + lo
  M[f < 0] <- M[f < 0] - 1
  M[f >= 1] <- M[f >= 1] + 1
  M
}

## The core: each value's 15-figure mantissa M (10^14 <= M < 10^15, so
## value = M x 10^(e - 14)), rounded half away from zero to `d` decimals.
## Returns the digits kept as a whole number q, the decimals they stand at
## (dq, which is d, or fewer when every figure is kept), the exponent e, and
## the sign. NaN, Inf and NA give NA; 0 gives q = 0.
#' @keywords internal
#' @noRd
ilm_disp_digits <- function(x, d) {
  n <- length(x)
  d <- rep_len(as.numeric(d), n)
  q <- rep(NA_real_, n)
  dq <- d
  e <- rep(0, n)
  ok <- is.finite(x)
  q[ok & x == 0] <- 0
  i <- which(ok & x != 0)
  if (length(i)) {
    a <- abs(x[i])
    ei <- floor(log10(a))
    M <- ilm_disp_mantissa(a, 14 - ei)
    ## log10() can be one off beside a power of ten
    for (step in 1:2) {
      hi <- M >= 1e15
      lo <- M < 1e14
      if (!any(hi | lo)) break
      ei[hi] <- ei[hi] + 1
      ei[lo] <- ei[lo] - 1
      fix <- hi | lo
      M[fix] <- ilm_disp_mantissa(a[fix], 14 - ei[fix])
    }
    j <- (14 - ei) - d[i]
    Q <- M
    DQ <- 14 - ei                         # j <= 0: every figure is kept
    gone <- j > 15                        # under half a unit of the last decimal
    Q[gone] <- 0
    DQ[gone] <- d[i][gone]
    r <- which(j >= 1 & j <= 15)
    if (length(r)) {
      P <- ILM_POW10[j[r] + 1]
      rem <- M[r] %% P
      Q[r] <- (M[r] - rem) / P + (rem >= 5 * ILM_POW10[j[r]])
      DQ[r] <- d[i][r]
    }
    q[i] <- Q
    dq[i] <- DQ
    e[i] <- ei
  }
  list(q = q, dq = dq, e = e, neg = ok & x < 0 & !is.na(q) & q > 0)
}

## the number q x 10^-dq, with its sign
#' @keywords internal
#' @noRd
ilm_disp_value <- function(z) {
  v <- ilm_disp_scale(z$q, -z$dq)
  ifelse(z$neg, -v, v)
}

## To `decimals` decimals, half away from zero, by the 15-figure form. NaN
## and Inf pass through; -0 becomes 0. A number the rounding would carry
## past the largest double (1.7976931348623157e308 to 2 figures is 1.8e308)
## is returned as it is: no double holds the rounded value, and Inf would
## break whatever is computed from it. Its text, written from the digits,
## still reads the rounded figures. A subnormal (below 2.2e-308) holds fewer
## than 15 figures, so its rounding is as good as its bits allow.
#' @keywords internal
#' @noRd
ilm_disp_round <- function(x, decimals) {
  z <- ilm_disp_digits(x, decimals)
  out <- ilm_disp_value(z)
  ## every figure kept (a value of 15 or more figures before the point, or
  ## asked for at its full precision): the value is itself
  pass <- !is.finite(x) | (!is.na(z$q) & z$q != 0 & z$dq < decimals) |
    (is.finite(x) & !is.finite(out))
  out[pass] <- x[pass]
  out[out == 0] <- 0
  out
}

## To `digits` significant figures, by the same rounding: as digits and
## decimals (for the text), or as the number
#' @keywords internal
#' @noRd
ilm_disp_signif_digits <- function(x, digits) {
  z <- ilm_disp_digits(x, digits - 1 - ilm_disp_digits(x, 0)$e)
  ## a carry into the next power of ten (9.995 to 3 figures is 10.0): the
  ## figures stay `digits`
  up <- !is.na(z$q) & z$q >= ILM_POW10[digits + 1]
  z$q[up] <- z$q[up] / 10
  z$dq[up] <- z$dq[up] - 1
  ## zero has no significant figures to show: it is written "0"
  z$dq[!is.na(z$q) & z$q == 0] <- 0
  z
}

#' @keywords internal
#' @noRd
ilm_disp_signif <- function(x, digits) {
  out <- ilm_disp_value(ilm_disp_signif_digits(x, digits))
  ## past the largest double, the value as it is (see ilm_disp_round())
  pass <- !is.finite(x) | !is.finite(out)
  out[pass] <- x[pass]
  out
}

## The text of q x 10^-dq, written from q's digits (never by printing the
## number), with `d` decimals: zeros added where dq < d, or dropped from the
## end with `drop_zeros`. Thousands grouped by `big_mark`.
#' @keywords internal
#' @noRd
ilm_disp_text <- function(z, d, big_mark = ",", drop_zeros = FALSE) {
  n <- length(z$q)
  d <- rep_len(d, n)
  out <- rep(NA_character_, n)
  for (k in which(!is.na(z$q))) {
    s <- sprintf("%.0f", z$q[k])
    dk <- z$dq[k]
    if (dk < 0) {
      if (z$q[k] != 0) s <- paste0(s, strrep("0", -dk))
      dk <- 0
    }
    if (nchar(s) <= dk) s <- paste0(strrep("0", dk - nchar(s) + 1), s)
    int <- substr(s, 1, nchar(s) - dk)
    frac <- if (dk > 0) substr(s, nchar(s) - dk + 1, nchar(s)) else ""
    if (d[k] > dk) frac <- paste0(frac, strrep("0", d[k] - dk))
    if (drop_zeros) frac <- sub("0+$", "", frac)
    if (nzchar(big_mark) && nchar(int) > 3)
      int <- gsub("(\\d)(?=(\\d{3})+$)", paste0("\\1", big_mark), int, perl = TRUE)
    out[k] <- paste0(if (z$neg[k]) "-", int, if (nzchar(frac)) paste0(".", frac))
  }
  out
}

## A quantity as a statement shows it, and the rule it was shown by.
## Returns list(text, rule, digits): `text` a character vector the length of
## `x`, NA where `x` is NA; `rule` and `digits` the rule applied (`digits`
## NULL for "p"). Each value's text depends on that value alone.
##   "signif"  `digits` significant figures, trailing zeros kept (51.0, 10.0,
##             2.50; Craig's item 256), or dropped with `trailing_zeros =
##             FALSE` (51, 10, 2.5). An estimate that happens to be whole
##             still shows its figures (51.0, and 19234 is 19,200): a count,
##             or a data value that is a whole number, takes "fixed" with 0
##             digits instead, as its caller declares -- whole numbers are
##             not guessed at here
##   "fixed"   `digits` decimals, zeros kept; a count, or a whole data
##             value, is "fixed" with 0
##   "percent" a proportion times 100 at `digits` decimals, with "%"
##   "p"       Craig's rule: three significant figures down to 0.001, and
##             "< 0.001" below (the text of the number only; the sentence
##             supplies "p = " or "p ")
##   "clock"   an hour (0 to 23) as a time of day: "21:00" at the start of a
##             stretch, and with `end = TRUE` its last minute, "02:59"
##   "ordinal" a day of the month (1 to 31) in English, "1st", "22nd",
##             "28th"
## Thousands are grouped from 1,000 up ("1,230"), the family's rule. Only
## quantities come through here: a year, an identifier or a code that happens
## to be numeric is data, written as it is ("2024", never "2,024").
## Infinite values read "Inf" and "-Inf". The list carries `digits` for
## signif, fixed and percent, `end` for a clock that ends a stretch, and
## `trailing_zeros = FALSE` for a signif that drops them, the exception (kept
## is the rule, and says nothing).
## Every text is written from the rounded digits (ilm_disp_text()), never by
## printing the rounded number.
#' @keywords internal
#' @noRd
ilm_disp <- function(x, rule = c("signif", "fixed", "percent", "p", "clock", "ordinal"),
                     digits = NULL, big_mark = ",", end = FALSE, trailing_zeros = TRUE) {
  rule <- match.arg(rule)
  if (rule %in% c("signif", "fixed", "percent") &&
      (is.null(digits) || length(digits) != 1L || digits < 0))
    stop("ilm_disp(): the \"", rule, "\" rule needs `digits`", call. = FALSE)
  x <- as.numeric(x)
  fin <- is.finite(x)
  txt <- rep(NA_character_, length(x))
  txt[!is.na(x) & !fin] <- ifelse(x[!is.na(x) & !fin] > 0, "Inf", "-Inf")
  if (any(fin)) {
    v <- x[fin]
    txt[fin] <- switch(rule,
      signif = {
        z <- ilm_disp_signif_digits(v, digits)
        ilm_disp_text(z, pmax(0, z$dq), big_mark, drop_zeros = !isTRUE(trailing_zeros))
      },
      fixed = ilm_disp_text(ilm_disp_digits(v, digits), digits, big_mark),
      percent = paste0(ilm_disp_text(ilm_disp_digits(100 * v, digits), digits, big_mark), "%"),
      p = {
        z <- ilm_disp_signif_digits(v, 3)
        ifelse(v < 0.001, "< 0.001", ilm_disp_text(z, pmax(0, z$dq), "", drop_zeros = TRUE))
      },
      clock = sprintf(if (isTRUE(end)) "%02d:59" else "%02d:00", as.integer(round(v))),
      ordinal = {
        d <- as.integer(round(v))
        paste0(d, ifelse(d %% 100 %in% 11:13, "th",
                         c("th", "st", "nd", "rd", rep("th", 6))[d %% 10 + 1]))
      })
  }
  c(list(text = txt, rule = rule),
    if (rule %in% c("signif", "fixed", "percent")) list(digits = as.integer(digits)),
    if (rule == "clock" && isTRUE(end)) list(end = TRUE),
    if (rule == "signif" && !isTRUE(trailing_zeros)) list(trailing_zeros = FALSE))
}

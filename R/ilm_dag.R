## ---------------------------------------------------------------------------
## Causal graphs: the structure, the adjustment sets it licenses, and the
## conditional independencies it claims.
##
## The graph algorithms are written here rather than taken from `dagitty`
## because dagitty imports V8, and a JavaScript engine is a heavy thing to
## require of someone who wants to fit a regression. dagitty stays in Suggests
## and the test suite pins these functions to it, which is the same arrangement
## the rest of the package uses for lme4, nlme and survreg: implement it, then
## make an independent implementation agree.
##
## The scope is a DAG with directed edges and optionally unobserved nodes. A
## bidirected edge is the same claim as an unobserved common cause, so
## `X <-> Y` is read as exactly that rather than refused.
## ---------------------------------------------------------------------------

## ---- construction ----------------------------------------------------------

#' A causal graph
#'
#' Holds a directed acyclic graph, which variables are unobserved, and which
#' are the exposure and the outcome. This is the object [ilm_adjust_sets()] and
#' [ilm_dag_test()] work from.
#'
#' Three input forms, all giving the same object:
#'
#' * a dagitty-style string, with or without the `dag { }` wrapper;
#' * a data frame of edges with columns `from` and `to`;
#' * a `dagitty` object, if that package is installed.
#'
#' A bidirected edge `X <-> Y` states that something unobserved causes both, so
#' it is stored as an unobserved node with an arrow into each. That keeps one
#' representation -- a DAG, possibly with unobserved nodes -- rather than two.
#'
#' @param x A dagitty-style string, a data frame of edges, or a `dagitty`
#'   object.
#' @param exposure,outcome Variable names. Optional if the string or object
#'   already declares them; required otherwise before an adjustment set can be
#'   found.
#' @param latent Names of unobserved variables, added to any the input declares.
#' @return An object of class `"ilm_dag"`.
#' @seealso [ilm_adjust_sets()] for what the graph licenses, [ilm_dag_test()]
#'   for whether the data agree with it.
#' @examples
#' g <- ilm_dag("dag {
#'   smoking [exposure]
#'   cancer  [outcome]
#'   smoking -> tar -> cancer
#'   genotype -> smoking
#'   genotype -> cancer
#' }")
#' g
#' @export
ilm_dag <- function(x, exposure = NULL, outcome = NULL, latent = NULL) {
  g <- if (inherits(x, "ilm_dag")) x
       else if (inherits(x, "dagitty")) ilm_dag_from_dagitty(x)
       else if (is.data.frame(x)) ilm_dag_from_edges(x)
       else if (is.character(x) && length(x) == 1L) ilm_dag_parse(x)
       else stop("`x` must be a dagitty-style string, a data frame of edges ",
                 "with columns `from` and `to`, or a dagitty object; it is ",
                 class(x)[1], call. = FALSE)

  if (!is.null(exposure)) g$exposure <- as.character(exposure)
  if (!is.null(outcome))  g$outcome  <- as.character(outcome)
  if (!is.null(latent))   g$latent   <- union(g$latent, as.character(latent))

  miss <- setdiff(c(g$exposure, g$outcome, g$latent), g$nodes)
  if (length(miss))
    stop("named variable(s) are not in the graph: ", paste(miss, collapse = ", "),
         ". The graph has: ", paste(g$nodes, collapse = ", "), ".", call. = FALSE)
  ## a cycle makes ancestry ill-defined, so it is caught at construction rather
  ## than as a hang or a wrong answer later
  cyc <- ilm_dag_cycle(g)
  if (!is.null(cyc))
    stop("the graph has a cycle, so it is not a DAG: ",
         paste(c(cyc, cyc[1]), collapse = " -> "), call. = FALSE)
  g
}

#' @keywords internal
#' @noRd
ilm_dag_new <- function(nodes, from, to, latent = character(),
                        exposure = character(), outcome = character()) {
  structure(list(nodes = unique(nodes),
                 edges = data.frame(from = from, to = to,
                                    stringsAsFactors = FALSE),
                 latent = unique(latent),
                 exposure = exposure, outcome = outcome),
            class = "ilm_dag")
}

## Parse the dagitty text form. Supported: an optional `dag { }` wrapper, node
## declarations carrying [exposure], [outcome], [latent] or [unobserved], edges
## written -> or <-, chains of them, and ; or a newline between statements.
#' @keywords internal
#' @noRd
ilm_dag_parse <- function(txt) {
  s <- gsub("^\\s*(dag|digraph)?\\s*\\w*\\s*\\{", "", txt)
  s <- sub("\\}\\s*$", "", s)
  ## a comment runs to the end of its line, so strip before flattening
  s <- gsub("#[^\n]*", "", s)
  stm <- unlist(strsplit(s, "[;\n]"))
  stm <- trimws(stm); stm <- stm[nzchar(stm)]
  if (!length(stm))
    stop("the graph is empty: no nodes or edges were found in the string",
         call. = FALSE)

  nodes <- character(); from <- character(); to <- character()
  latent <- character(); expo <- character(); outc <- character()

  for (st in stm) {
    if (grepl("<->", st, fixed = TRUE)) {
      ## a bidirected edge is an unobserved common cause; name it so the user
      ## can see it in the adjustment sets and in any error message
      p <- trimws(strsplit(st, "<->", fixed = TRUE)[[1]])
      p <- ilm_dag_bare(p)
      if (length(p) != 2L)
        stop("cannot read the bidirected edge: ", st, call. = FALSE)
      u <- paste0("U_", p[1], "_", p[2])
      nodes <- c(nodes, p, u); latent <- c(latent, u)
      from <- c(from, u, u); to <- c(to, p[1], p[2])
      next
    }
    if (grepl("--", st, fixed = TRUE))
      stop("undirected edges (--) have no causal direction, so no adjustment ",
           "set follows from them: ", st, call. = FALSE)

    if (grepl("->|<-", st)) {
      ## split on arrows but keep them, so a chain a -> b <- c reads correctly
      tk <- regmatches(st, gregexpr("->|<-|[^\\s>< -]+", st, perl = TRUE))[[1]]
      tk <- trimws(tk); tk <- tk[nzchar(tk)]
      ops <- tk[tk %in% c("->", "<-")]
      vs  <- ilm_dag_bare(tk[!tk %in% c("->", "<-")])
      if (length(vs) != length(ops) + 1L)
        stop("cannot read the edge statement: ", st, call. = FALSE)
      nodes <- c(nodes, vs)
      for (k in seq_along(ops)) {
        if (ops[k] == "->") { from <- c(from, vs[k]); to <- c(to, vs[k + 1L]) }
        else                { from <- c(from, vs[k + 1L]); to <- c(to, vs[k]) }
      }
      next
    }

    ## otherwise a bare node declaration, possibly with attributes
    nm <- ilm_dag_bare(sub("\\[.*$", "", st))
    if (!length(nm) || !nzchar(nm[1])) next
    nodes <- c(nodes, nm[1])
    at <- regmatches(st, regexpr("\\[[^]]*\\]", st))
    if (length(at)) {
      at <- tolower(at)
      if (grepl("exposure", at, fixed = TRUE)) expo <- c(expo, nm[1])
      if (grepl("outcome", at, fixed = TRUE))  outc <- c(outc, nm[1])
      if (grepl("latent", at, fixed = TRUE) ||
          grepl("unobserved", at, fixed = TRUE)) latent <- c(latent, nm[1])
    }
  }
  ilm_dag_new(nodes, from, to, latent, unique(expo), unique(outc))
}

## strip quotes and surrounding space from variable names
#' @keywords internal
#' @noRd
ilm_dag_bare <- function(v) {
  v <- trimws(as.character(v))
  v <- gsub('^["\']|["\']$', "", v)
  v[nzchar(v)]
}

#' @keywords internal
#' @noRd
ilm_dag_from_edges <- function(d) {
  if (!all(c("from", "to") %in% names(d)))
    stop("an edge data frame needs columns `from` and `to`; it has: ",
         paste(names(d), collapse = ", "), call. = FALSE)
  f <- ilm_dag_bare(d$from); t <- ilm_dag_bare(d$to)
  if (length(f) != length(t) || !length(f))
    stop("`from` and `to` must be the same non-zero length", call. = FALSE)
  ilm_dag_new(c(f, t), f, t)
}

#' @keywords internal
#' @noRd
ilm_dag_from_dagitty <- function(x) {
  if (!requireNamespace("dagitty", quietly = TRUE))
    stop("the object is a dagitty graph but the dagitty package is not ",
         "installed", call. = FALSE)
  g <- ilm_dag_parse(paste(format(x), collapse = "\n"))
  g$exposure <- union(g$exposure, dagitty::exposures(x))
  g$outcome  <- union(g$outcome,  dagitty::outcomes(x))
  g$latent   <- union(g$latent,   dagitty::latents(x))
  g
}

#' @export
print.ilm_dag <- function(x, ...) {
  cat("<ilm_dag>", length(x$nodes), "variables,", nrow(x$edges), "edges\n")
  cat("  exposure: ", if (length(x$exposure)) paste(x$exposure, collapse = ", ")
                      else "(not set)", "\n", sep = "")
  cat("  outcome:  ", if (length(x$outcome)) paste(x$outcome, collapse = ", ")
                      else "(not set)", "\n", sep = "")
  if (length(x$latent))
    cat("  unobserved: ", paste(x$latent, collapse = ", "), "\n", sep = "")
  cat("  edges:\n")
  for (i in seq_len(nrow(x$edges)))
    cat("    ", x$edges$from[i], " -> ", x$edges$to[i], "\n", sep = "")
  invisible(x)
}

## ---- graph primitives ------------------------------------------------------

## The smallest set of OBSERVED variables that d-separates a from b, or NULL if
## none does.
##
## Building the set from the parents of a and b and then rejecting it when a
## parent is unobserved is wrong, and was: the pair is often separated by a set
## that does not contain that parent, most commonly by the empty set. Any
## minimal separator lies inside the ancestors of the two variables, so the
## search is over subsets of the observed ancestors, shortest first -- which
## also makes the answer a minimum rather than merely minimal, so the test it
## produces holds as few things fixed as possible.
#' @keywords internal
#' @noRd
ilm_dag_sep_obs <- function(g, a, b, obs, max_pool = 12L) {
  if (ilm_dsep(g, a, b, character())) return(character())
  pool <- setdiff(intersect(ilm_dag_ancestors(g, c(a, b)), obs), c(a, b))
  if (!length(pool)) return(NULL)
  if (length(pool) > max_pool) {
    ## too many to enumerate: fall back to the ancestor set itself, reduced
    ## greedily. Sound but possibly not minimum, and never silently so.
    if (!ilm_dsep(g, a, b, pool)) return(NULL)
    repeat {
      drop <- NULL
      for (v in pool) if (ilm_dsep(g, a, b, setdiff(pool, v))) { drop <- v; break }
      if (is.null(drop)) return(pool)
      pool <- setdiff(pool, drop)
    }
  }
  for (k in seq_len(length(pool))) {
    cmb <- utils::combn(pool, k, simplify = FALSE)
    for (z in cmb) if (ilm_dsep(g, a, b, z)) return(z)
  }
  NULL
}

#' @keywords internal
#' @noRd
ilm_dag_parents <- function(g, v) unique(g$edges$from[g$edges$to %in% v])

#' @keywords internal
#' @noRd
ilm_dag_children <- function(g, v) unique(g$edges$to[g$edges$from %in% v])

## reachability by repeated expansion; graphs here are small enough that the
## simple fixed point is not worth replacing
#' @keywords internal
#' @noRd
ilm_dag_reach <- function(g, v, step) {
  out <- unique(as.character(v)); repeat {
    nx <- union(out, step(g, out))
    if (length(nx) == length(out)) return(out)
    out <- nx
  }
}

#' @keywords internal
#' @noRd
ilm_dag_ancestors <- function(g, v) ilm_dag_reach(g, v, ilm_dag_parents)

#' @keywords internal
#' @noRd
ilm_dag_descendants <- function(g, v) ilm_dag_reach(g, v, ilm_dag_children)

## Returns a cycle as a character vector, or NULL. Kahn's algorithm: peel off
## nodes with no remaining parent, and whatever will not peel lies on a cycle.
#' @keywords internal
#' @noRd
ilm_dag_cycle <- function(g) {
  left <- g$nodes; e <- g$edges
  repeat {
    has_parent <- left %in% e$to
    if (all(has_parent)) break
    drop <- left[!has_parent]
    if (!length(drop)) break
    left <- setdiff(left, drop)
    e <- e[e$from %in% left & e$to %in% left, , drop = FALSE]
    if (!length(left)) return(NULL)
  }
  if (!length(left)) return(NULL)
  ## walk forward from any survivor until a node repeats
  path <- left[1]; cur <- left[1]
  repeat {
    nxt <- e$to[e$from == cur]
    nxt <- nxt[nxt %in% left]
    if (!length(nxt)) return(left[1])
    cur <- nxt[1]
    if (cur %in% path) return(path[which(path == cur):length(path)])
    path <- c(path, cur)
  }
}

#' Are two variables d-separated given a conditioning set?
#'
#' The graphical criterion behind every claim a DAG makes about the data. Tested
#' by the moralisation route: take the ancestral subgraph of everything
#' mentioned, marry the parents, drop the directions, remove the conditioning
#' set, and ask whether a path remains.
#'
#' @param g An [ilm_dag()].
#' @param x,y Variable names.
#' @param z Conditioning set; may be empty.
#' @return `TRUE` if `x` and `y` are d-separated given `z`.
#' @examples
#' g <- ilm_dag("dag { z -> x -> y ; z -> y }")
#' ilm_dsep(g, "x", "y", character())   # FALSE, x causes y
#' ilm_dsep(g, "x", "y", "z")           # still FALSE, the direct edge remains
#' @export
ilm_dsep <- function(g, x, y, z = character()) {
  g <- ilm_dag(g)
  z <- as.character(z)
  all_v <- c(x, y, z)
  miss <- setdiff(all_v, g$nodes)
  if (length(miss))
    stop("variable(s) not in the graph: ", paste(miss, collapse = ", "),
         call. = FALSE)
  if (length(intersect(c(x, y), z)))
    stop("`x` and `y` cannot also be in the conditioning set", call. = FALSE)

  keep <- ilm_dag_ancestors(g, all_v)
  e <- g$edges[g$edges$from %in% keep & g$edges$to %in% keep, , drop = FALSE]
  ## moralise: every pair of parents of a common child becomes adjacent
  adj <- list()
  push <- function(a, b) {
    adj[[a]] <<- union(adj[[a]], b); adj[[b]] <<- union(adj[[b]], a)
  }
  for (v in keep) adj[[v]] <- character()
  for (i in seq_len(nrow(e))) push(e$from[i], e$to[i])
  for (v in unique(e$to)) {
    pa <- e$from[e$to == v]
    if (length(pa) > 1L)
      for (i in 1:(length(pa) - 1L)) for (j in (i + 1L):length(pa))
        push(pa[i], pa[j])
  }
  ## remove the conditioning set and look for any remaining path
  open <- setdiff(keep, z)
  if (!(x %in% open) || !(y %in% open)) return(TRUE)
  seen <- x; stack <- x
  while (length(stack)) {
    cur <- stack[1]; stack <- stack[-1]
    nb <- intersect(adj[[cur]], open)
    nb <- setdiff(nb, seen)
    if (y %in% nb) return(FALSE)
    seen <- c(seen, nb); stack <- c(stack, nb)
  }
  TRUE
}

## ---- adjustment sets -------------------------------------------------------

## The back-door criterion (Pearl 1995). Z is admissible for the effect of X on
## Y when no member of Z is a descendant of X, and Z blocks every path into X:
## the second half is d-separation of X and Y in the graph with X's OUTGOING
## edges deleted, which is what turns a path question into a test we already
## have.
#' @keywords internal
#' @noRd
ilm_backdoor_ok <- function(g, x, y, z) {
  if (length(intersect(z, ilm_dag_descendants(g, x)))) return(FALSE)
  gx <- g
  gx$edges <- g$edges[!(g$edges$from == x), , drop = FALSE]
  ilm_dsep(gx, x, y, z)
}

#' Minimal adjustment sets for an exposure effect
#'
#' Every set of observed variables that, adjusted for, identifies the total
#' effect of `exposure` on `outcome`, with no redundant member. An empty result
#' is a finding rather than a failure: it says the effect is not identifiable
#' from the variables available, which is the most useful thing a DAG can tell
#' you before any model is fitted.
#'
#' Several minimal sets are common and worth having. Adjusting for each in turn
#' should give the same exposure effect if the graph is right, so the spread
#' across them is a sensitivity analysis that costs only compute --
#' [ilm_dag_model()] runs it automatically.
#'
#' This is the back-door criterion, which is sufficient but not quite complete:
#' the full adjustment criterion admits a small number of further sets
#' containing descendants of the exposure that lie off every causal path. Where
#' they differ, the sets returned here are a subset of the valid ones, never a
#' superset, so a set reported here is always admissible.
#'
#' @param g An [ilm_dag()].
#' @param exposure,outcome Variable names; taken from the graph if it declares
#'   them.
#' @param observed Variables available for adjustment. Defaults to every node
#'   the graph does not mark unobserved. Pass `names(data)` to restrict to what
#'   was actually measured.
#' @param max_size Largest set to look for. The default searches every subset
#'   of the candidate pool when that is affordable and warns when it is not.
#' @return A list of character vectors, one per minimal set, shortest first.
#'   `character(0)` as the only element means the empty set suffices -- no
#'   adjustment is needed. A zero-length list means no admissible set exists.
#' @references
#' Pearl, J. (1995). Causal diagrams for empirical research. Biometrika 82(4).
#' @examples
#' g <- ilm_dag("dag {
#'   smoking [exposure]
#'   cancer  [outcome]
#'   smoking -> cancer
#'   genotype -> smoking
#'   genotype -> cancer
#' }")
#' ilm_adjust_sets(g)
#' ## and with the confounder unmeasured, nothing identifies the effect
#' ilm_adjust_sets(g, observed = c("smoking", "cancer"))
#' @export
ilm_adjust_sets <- function(g, exposure = NULL, outcome = NULL,
                            observed = NULL, max_size = NULL) {
  g <- ilm_dag(g)
  x <- if (is.null(exposure)) g$exposure else as.character(exposure)
  y <- if (is.null(outcome))  g$outcome  else as.character(outcome)
  if (length(x) != 1L || length(y) != 1L)
    stop("one `exposure` and one `outcome` are needed; the graph declares ",
         if (!length(g$exposure)) "no exposure" else
           paste0("exposure ", paste(g$exposure, collapse = "/")),
         " and ",
         if (!length(g$outcome)) "no outcome" else
           paste0("outcome ", paste(g$outcome, collapse = "/")),
         ". Pass them to ilm_dag() or to this function.", call. = FALSE)
  miss <- setdiff(c(x, y), g$nodes)
  if (length(miss))
    stop("variable(s) not in the graph: ", paste(miss, collapse = ", "),
         call. = FALSE)
  if (identical(x, y))
    stop("`exposure` and `outcome` are the same variable", call. = FALSE)

  obs <- if (is.null(observed)) setdiff(g$nodes, g$latent)
         else intersect(as.character(observed), g$nodes)
  ## the exposure's own descendants can never adjust for confounding, and
  ## conditioning on them removes part of the effect being estimated
  pool <- setdiff(obs, c(x, y, ilm_dag_descendants(g, x)))

  np <- length(pool)
  cap <- if (!is.null(max_size)) min(max_size, np) else np
  ## full enumeration is 2^np d-separation tests; past about sixteen candidates
  ## that stops being instant, so the search is capped and says so
  if (is.null(max_size) && np > 16L) {
    cap <- 3L
    warning("the graph offers ", np, " candidate variables, so every subset ",
            "is not searchable; looking at sets of up to ", cap,
            " variables. Pass `max_size` to change that, or `observed` to ",
            "narrow the pool.", call. = FALSE)
  }

  found <- list()
  for (k in 0:cap) {
    cmb <- if (k == 0L) list(character(0))
           else utils::combn(pool, k, simplify = FALSE)
    for (z in cmb) {
      ## a superset of an admissible set is admissible but not minimal
      if (length(found) &&
          any(vapply(found, function(f) all(f %in% z), TRUE))) next
      if (ilm_backdoor_ok(g, x, y, z)) found[[length(found) + 1L]] <- z
    }
  }
  attr(found, "exposure") <- x
  attr(found, "outcome") <- y
  attr(found, "pool") <- pool
  attr(found, "unobserved") <- setdiff(g$nodes, obs)
  found
}

## ---- what the graph claims about the data ----------------------------------

#' Conditional independencies the graph implies
#'
#' One testable claim per missing edge: two variables with no arrow between them
#' are independent given some set of the others. If the graph is right the data
#' should agree; where they do not, either the graph is wrong or the data are
#' not what they are taken to be.
#'
#' The conditioning set is minimal -- no member can be dropped without the
#' claim failing. Conditioning on the parents of both variables always works and
#' is the usual textbook basis, but it is often far larger than it needs to be,
#' and every unnecessary covariate costs power in the test and invites a
#' modelling error of its own. Where the two variables are already independent
#' with nothing held fixed, the set is empty.
#'
#' @param g An [ilm_dag()].
#' @param observed Restrict to claims every variable of which was measured.
#' @return A data frame with `x`, `y` and a list column `given`.
#' @seealso [ilm_dag_test()], which tests these against data.
#' @examples
#' g <- ilm_dag("dag { z -> x -> y ; z -> w }")
#' ilm_dag_implied(g)
#' @export
ilm_dag_implied <- function(g, observed = NULL) {
  g <- ilm_dag(g)
  obs <- if (is.null(observed)) setdiff(g$nodes, g$latent)
         else intersect(as.character(observed), g$nodes)
  nd <- intersect(g$nodes, obs)
  out <- list()
  if (length(nd) >= 2L) {
    pr <- utils::combn(nd, 2L, simplify = FALSE)
    for (p in pr) {
      a <- p[1]; b <- p[2]
      adjacent <- any((g$edges$from == a & g$edges$to == b) |
                      (g$edges$from == b & g$edges$to == a))
      if (adjacent) next
      ## NULL means no set of measured variables makes the two independent, so
      ## the graph makes no claim here that the data could check
      z <- ilm_dag_sep_obs(g, a, b, obs)
      if (is.null(z)) next
      out[[length(out) + 1L]] <- list(x = a, y = b, given = z)
    }
  }
  if (!length(out))
    return(structure(data.frame(x = character(), y = character(),
                                stringsAsFactors = FALSE),
                     given = list()))
  d <- data.frame(x = vapply(out, `[[`, "", "x"),
                  y = vapply(out, `[[`, "", "y"),
                  stringsAsFactors = FALSE)
  d$given <- lapply(out, `[[`, "given")
  d
}

## ---- testing the graph against the data ------------------------------------

## What kind of variable is this, in the sense of which likelihood to test it
## under? Deliberately conservative: anything it cannot place confidently comes
## back "other" and the claim is skipped with a reason rather than tested under
## a guess.
#' @keywords internal
#' @noRd
ilm_var_kind <- function(v) {
  if (is.logical(v)) return("binary")
  if (is.factor(v) || is.character(v)) {
    k <- length(unique(stats::na.omit(as.character(v))))
    return(if (k == 2L) "binary" else if (k < 2L) "constant" else "nominal")
  }
  if (!is.numeric(v)) return("other")
  u <- unique(stats::na.omit(v))
  if (length(u) < 2L) return("constant")
  if (length(u) == 2L) return("binary")
  if (all(u >= 0) && all(abs(u - round(u)) < 1e-8)) return("count")
  "continuous"
}

#' @keywords internal
#' @noRd
ilm_ci_fit <- function(y, kind, rhs, d) {
  fam <- switch(kind, continuous = stats::gaussian(),
                binary = stats::binomial(), count = stats::poisson())
  f <- stats::reformulate(if (length(rhs)) rhs else "1", response = y)
  suppressWarnings(stats::glm(f, data = d, family = fam))
}

#' Test what the graph claims against the data
#'
#' Every missing edge in a DAG is a testable claim: those two variables should
#' be independent given the conditioning set. This fits that claim and reports
#' whether the data agree. It is the step that turns a DAG from an assumption
#' into something the data can argue with.
#'
#' A failed claim does not say which of the two is at fault. The graph may be
#' missing an arrow, or the data may not measure what the graph's node names
#' suppose -- a mismeasured covariate, a selected sample, a variable that means
#' something different from what it is called. Both are worth knowing before
#' any effect estimate is taken seriously.
#'
#' Each claim is tested by a likelihood-ratio test of the exposure term in a
#' generalised linear model for one of the pair given the other and the
#' conditioning set, with the family taken from the response's type. That is a
#' screen rather than a full model: it assumes the conditional mean is linear on
#' the link scale, so a curved dependence can pass. It will not test a claim
#' between two variables that are both multi-level factors, and says so in the
#' `note` column rather than quietly dropping it.
#'
#' **Two thresholds, not one.** With a few thousand rows, a correlation of 0.03
#' is significant and means nothing, so a claim counts as contradicted only when
#' it is both statistically significant after adjustment and larger than
#' `min_effect` in partial correlation. This is the stance [ilm_variogram()]
#' takes for the same reason.
#'
#' @param g An [ilm_dag()].
#' @param data A data frame holding the graph's observed variables.
#' @param min_effect Smallest partial correlation worth calling a contradiction.
#' @param alpha Level for the adjusted p-values.
#' @param adjust Multiplicity adjustment across claims, passed to
#'   [stats::p.adjust()]. A DAG can imply dozens of claims and testing them all
#'   at 0.05 guarantees false alarms, so the default is Holm.
#' @param verbose Print progress.
#' @return A data frame with one row per claim: `x`, `y`, `given`, `n`,
#'   `estimate` (partial correlation), `p_value`, `p_adj`, `verdict` and `note`.
#'   The overall verdict is attached as the `verdict` attribute.
#' @seealso [ilm_dag_implied()] for the claims themselves.
#' @examples
#' set.seed(1)
#' n <- 400
#' z <- rnorm(n); x <- 0.6 * z + rnorm(n); y <- 0.5 * x + 0.4 * z + rnorm(n)
#' d <- data.frame(x = x, y = y, z = z, w = rnorm(n))
#' ## w is unconnected, so every claim involving it should hold
#' g <- ilm_dag("dag { x [exposure] ; y [outcome] ; z -> x -> y ; z -> y }")
#' ilm_dag_test(g, cbind(d["x"], d["y"], d["z"]), verbose = FALSE)
#' @export
ilm_dag_test <- function(g, data, min_effect = 0.1, alpha = 0.05,
                         adjust = "holm", verbose = TRUE) {
  g <- ilm_dag(g)
  if (!is.data.frame(data))
    stop("`data` must be a data frame; it is ", class(data)[1], call. = FALSE)
  obs <- intersect(g$nodes, names(data))
  absent <- setdiff(setdiff(g$nodes, g$latent), names(data))
  if (length(absent) && verbose)
    message("  variables in the graph but not in the data, treated as ",
            "unobserved: ", paste(absent, collapse = ", "))

  im <- ilm_dag_implied(g, observed = obs)
  if (!nrow(im)) {
    if (verbose)
      message("  the graph implies no conditional independence that these ",
              "variables can test")
    out <- data.frame(x = character(), y = character(), n = integer(),
                      estimate = numeric(), p_value = numeric(),
                      p_adj = numeric(), verdict = character(),
                      note = character(), stringsAsFactors = FALSE)
    out$given <- list()
    return(structure(out, verdict = "UNTESTED"))
  }
  if (verbose)
    message("  testing ", nrow(im), " conditional independenc",
            if (nrow(im) == 1L) "y" else "ies", " implied by the graph")

  rows <- lapply(seq_len(nrow(im)), function(i) {
    a <- im$x[i]; b <- im$y[i]; z <- im$given[[i]]
    ka <- ilm_var_kind(data[[a]]); kb <- ilm_var_kind(data[[b]])
    ## the response has to be something a GLM can take; a nominal factor can
    ## be a predictor but not a response here
    resp <- if (ka %in% c("continuous", "binary", "count")) a
            else if (kb %in% c("continuous", "binary", "count")) b else NA
    note <- ""
    if (is.na(resp))
      note <- sprintf("not tested: %s is %s and %s is %s", a, ka, b, kb)
    pred <- if (identical(resp, a)) b else a

    d <- data[, unique(c(a, b, z)), drop = FALSE]
    d <- d[stats::complete.cases(d), , drop = FALSE]
    if (nzchar(note) || nrow(d) < 10L) {
      if (!nzchar(note)) note <- "not tested: fewer than 10 complete rows"
      return(data.frame(x = a, y = b, n = nrow(d), estimate = NA_real_,
                        p_value = NA_real_, verdict = "UNTESTED", note = note,
                        stringsAsFactors = FALSE))
    }
    kr <- ilm_var_kind(d[[resp]])
    fit0 <- try(ilm_ci_fit(resp, kr, z, d), silent = TRUE)
    fit1 <- try(ilm_ci_fit(resp, kr, c(pred, z), d), silent = TRUE)
    if (inherits(fit0, "try-error") || inherits(fit1, "try-error"))
      return(data.frame(x = a, y = b, n = nrow(d), estimate = NA_real_,
                        p_value = NA_real_, verdict = "UNTESTED",
                        note = "not tested: the model would not fit",
                        stringsAsFactors = FALSE))
    an <- suppressWarnings(stats::anova(fit0, fit1,
            test = if (kr == "continuous") "F" else "Chisq"))
    pcol <- grep("^Pr\\(", names(an), value = TRUE)
    pv <- if (length(pcol)) an[[pcol[1]]][2] else NA_real_

    ## a magnitude on one scale whatever the family: the correlation left
    ## between the two once the conditioning set is taken out of both
    pc <- tryCatch({
      ra <- stats::residuals(ilm_ci_fit(resp, kr, z, d), type = "pearson")
      kp <- ilm_var_kind(d[[pred]])
      rb <- if (kp %in% c("continuous", "binary", "count"))
              stats::residuals(ilm_ci_fit(pred, kp, z, d), type = "pearson")
            else {
              mm <- stats::model.matrix(stats::reformulate(pred, response = NULL),
                                        d)[, -1, drop = FALSE]
              if (!length(z)) mm[, 1] else
                stats::residuals(stats::lm(mm[, 1] ~ ., data = d[z]))
            }
      unname(stats::cor(ra, rb, use = "complete.obs"))
    }, error = function(e) NA_real_)

    data.frame(x = a, y = b, n = nrow(d), estimate = pc, p_value = pv,
               verdict = NA_character_, note = "", stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$p_adj <- stats::p.adjust(out$p_value, method = adjust)
  ## contradicted needs BOTH: significant after adjustment, and big enough to
  ## act on. Either alone is not a finding.
  sig <- !is.na(out$p_adj) & out$p_adj < alpha
  big <- !is.na(out$estimate) & abs(out$estimate) >= min_effect
  out$verdict <- ifelse(out$verdict %in% "UNTESTED", "UNTESTED",
                 ifelse(sig & big, "FAIL",
                 ifelse(sig & !big, "WARN", "OK")))
  out$given <- im$given
  out <- out[, c("x", "y", "given", "n", "estimate", "p_value", "p_adj",
                 "verdict", "note")]
  rownames(out) <- NULL

  nf <- sum(out$verdict == "FAIL"); nw <- sum(out$verdict == "WARN")
  overall <- if (nf) "FAIL" else if (nw) "WARN" else
             if (all(out$verdict == "UNTESTED")) "UNTESTED" else "OK"
  if (verbose) {
    if (overall == "OK")
      message("  the data are consistent with the graph on every claim tested")
    else if (overall == "WARN")
      message("  ", nw, " claim(s) reached significance but stayed under ",
              "min_effect = ", min_effect, "; the graph is not contradicted ",
              "in any way worth acting on")
    else if (overall == "FAIL")
      message("  ", nf, " claim(s) CONTRADICTED: the graph asserts an ",
              "independence the data do not show. Either an arrow is missing, ",
              "or a variable does not measure what its name supposes. ",
              "Effects estimated from this graph inherit that.")
  }
  structure(out, verdict = overall, alpha = alpha, min_effect = min_effect)
}

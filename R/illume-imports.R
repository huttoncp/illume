#' Internal imports and global-variable declarations
#'
#' `RTMB` is imported wholesale rather than selectively. The likelihood is
#' written with its automatic-differentiation types, and the arithmetic,
#' matrix and indexing operations inside a taped function dispatch on those
#' types. Importing only a few named functions would leave those operator
#' methods unavailable. Distribution functions that RTMB also provides
#' (`pnorm`, `qnorm`, `pchisq`, `qchisq`, `cov2cor`, `sdreport`) are left to
#' RTMB rather than taken from stats: its versions accept both ordinary
#' numbers and AD types, so masking them would break code inside the tape.
#'
#' @import RTMB
#' @importFrom TMB openmp
#' @importFrom splines splineDesign
#' @importFrom stats nlminb contr.sum quantile rnorm runif sd cor terms
#'   formula model.frame model.matrix delete.response reformulate na.omit
#'   setNames AIC BIC logLik nobs predict coef vcov ppoints lowess
#'   complete.cases aggregate .getXlevels printCoefmat getCall qqnorm
#'   qqline density
#' @importFrom graphics par plot points lines abline arrows legend axis title
#'   mtext polygon text plot.new segments strwidth
#' @importFrom grDevices adjustcolor hcl.colors
#' @importFrom utils head tail
#' @importFrom collapse fcount fmatch fmean fndistinct fnobs fquantile
#'   fsd fsum fvar group
#' @importFrom illumex ilm_glrm
#' @name illume-imports
#' @noRd
NULL

## Variables that RTMB::getAll() brings into the likelihood's scope, plus
## symbols used in non-standard evaluation.  R CMD check cannot see these,
## because they are never assigned in the ordinary way.
utils::globalVariables(c(
  "X", "Ycount", "wrow", "Tct", "grp", "Zl", "kind", "bas", "b_idx", "t_idx",
  "nlk", "dk", "wk", "npc", "ty", "rk", "dcor", "K", "C", "has_ar",
  "ar_idx", "idx1", "idx_t", "idx_lag", "n_g", "Tt", "qk", "brng", "lc_idx",
  "is_car", "is_rw", "ar_gap", "nre_ar", "has_dm", "disp_mu", "Zdisp", "gamma", "mu_pow",
  "has_rp", "Drp",
  ## the boundary-avoiding penalty: its weight, and which terms it applies to
  "pen_re", "pen_k",
  ## the zero part and the ordered-response block, which live in the TMB
  ## data and parameter lists rather than in this frame
  "Zzi", "gzi", "has_zi", "zi_hurdle", "zi_cont", "i_zero", "i_pos",
  "yzero",
  "has_ord", "zeta_raw", "ord_L", "ord_iu", "ord_il", "ord_mu", "ord_ml",
  "obs_lev", "obs_row",
  "beta", "theta", "bvec", "lchol", "lchol_ar", "rho_raw", "B_ar", "B",
  "logdisp", "n_disp", "yobs",
  "rl_re"
))

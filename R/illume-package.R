#' illume: Frequentist Regression and Mixed Models with Built-In Diagnostics
#'
#' A frequentist regression workflow aimed at **inference** rather than
#' prediction. One engine fits gaussian, binomial, Poisson, negative binomial
#' and multinomial models, so every model returns the same object and shares the
#' same methods and diagnostics.
#'
#' @section What problem this solves:
#' Fitting a range of models in R usually means a range of packages, each with
#' its own object, its own `anova()` semantics, its own `predict()` arguments and
#' its own idea of what a diagnostic plot is. `illume` fits them with one engine,
#' so the interface does not change as the model gets more complicated -- from
#' `y ~ x` through to a multinomial mixed model with smooths and correlated
#' errors.
#'
#' Where an exact small-sample reference exists it is used: a gaussian model
#' with no random effects reports t and F tests agreeing with [stats::lm()],
#' rather than the large-sample approximations the Laplace machinery would
#' otherwise give.
#'
#' The **multinomial** family closes a specific gap. [nnet::multinom()] handles
#' unordered categorical outcomes but has no random effects; [lme4::glmer()] and
#' `glmmTMB` have rich random-effect machinery but no multinomial family; `brms`
#' fits the full model but is far slower and answers a different inferential
#' question.
#'
#' @section Families:
#' `"gaussian"`, `"binomial"`, `"poisson"`, `"nbinom"` and `"multinomial"`; see
#' [lum_family()]. Zero-inflation, hurdle models and Tweedie are deliberately out
#' of scope, because `glmmTMB` covers those well.
#'
#' @section How the model is written:
#' With `J` outcome categories the model has `C = J - 1` free linear predictors.
#' `illume` uses **sum-to-zero** coding across categories, so each coefficient is
#' that category's deviation from the average across categories, *not* a contrast
#' against a baseline category. This differs from [nnet::multinom()] and `brms`,
#' which use baseline-category coding, and matters when reading output: see
#' [lum_coef_table()].
#'
#' @section Why the diagnostics are prominent:
#' This model class fails *quietly*. A fit can return entirely plausible
#' coefficients while its covariance matrix is unusable, so the standard errors
#' are meaningless even though nothing looks wrong. Every fit therefore carries a
#' table of checks (see `fit$checks`), each with a status, a plain-language
#' reason, and a suggested remedy. [summary.lum_model()] prints them.
#'
#' @section Main entry points:
#' \describe{
#'   \item{[lum_model()]}{fit a model from a formula, for any family}
#'   \item{[lum_family()]}{what each response distribution assumes}
#'   \item{[summary.lum_model()]}{coefficient table plus assumption checks}
#'   \item{[lum_anova()]}{Type II/III analysis of deviance for fixed effects}
#'   \item{[predict.lum_model()]}{category probabilities, conditional or population-averaged}
#'   \item{[lum_appraise()]}{residual diagnostic panels}
#'   \item{[lum_pb_lrt()]}{parametric-bootstrap likelihood-ratio test}
#' }
#'
#' @references
#' Agresti, A. (2013). *Categorical Data Analysis*, 3rd ed. Wiley.
#' (Chapter 8 covers baseline-category logit models.)
#'
#' Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. M. (2016).
#' TMB: Automatic differentiation and Laplace approximation.
#' *Journal of Statistical Software*, 70(5), 1--21. \doi{10.18637/jss.v070.i05}
#'
#' Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting linear
#' mixed-effects models using lme4. *Journal of Statistical Software*, 67(1),
#' 1--48. \doi{10.18637/jss.v067.i01}
#'
#' Brooks, M. E., Kristensen, K., van Benthem, K. J., Magnusson, A., Berg, C. W.,
#' Nielsen, A., Skaug, H. J., Maechler, M., & Bolker, B. M. (2017). glmmTMB
#' balances speed and flexibility among packages for zero-inflated generalized
#' linear mixed modeling. *The R Journal*, 9(2), 378--400.
#'
#' Wood, S. N. (2017). *Generalized Additive Models: An Introduction with R*,
#' 2nd ed. Chapman & Hall/CRC.
#'
#' @keywords internal
"_PACKAGE"

## Generate a grouped pkgdown reference index from the actual Rd topics, so
## that no topic can be left out of the index by hand.
## Run from the package root.

rd <- list.files("man", pattern = "[.]Rd$")
topic <- sub("[.]Rd$", "", rd)
## pkgdown indexes by the Rd \name, which is what the filename encodes except
## where roxygen escaped a character; read the name to be certain
nm <- vapply(file.path("man", rd), function(f) {
  x <- readLines(f, warn = FALSE)
  i <- grep("^\\\\name\\{", x)[1]
  if (is.na(i)) NA_character_ else sub("^\\\\name\\{(.*)\\}.*$", "\\1", x[i])
}, character(1))
topic <- ifelse(is.na(nm), topic, nm)
topic <- setdiff(topic, "illume-package")

## ---- groups, in the order they should appear -------------------------------
## Each entry: title, a one-line description, and a predicate over topic names.
grp <- list(
  list("Fitting models",
       "One engine, every family. The formula is lme4's and the object that comes back is the same whatever was fitted.",
       function(x) x %in% c("ilm_model", "ilm_fit", "ilm_family",
         "ilm_model_formula", "ilm_surv", "ilm_censor", "ilm_ar1", "ilm_car1",
         "ilm_fourier", "ilm_cyclic", "ilm_squeeze", "ilm_thresholds")),

  list("Other designs",
       "Identification strategies and sampling designs that are not a single regression.",
       function(x) grepl("^ilm_(iv|did|rdd|design|svy|mediate)", x)),

  list("Diagnostics",
       "Each check reports whether an assumption is consistent with the data, and names a remedy that exists in this package when it is not.",
       function(x) grepl("^ilm_(check_|appraise|rqr|binned|calibration|consistency|variogram|re_mahalanobis|scores|gauss_check)", x)),

  list("ANOVA",
       "Factorial and repeated-measures designs, specified by naming columns rather than by writing a formula with an error term.",
       function(x) grepl("^ilm_aov", x)),

  list("Inference and interpretation",
       "What the model says, on a scale someone can read.",
       function(x) grepl("^ilm_(anova|effects|emmeans|contrast|trends|ame|robust|vcov_cluster|denom_df|pb_lrt|rp_lrt|boot_|coef_table|se_fixef|zi_|scenario|interpret|translate)", x)),

  list("Design and power",
       "Before the data exist.",
       function(x) grepl("^ilm_power|^plot[.]ilm_power$", x)),

  list("Describing data",
       "Descriptive statistics, counts, and the things that are wrong with a data frame before any model sees it.",
       function(x) grepl("^ilm_(describe|counts|dupes|copies|wash_df|recode_errors|frame_issues)", x) &&
                   !grepl("_na", x)),

  list("Outliers and anomalies",
       "A value extreme for its own column, against a row implausible as a combination.",
       function(x) grepl("^ilm_(outliers|anomaly)", x)),

  list("Structure: reduce, cluster, profile",
       "The few directions a set of correlated columns shares, the groups in that space, and what distinguishes them.",
       function(x) grepl("^ilm_(reduce|cluster|profile|glrm)", x)),

  list("Missing data",
       "Diagnosing it, filling it in honestly, and pooling across the imputations.",
       function(x) grepl("^ilm_(impute|mi_pool)", x) || grepl("_na$|_na_all$", x)),

  list("Causal models",
       "A DAG, what it implies, and what it licenses you to say.",
       function(x) grepl("^ilm_(dag|adjust_sets|dsep)", x)),

  list("Plots",
       "Built on tinyplot, named for what they show.",
       function(x) grepl("^ilm_(plot|pick_geom|geom_spec)", x)),

  list("Simulation",
       "Data with a known structure, and draws from a fitted model.",
       function(x) grepl("^ilm_(sim|simulate|fitted|survival)$", x)),

  list("Working with other packages",
       "Registration shims and the methods that let the wider ecosystem dispatch on an illume fit.",
       function(x) grepl("^ilm_register|[.]ilm_model$", x))
)

assigned <- character(0)
lines <- c("reference:")
for (g in grp) {
  hit <- Filter(function(x) isTRUE(tryCatch(g[[3]](x), error = function(e) FALSE)),
                setdiff(topic, assigned))
  hit <- sort(hit)
  if (!length(hit)) next
  assigned <- c(assigned, hit)
  lines <- c(lines,
             paste0("- title: \"", g[[1]], "\""),   # quoted: titles may hold a colon
             paste0("  desc: >"),
             paste0("    ", g[[2]]),
             "  contents:",
             paste0("  - ", hit))
}

left <- sort(setdiff(topic, assigned))
if (length(left)) {
  lines <- c(lines, "- title: \"Shared parameters and print methods\"",
             "  desc: >",
             "    Documentation shared across functions, and methods you call by printing rather than by name.",
             "  contents:",
             paste0("  - ", left))
}

head <- c(
  "url: https://huttoncp.github.io/illume/",
  "template:",
  "  bootstrap: 5",
  "  bslib:",
  "    primary: \"#2a6f97\"",
  "",
  "home:",
  "  title: Exploration and frequentist inference in one toolkit",
  "",
  "navbar:",
  "  structure:",
  "    left:  [intro, reference, articles, news]",
  "    right: [search, github]",
  "",
  "articles:",
  "- title: Start here",
  "  navbar: ~",
  "  contents:",
  "  - illume",
  "  - workflow",
  "- title: Modelling",
  "  navbar: Modelling",
  "  contents:",
  "  - regression-models",
  "  - causal-models",
  "  - effect-size-and-power",
  "  - anova",
  "- title: Exploration",
  "  navbar: Exploration",
  "  contents:",
  "  - exploring-data",
  "  - profiling",
  "  - anomaly-detection",
  "  - missing-data",
  "")

writeLines(c(head, lines), "_pkgdown.yml")
cat("topics:", length(topic), " assigned:", length(assigned),
    " unassigned:", length(left), "\n")
if (length(left)) cat("UNASSIGNED:", paste(left, collapse = ", "), "\n")

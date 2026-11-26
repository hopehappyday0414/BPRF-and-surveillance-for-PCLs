###############################################################################
# BPRF_binary_adjusted_two_stage_IHME_dichotomous_final.R
# Adjusted-only IHME dichotomous pipeline: each stage runs CovFinder and
# MR-BRT on the full eligible dataset with the same stage-specific inlier_pct;
# BPRF/ROS are calculated directly from that fitted stage model.
# no conventional meta-analysis/I-squared;
# manuscript Egger=S4 all-observation (pre-exclusion) result;
# w_soln >= 0.1 defines the post-trimming diagnostic set.
#
# Strategy B: structural/estimability exclusions only, without the stricter
# floor(k/10) or k-5 covariate-capacity limits, plus gamma-boundary tolerance.
# Complete PCN binary BPRF analysis for two adjusted model scenarios.
#
# TWO PRESPECIFIED ADJUSTED MODEL STAGES
# --------------------------------
#   Stage 3 / S3: intercept random effect + selected cv fixed effects,
#     trim = 1.0
#     Purpose: CovFinder and MR-BRT both use inlier_pct=1.0.
#      Candidate cv set includes cv_adj and excludes cov_adj*, cv_adj_L0,
#      and cv_outcome.
#
#   Stage 4 / S4: intercept random effect + selected cv fixed effects,
#     trim = 0.9
#     Purpose: CovFinder and MR-BRT both use inlier_pct=0.9. This is the
#     prespecified conservative primary scenario and is scored directly.
#
# Notes:
#   - cv_mri, cv_ct, cv_eus, and cv_mixed are treated as mutually exclusive.
#     If CovFinder selects more than one, only the one with the strongest
#     univariate absolute association is retained.
#   - The final adjusted effect is predicted at selected covariates = 0, so it
#     represents the expected effect for an ideal reference study after bias
#     covariate adjustment.
#   - Final rating uses ROS = ln(BPRF) / 2 for dichotomous risks.
#   - cv_disease_IPMN and cv_disease_MCN are candidate covariates for all
#     indicators.
#   - cv_reader_quality is eligible only for imaging indicators. It is
#     excluded for clinical, biochemical, demographic, and composite
#     indicators because reader-based image-assessment quality cannot explain
#     ascertainment of those non-imaging exposures.
#   - cv_subgroup_imaging and cv_subgroup_pathology are candidate covariates
#     only when the indicator name contains "mpd".
#   - Reported model uncertainty intervals use the 2.5th and 97.5th
#     percentiles (95% UI).
#   - BPRF/ROS scoring remains based on the official IHME
#     DichotomousScorelator 5th/95th-percentile tails:
#       gamma_ub = gamma + 2 * sqrt(gamma_var).
#   - If gamma is on the zero boundary and limetr asserts gamma > 0 for Fisher
#     information, the fallback evaluates Fisher at gamma = 1e-10 and records
#     this positive-side boundary approximation in QC.
#   - ROS = score_log / 2.
#   - Two MR-BRT interval families are kept separate. The published BPRF label
#     "95% UI with gamma" corresponds to (b), not (a):
#       (a) no_gamma_OR_*: beta-only 95% CI from each MR-BRT score model;
#       (b) MRBRT_with_gamma_UI95_* / outer_OR_*: conservative 95% UI using
#           gamma_ub = gamma + 2*sqrt(gamma_var).
#     The inner_OR_* interval includes only the fitted gamma point estimate
#     and is retained as a diagnostic; when gamma=0 it is expected to equal
#     the no-gamma interval and must not be reported as "with gamma".
###############################################################################
#
# SCORE/PLOT PROVENANCE LOCK
# --------------------------
# The BPRF score, ROS, star rating, and forest-plot annotations in this script
# are intentionally inherited from the user's three previous validated scripts.
#
# Canonical quantities:
#   score_log       = direction-reflected outer-tail score used for evidence
#   BPRF_OR         = exp(score_log)
#   ROS             = score_log / 2
#   ROS_rounded     = round(ROS, 4)
#   stars           = legacy thresholds based on ROS
#
# Plot-only directional quantities:
#   score_curve_log = outer 5th percentile when beta>0, outer 95th percentile
#                     when beta<0, WITHOUT reflection
#   BPRF_curve_OR   = exp(score_curve_log)
#
# Therefore:
#   - evidence grading uses score_log -> ROS -> stars
#   - the BPRF vertical line on the original OR axis uses score_curve_log
#   - the plot title annotates ROS_rounded and stars exactly as before
#   - no alternative/new evidence score is introduced
###############################################################################

rm(list = ls())

library(reticulate)
library(data.table)
library(dplyr)

mrtool <- import("mrtool")
evidence_score <- import("mrtool.evidence_score.scorelator")
other_sampling <- import("mrtool.core.other_sampling")

# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------

# BATCH MODE.
# All files matching *_all_extraction.csv in cleaned_data will be analyzed.
n_samples <- 1000L

project_dir <- "/home/hopehappy/PCN_BPRF"
in_dir <- file.path(project_dir, "cleaned_data")

save_dir <- file.path(project_dir, "results", "batch_v7_correct_with_gamma")
fig_dir <- file.path(project_dir, "figures", "batch_v7_correct_with_gamma")
dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Egger interpretation is formally enabled only when >=10 unique studies.
MIN_K_EGGER <- 10L
EGGER_ALPHA <- 0.05

# Continue to the next indicator if one indicator fails.
CONTINUE_ON_ERROR <- TRUE

mutually_exclusive_imaging_covs <- c("cv_mri", "cv_ct", "cv_eus", "cv_mixed")
core_bias_covs <- c(
  # Covariates potentially applicable to both imaging and non-imaging features.
  "cv_feature_in_selection",
  "cv_consecutive_sampling",
  "cv_studydesign",
  "cv_multicenter",
  "cv_adj",
  "cv_population"
  # IMPORTANT: surgical indication/criteria are intentionally excluded.
)
imaging_quality_covs <- "cv_reader_quality"
imaging_indicator_patterns <- c(
  "mpd", "mural", "nodule", "cyst_size", "cyst_growth", "growth_rate",
  "thick_wall", "wall", "sept", "duct", "location", "lymphaden",
  "solid_component", "enhanc", "calcification"
)
disease_scope_covs <- c(
  "cv_disease_IPMN",
  "cv_disease_MCN"
)
mpd_subgroup_covs <- c(
  "cv_subgroup_imaging",
  "cv_subgroup_pathology"
)

scenarios <- data.table(
  model_stage = c(3L, 4L),
  scenario_id = c(
    "S3_cvFE_trim1",
    "S4_cvFE_trim09"
  ),
  trm = c(1.0, 0.9),
  use_covfinder = c(TRUE, TRUE),
  covariate_mode = c("discover_by_stage", "discover_by_stage"),
  model_stage_label = c(
    "Stage 3: covariate-adjusted, no trimming",
    "Stage 4: covariate-adjusted, 10% trimming"
  ),
  description = c(
    "CovFinder + adjusted MR-BRT; both use trim = 1.0",
    "CovFinder + adjusted MR-BRT; both use trim = 0.9"
  )
)

# ---------------------------------------------------------------------------
# FINAL SCORE CALCULATION AND EVIDENCE RATING
# -------------------------------------------
# Calculation chain used throughout this script:
#   1. Estimate the reference-study log effect beta and its variance.
#   2. Estimate between-study variance gamma and Fisher variance gamma_var.
#   3. Protect against underestimated heterogeneity using:
#          gamma_ub = gamma + 2 * sqrt(gamma_var)
#   4. Construct central 95% UIs for reporting:
#          inner = qnorm(c(0.025, 0.975), beta,
#                        sqrt(beta_var + gamma))
#          outer = qnorm(c(0.025, 0.975), beta,
#                        sqrt(beta_var + gamma_ub))
#   5. Separately construct the IHME scoring tails:
#          bprf_bounds = qnorm(c(0.05, 0.95), beta,
#                              sqrt(beta_var + gamma_ub))
#      Select the tail closest to the null and reflect protective directions:
#          score_log = bprf_bounds[1] when beta > 0
#          score_log = -bprf_bounds[2] when beta < 0
#   6. Convert to the Burden of Proof risk function and risk-outcome score:
#          BPRF_OR = exp(score_log)
#          ROS = score_log / 2
#   7. Convert ROS to the prespecified 1-5 star evidence category.
#
# The displayed inner and outer intervals are central 95% UIs. They are not
# used as BPRF endpoints. BPRF scoring uses the separate 5th/95th-percentile
# bounds. BPRF_curve_OR is plot-only and remains on the original direction.
# ---------------------------------------------------------------------------

ros_from_bprf <- function(BPRF) {
  BPRF <- suppressWarnings(as.numeric(BPRF))
  if (!is.finite(BPRF) || BPRF <= 0) return(NA_real_)
  log(BPRF) / 2
}

classify_ros <- function(ROS) {
  if (!is.finite(ROS)) {
    return(list(
      stars = NA_integer_,
      star_label = NA_character_,
      interpretation = "Not estimable"
    ))
  }

  # One star: the conservative BPRF bound crosses the null (BPRF_OR < 1).
  if (ROS < 0) {
    return(list(
      stars = 1L,
      star_label = "\u2605",
      interpretation = "Insufficient or very weak evidence of association"
    ))
  }
  # Two stars: 0 <= ROS <= 0.14.
  if (ROS <= 0.14) {
    return(list(
      stars = 2L,
      star_label = "\u2605\u2605",
      interpretation = "Weak evidence of harmful association"
    ))
  }
  # Three stars: 0.14 < ROS <= 0.41.
  if (ROS <= 0.41) {
    return(list(
      stars = 3L,
      star_label = "\u2605\u2605\u2605",
      interpretation = "Moderate evidence of harmful association"
    ))
  }
  # Four stars: 0.41 < ROS <= 0.62.
  if (ROS <= 0.62) {
    return(list(
      stars = 4L,
      star_label = "\u2605\u2605\u2605\u2605",
      interpretation = "Strong evidence of harmful association"
    ))
  }
  # Five stars: ROS > 0.62.
  list(
    stars = 5L,
    star_label = "\u2605\u2605\u2605\u2605\u2605",
    interpretation = "Very strong evidence of harmful association"
  )
}

add_ros_rating <- function(dt, bprf_col = "BPRF_OR", score_log_col = NULL) {
  if (!is.null(score_log_col) && score_log_col %in% names(dt)) {
    dt[, ROS := {
      score_log <- suppressWarnings(as.numeric(get(score_log_col)))
      ifelse(is.finite(score_log), score_log / 2, NA_real_)
    }]
  } else {
    dt[, ROS := vapply(get(bprf_col), ros_from_bprf, numeric(1))]
  }
  dt[, ROS_rounded := round(ROS, 4)]

  rating <- lapply(dt$ROS, classify_ros)
  dt[, stars := vapply(rating, `[[`, integer(1), "stars")]
  dt[, star_label := vapply(rating, `[[`, character(1), "star_label")]
  dt[, interpretation := vapply(rating, `[[`, character(1), "interpretation")]
  dt
}

star_label <- function(stars) {
  if (is.na(stars)) return(NA_character_)
  paste0(stars, " star", ifelse(stars == 1L, "", "s"))
}

# Enforce exact consistency with the previous score/plot convention.
validate_legacy_score_fields <- function(score_dt, tolerance = 1e-10) {
  if (nrow(score_dt) < 1) stop("score_dt is empty")

  required <- c(
    "score_log", "BPRF_OR", "score_curve_log", "BPRF_curve_OR",
    "ROS", "ROS_rounded", "stars"
  )
  missing <- setdiff(required, names(score_dt))
  if (length(missing) > 0) {
    stop("Legacy score fields missing: ", paste(missing, collapse = ", "))
  }

  x <- score_dt[1]

  checks <- c(
    BPRF_OR_from_score_log =
      is.finite(x$score_log) &&
      is.finite(x$BPRF_OR) &&
      abs(x$BPRF_OR - exp(x$score_log)) <= tolerance,

    ROS_from_score_log =
      is.finite(x$score_log) &&
      is.finite(x$ROS) &&
      abs(x$ROS - x$score_log / 2) <= tolerance,

    ROS_rounded_from_ROS =
      is.finite(x$ROS) &&
      is.finite(x$ROS_rounded) &&
      abs(x$ROS_rounded - round(x$ROS, 4)) <= tolerance,

    BPRF_curve_from_score_curve_log =
      is.finite(x$score_curve_log) &&
      is.finite(x$BPRF_curve_OR) &&
      abs(x$BPRF_curve_OR - exp(x$score_curve_log)) <= tolerance
  )

  if (!all(checks)) {
    stop(
      "Legacy BPRF/ROS score consistency check failed: ",
      paste(names(checks)[!checks], collapse = ", ")
    )
  }

  invisible(TRUE)
}

safe_quantile <- function(x, probs) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

make_study_label <- function(dt) {
  if (all(c("author", "year") %in% names(dt))) {
    return(paste0(substr(as.character(dt$author), 1, 24), " (", dt$year, ")"))
  }
  as.character(dt$pmid)
}

save_dill <- function(object, filename) {
  dill <- import("dill")
  builtins <- import_builtins()
  f <- builtins$open(filename, "wb")
  dill$dump(object, f)
  f$close()
}

safe_gamma_sd <- function(model) {
  gamma <- model$gamma_soln
  out <- tryCatch({
    gamma_fisher <- as.matrix(model$lt$get_gamma_fisher(gamma))
    gamma_var <- solve(gamma_fisher)
    as.numeric(sqrt(diag(gamma_var))[1])
  }, error = function(e) NA_real_)
  if (!is.finite(out) || out < 0) NA_real_ else out
}

manual_beta_var <- function(model) {
  obs_se <- as.numeric(model$data$obs_se)
  gamma <- max(0, as.numeric(model$gamma_soln)[1])
  w <- as.numeric(model$w_soln)
  if (length(w) != length(obs_se)) w <- rep(1, length(obs_se))
  w[!is.finite(w) | w < 0] <- 0

  keep <- is.finite(obs_se) & obs_se > 0 & w > 0
  if (sum(keep) < 2) {
    keep <- is.finite(obs_se) & obs_se > 0
    w <- rep(1, length(obs_se))
  }

  1 / sum(w[keep] / (obs_se[keep]^2 + gamma))
}

ihme_formula_score <- function(model, indicator_name, trigger_error = NA_character_,
                               gamma_boundary_epsilon = 1e-10,
                               gamma_boundary_tolerance = 1.01e-8) {
  cov_index0 <- as.integer(model$get_cov_model_index("intercept"))
  x_ids <- as.integer(model$x_vars_indices[[cov_index0 + 1L]]) + 1L
  z_ids <- as.integer(model$z_vars_indices[[cov_index0 + 1L]]) + 1L
  if (length(x_ids) < 1 || length(z_ids) < 1) {
    stop("intercept beta/gamma index is unavailable")
  }

  beta <- as.numeric(model$beta_soln)[x_ids[1]]
  gamma_raw <- as.numeric(model$gamma_soln)[z_ids[1]]
  gamma_tol <- gamma_boundary_tolerance
  if (!is.finite(gamma_raw) || gamma_raw < -gamma_tol) {
    stop("materially negative or nonfinite gamma; refit is required")
  }
  gamma <- if (gamma_raw < 0) 0 else gamma_raw

  specs <- other_sampling$extract_simple_lme_specs(model)
  beta_hessian <- as.matrix(other_sampling$extract_simple_lme_hessian(specs))
  beta_cov <- solve(beta_hessian)
  beta_var <- as.numeric(beta_cov[x_ids[1], x_ids[1]])

  # Some limetr versions assert gamma > 0 inside get_gamma_fisher().  When the
  # fitted variance is exactly on the non-negative boundary, evaluate Fisher
  # information at a documented positive-side epsilon.  This estimates the
  # boundary limit; it does NOT set gamma_sd to zero.
  gamma_for_fisher <- as.numeric(model$gamma_soln)
  fisher_boundary_used <- gamma_for_fisher[z_ids[1]] <= gamma_boundary_epsilon
  if (fisher_boundary_used) {
    gamma_for_fisher[z_ids[1]] <- gamma_boundary_epsilon
  }
  gamma_fisher <- as.matrix(model$lt$get_gamma_fisher(gamma_for_fisher))
  gamma_cov <- solve(gamma_fisher)
  gamma_var <- as.numeric(gamma_cov[z_ids[1], z_ids[1]])
  if (!is.finite(beta_var) || beta_var < 0 ||
      !is.finite(gamma_var) || gamma_var < 0) {
    stop("beta/gamma uncertainty matrix is not estimable")
  }

  gamma_sd <- sqrt(gamma_var)
  gamma_ub <- gamma + 2 * gamma_sd
  # Three distinct central 95% intervals are retained:
  #   no_gamma: uncertainty in the pooled mean beta only;
  #   inner: beta uncertainty plus the fitted gamma;
  #   outer: beta uncertainty plus the Fisher upper quantile of gamma.
  no_gamma <- qnorm(c(0.025, 0.975), beta, sqrt(beta_var))
  inner <- qnorm(c(0.025, 0.975), beta, sqrt(beta_var + gamma))
  outer <- qnorm(c(0.025, 0.975), beta, sqrt(beta_var + gamma_ub))
  # BPRF is a conservative tail estimate, not a 95% interval endpoint.
  bprf_bounds <- qnorm(c(0.05, 0.95), beta, sqrt(beta_var + gamma_ub))
  score <- if (beta > 0) bprf_bounds[1] else -bprf_bounds[2]

  out <- data.table(
    score_method = "IHME_formula_fallback_after_Scorelator_error",
    score_log = score,
    BPRF_OR = exp(score),
    score_curve_log = if (beta > 0) bprf_bounds[1] else bprf_bounds[2],
    BPRF_curve_OR = exp(if (beta > 0) bprf_bounds[1] else bprf_bounds[2]),
    score_quantile = if (beta > 0) 0.05 else 0.95,
    beta = beta, beta_var = beta_var, beta_sd = sqrt(beta_var),
    gamma = gamma, gamma_var = gamma_var, gamma_sd = gamma_sd,
    gamma_ub = gamma_ub,
    no_gamma_lo = no_gamma[1], no_gamma_hi = no_gamma[2],
    no_gamma_OR_lo = exp(no_gamma[1]), no_gamma_OR_hi = exp(no_gamma[2]),
    inner_lo = inner[1], inner_hi = inner[2],
    outer_lo = outer[1], outer_hi = outer[2],
    inner_OR_lo = exp(inner[1]), inner_OR_hi = exp(inner[2]),
    outer_OR_lo = exp(outer[1]), outer_OR_hi = exp(outer[2]),
    bprf_q05 = bprf_bounds[1], bprf_q95 = bprf_bounds[2],
    bprf_OR_q05 = exp(bprf_bounds[1]), bprf_OR_q95 = exp(bprf_bounds[2]),
    score_direction = if (beta > 0) "harmful" else "protective_reflected",
    gamma_status = if (gamma_raw < 0) {
      "tiny_negative_clamped_to_zero_fisher_positive_limit"
    } else if (gamma == 0) {
      "zero_boundary_fisher_positive_limit"
    } else {
      "positive_with_estimable_uncertainty"
    },
    gamma_raw = gamma_raw,
    gamma_fisher_boundary_used = fisher_boundary_used,
    gamma_fisher_epsilon = ifelse(
      fisher_boundary_used, gamma_boundary_epsilon, NA_real_
    ),
    scorelator_trigger_error = trigger_error,
    formula_fallback_error = NA_character_,
    score_status = "ok_formula_fallback",
    score_error = trigger_error
  )
  add_ros_rating(out, score_log_col = "score_log")
}


# ---------------------------------------------------------------------------
# Direct model-matrix + Fisher fallback for complex fixed-effect MR-BRT models
# ---------------------------------------------------------------------------
#
# This fallback does NOT change the fitted MR-BRT model and does NOT remove any
# selected covariates. It only replaces uncertainty extraction when
# DichotomousScorelator and extract_simple_lme_* assert on a complex model.
#
# beta variance:
#   Var(beta) = (X' W X)^(-1), where
#   W_i = robust_weight_i / (SE_i^2 + gamma)
#
# For score-only refits after trimming, robust weights are 1 by construction.
#
# gamma variance:
#   extracted from the fitted limetr object's gamma Fisher information.
#   If gamma is on the boundary, Fisher is evaluated from the positive side.
# ---------------------------------------------------------------------------

ihme_direct_matrix_fisher_score <- function(
    model,
    indicator_name,
    selected_covs = character(0),
    trigger_error = NA_character_,
    gamma_boundary_grid = c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6),
    tolerance = 1e-10
) {

  beta_soln <- as.numeric(model$beta_soln)
  gamma_soln <- as.numeric(model$gamma_soln)

  if (length(beta_soln) < 1 || length(gamma_soln) < 1) {
    stop("beta/gamma solution vectors are unavailable")
  }

  beta <- beta_soln[1]
  gamma_raw <- gamma_soln[1]

  if (!is.finite(beta)) stop("intercept beta is non-finite")
  if (!is.finite(gamma_raw) || gamma_raw < -1e-7) {
    stop("gamma is materially negative or non-finite")
  }

  gamma <- max(0, gamma_raw)

  obs_se <- as.numeric(model$data$obs_se)
  n <- length(obs_se)

  if (n < 2 || any(!is.finite(obs_se)) || any(obs_se <= 0)) {
    stop("invalid observation standard errors")
  }

  # Reconstruct the fixed-effect design matrix in the same order as
  # fit_mrbrt_model(): intercept first, then selected fixed-effect covariates.
  model_df <- as.data.table(model$data$to_df())

  if (length(selected_covs) == 0) {
    X <- matrix(1, nrow = n, ncol = 1)
    colnames(X) <- "intercept"
  } else {
    missing_covs <- setdiff(selected_covs, names(model_df))
    if (length(missing_covs) > 0) {
      stop(
        "selected covariates missing from model data: ",
        paste(missing_covs, collapse = ", ")
      )
    }

    X <- cbind(
      intercept = 1,
      as.matrix(model_df[, ..selected_covs])
    )
  }

  if (nrow(X) != n) stop("design matrix and observation count are misaligned")
  if (qr(X)$rank < ncol(X)) stop("fixed-effect design matrix is rank deficient")

  robust_w <- tryCatch(
    as.numeric(model$w_soln),
    error = function(e) rep(1, n)
  )

  if (length(robust_w) != n) robust_w <- rep(1, n)
  robust_w[!is.finite(robust_w) | robust_w < 0] <- 0

  # Use the robust weights from the stage-specific MR-BRT model. This supports
  # direct Fisher fallback for both trimmed and untrimmed adjusted scenarios.
  keep <- is.finite(obs_se) & obs_se > 0 & robust_w > 0
  if (sum(keep) <= ncol(X)) {
    stop("insufficient positive-weight observations for beta covariance")
  }

  v <- obs_se[keep]^2 + gamma
  W <- robust_w[keep] / v
  Xk <- X[keep, , drop = FALSE]

  XtWX <- crossprod(Xk, Xk * W)

  beta_cov <- tryCatch(
    solve(XtWX),
    error = function(e) {
      # QR-based generalized inverse is deliberately not used here because a
      # singular matrix means the adjusted intercept variance is not identifiable.
      stop("beta information matrix is singular: ", conditionMessage(e))
    }
  )

  beta_var <- as.numeric(beta_cov[1, 1])
  if (!is.finite(beta_var) || beta_var < 0) {
    stop("intercept beta variance is not estimable")
  }

  # Extract gamma Fisher information. For boundary gamma, test a short grid of
  # positive-side epsilon values because limetr versions differ in their
  # strict positivity tolerance.
  gamma_fisher <- NULL
  fisher_epsilon_used <- NA_real_
  fisher_errors <- character(0)

  gamma_candidates <- if (gamma_raw <= max(gamma_boundary_grid)) {
    gamma_boundary_grid
  } else {
    gamma_raw
  }

  for (g_eval in gamma_candidates) {
    gamma_eval <- gamma_soln
    gamma_eval[1] <- max(g_eval, .Machine$double.eps)

    gf <- tryCatch(
      as.matrix(model$lt$get_gamma_fisher(gamma_eval)),
      error = function(e) {
        fisher_errors <<- c(
          fisher_errors,
          paste0("gamma=", format(gamma_eval[1]), ": ", conditionMessage(e))
        )
        NULL
      }
    )

    if (!is.null(gf) &&
        all(dim(gf) >= c(1, 1)) &&
        is.finite(gf[1, 1]) &&
        gf[1, 1] > 0) {
      gamma_fisher <- gf
      fisher_epsilon_used <- if (gamma_raw <= max(gamma_boundary_grid)) {
        gamma_eval[1]
      } else {
        NA_real_
      }
      break
    }
  }

  if (is.null(gamma_fisher)) {
    stop(
      "gamma Fisher information is not estimable",
      if (length(fisher_errors) > 0) {
        paste0(" | attempts: ", paste(fisher_errors, collapse = " || "))
      } else {
        ""
      }
    )
  }

  gamma_cov <- tryCatch(
    solve(gamma_fisher),
    error = function(e) stop(
      "gamma Fisher matrix is singular: ", conditionMessage(e)
    )
  )

  gamma_var <- as.numeric(gamma_cov[1, 1])

  if (!is.finite(gamma_var) || gamma_var < 0) {
    stop("gamma variance from Fisher information is not estimable")
  }

  gamma_sd <- sqrt(gamma_var)
  gamma_ub <- gamma + 2 * gamma_sd

  no_gamma <- qnorm(
    c(0.025, 0.975),
    mean = beta,
    sd = sqrt(beta_var)
  )

  inner <- qnorm(
    c(0.025, 0.975),
    mean = beta,
    sd = sqrt(beta_var + gamma)
  )

  outer <- qnorm(
    c(0.025, 0.975),
    mean = beta,
    sd = sqrt(beta_var + gamma_ub)
  )

  bprf_bounds <- qnorm(
    c(0.05, 0.95),
    mean = beta,
    sd = sqrt(beta_var + gamma_ub)
  )

  score <- if (beta > 0) bprf_bounds[1] else -bprf_bounds[2]
  score_curve_log <- if (beta > 0) bprf_bounds[1] else bprf_bounds[2]

  out <- data.table(
    score_method = "IHME_direct_design_matrix_plus_gamma_Fisher_fallback",
    score_log = score,
    BPRF_OR = exp(score),
    score_curve_log = score_curve_log,
    BPRF_curve_OR = exp(score_curve_log),
    score_quantile = if (beta > 0) 0.05 else 0.95,

    beta = beta,
    beta_var = beta_var,
    beta_sd = sqrt(beta_var),

    gamma = gamma,
    gamma_var = gamma_var,
    gamma_sd = gamma_sd,
    gamma_ub = gamma_ub,

    no_gamma_lo = no_gamma[1],
    no_gamma_hi = no_gamma[2],
    no_gamma_OR_lo = exp(no_gamma[1]),
    no_gamma_OR_hi = exp(no_gamma[2]),

    inner_lo = inner[1],
    inner_hi = inner[2],
    outer_lo = outer[1],
    outer_hi = outer[2],

    inner_OR_lo = exp(inner[1]),
    inner_OR_hi = exp(inner[2]),
    outer_OR_lo = exp(outer[1]),
    outer_OR_hi = exp(outer[2]),
    bprf_q05 = bprf_bounds[1],
    bprf_q95 = bprf_bounds[2],
    bprf_OR_q05 = exp(bprf_bounds[1]),
    bprf_OR_q95 = exp(bprf_bounds[2]),

    score_direction = if (beta > 0) "harmful" else "protective_reflected",

    gamma_status = if (gamma_raw <= 0) {
      "boundary_positive_side_direct_fisher"
    } else {
      "positive_direct_fisher"
    },

    gamma_raw = gamma_raw,
    gamma_fisher_boundary_used = is.finite(fisher_epsilon_used),
    gamma_fisher_epsilon = fisher_epsilon_used,

    scorelator_trigger_error = trigger_error,
    formula_fallback_error = trigger_error,
    score_status = "ok_direct_matrix_fisher_fallback",
    score_error = trigger_error
  )

  add_ros_rating(out, score_log_col = "score_log")
}


ihme_bprf_score <- function(model, indicator_name, selected_covs = character(0), tolerance = 1e-7) {
  tryCatch({
    # Primary path: obtain beta, beta variance, gamma and gamma variance from
    # the official IHME DichotomousScorelator object.
    sc <- evidence_score$DichotomousScorelator(
      model,
      cov_name = "intercept",
      name = indicator_name
    )
    beta <- as.numeric(sc$beta)
    beta_var <- as.numeric(sc$beta_var)
    gamma <- max(0, as.numeric(sc$gamma))
    gamma_var <- as.numeric(sc$gamma_var)

    if (!is.finite(beta) || !is.finite(beta_var) || beta_var < 0) {
      stop("beta or beta variance is not estimable")
    }
    if (!is.finite(gamma_var) || gamma_var < 0) {
      stop("gamma uncertainty is not estimable; it must not be replaced by zero")
    }

    # Heterogeneity uncertainty: use the upper approximation gamma + 2 SD.
    gamma_sd <- sqrt(gamma_var)
    gamma_ub <- gamma + 2 * gamma_sd
    inner_sd <- sqrt(beta_var + gamma)
    outer_sd <- sqrt(beta_var + gamma_ub)
    # Central 95% intervals from the SAME score model. The no-gamma interval
    # represents uncertainty in beta alone; inner includes fitted gamma; outer
    # uses the Fisher upper quantile of gamma and is the conservative interval.
    no_gamma <- qnorm(c(0.025, 0.975), mean = beta, sd = sqrt(beta_var))
    inner <- qnorm(c(0.025, 0.975), mean = beta, sd = inner_sd)
    outer <- qnorm(c(0.025, 0.975), mean = beta, sd = outer_sd)
    # Separate 5th/95th-percentile bounds used only for BPRF/ROS scoring.
    bprf_bounds <- qnorm(c(0.05, 0.95), mean = beta, sd = outer_sd)

    direction <- if (beta > 0) "harmful" else "protective_reflected"
    # Select the outer tail closest to the null. For beta < 0, reflect the
    # protective direction so that evidence strength uses a common scale.
    score_manual <- if (beta > 0) bprf_bounds[1] else -bprf_bounds[2]
    # Official Scorelator score. The following check requires agreement with
    # the independently reconstructed IHME formula above.
    score <- as.numeric(sc$get_score(use_gamma_ub = TRUE))
    if (!is.finite(score) || abs(score - score_manual) > tolerance) {
      stop(sprintf(
        "Scorelator/manual IHME check failed: tool=%s, manual=%s",
        format(score), format(score_manual)
      ))
    }

    out <- data.table(
      score_method = "IHME_DichotomousScorelator_gamma_plus_2SD",
      score_log = score,
      BPRF_OR = exp(score),
      score_curve_log = if (beta > 0) bprf_bounds[1] else bprf_bounds[2],
      BPRF_curve_OR = exp(if (beta > 0) bprf_bounds[1] else bprf_bounds[2]),
      score_quantile = if (beta > 0) 0.05 else 0.95,
      beta = beta,
      beta_var = beta_var,
      beta_sd = sqrt(beta_var),
      gamma = gamma,
      gamma_var = gamma_var,
      gamma_sd = gamma_sd,
      gamma_ub = gamma_ub,
      no_gamma_lo = no_gamma[1],
      no_gamma_hi = no_gamma[2],
      no_gamma_OR_lo = exp(no_gamma[1]),
      no_gamma_OR_hi = exp(no_gamma[2]),
      inner_lo = inner[1],
      inner_hi = inner[2],
      outer_lo = outer[1],
      outer_hi = outer[2],
      inner_OR_lo = exp(inner[1]),
      inner_OR_hi = exp(inner[2]),
      outer_OR_lo = exp(outer[1]),
      outer_OR_hi = exp(outer[2]),
      bprf_q05 = bprf_bounds[1],
      bprf_q95 = bprf_bounds[2],
      bprf_OR_q05 = exp(bprf_bounds[1]),
      bprf_OR_q95 = exp(bprf_bounds[2]),
      score_direction = direction,
      gamma_status = ifelse(
        gamma == 0,
        "zero_boundary_with_estimable_uncertainty",
        "positive_with_estimable_uncertainty"
      ),
      gamma_raw = gamma,
      gamma_fisher_boundary_used = FALSE,
      gamma_fisher_epsilon = NA_real_,
      scorelator_trigger_error = NA_character_,
      formula_fallback_error = NA_character_,
      score_status = "ok",
      score_error = NA_character_
    )
    # Final chain: score_log -> BPRF_OR -> ROS -> star category.
    add_ros_rating(out)
  }, error = function(e) {
    fallback_error <- NA_character_
    fallback <- tryCatch(
      ihme_formula_score(model, indicator_name, conditionMessage(e)),
      error = function(e2) {
        fallback_error <<- conditionMessage(e2)
        NULL
      }
    )
    if (!is.null(fallback)) return(fallback)

    # Third-level fallback: keep the fitted model and selected covariates, but
    # extract beta/gamma uncertainty directly from the design matrix and the
    # limetr gamma Fisher information.
    direct_fallback_error <- NA_character_
    direct_fallback <- tryCatch(
      ihme_direct_matrix_fisher_score(
        model = model,
        indicator_name = indicator_name,
        selected_covs = selected_covs,
        trigger_error = paste0(
          "Scorelator: ", conditionMessage(e),
          " | simple-LME fallback: ", fallback_error
        )
      ),
      error = function(e_direct) {
        direct_fallback_error <<- conditionMessage(e_direct)
        NULL
      }
    )

    if (!is.null(direct_fallback)) return(direct_fallback)

    beta_raw <- tryCatch(
      as.numeric(model$beta_soln)[1], error = function(e3) NA_real_
    )
    gamma_raw <- tryCatch(
      as.numeric(model$gamma_soln)[1], error = function(e3) NA_real_
    )
    gamma_fisher_11 <- tryCatch({
      gamma_eval <- as.numeric(model$gamma_soln)
      gamma_eval[1] <- max(1e-10, gamma_eval[1])
      as.numeric(as.matrix(model$lt$get_gamma_fisher(gamma_eval))[1, 1])
    }, error = function(e3) NA_real_)

    out <- data.table(
      score_method = "IHME_DichotomousScorelator_gamma_plus_2SD",
      score_log = NA_real_,
      BPRF_OR = NA_real_,
      beta = beta_raw,
      beta_var = NA_real_,
      beta_sd = NA_real_,
      gamma = gamma_raw,
      gamma_var = NA_real_,
      gamma_sd = NA_real_,
      gamma_ub = NA_real_,
      inner_lo = NA_real_,
      inner_hi = NA_real_,
      outer_lo = NA_real_,
      outer_hi = NA_real_,
      inner_OR_lo = NA_real_,
      inner_OR_hi = NA_real_,
      outer_OR_lo = NA_real_,
      outer_OR_hi = NA_real_,
      bprf_q05 = NA_real_,
      bprf_q95 = NA_real_,
      bprf_OR_q05 = NA_real_,
      bprf_OR_q95 = NA_real_,
      score_direction = NA_character_,
      gamma_status = "not_estimable",
      gamma_raw = gamma_raw,
      gamma_fisher_11 = gamma_fisher_11,
      gamma_fisher_boundary_used = is.finite(gamma_raw) && gamma_raw <= 1.01e-8,
      gamma_fisher_epsilon = ifelse(
        is.finite(gamma_raw) && gamma_raw <= 1.01e-8, 1e-10, NA_real_
      ),
      scorelator_trigger_error = conditionMessage(e),
      formula_fallback_error = fallback_error,
      score_status = "failed",
      score_error = paste0(
        "Scorelator: ", conditionMessage(e),
        " | simple-LME fallback: ", fallback_error,
        " | direct matrix/Fisher fallback: ", direct_fallback_error
      )
    )
    add_ros_rating(out)
  })
}

exclude_unwanted_covariates <- function(covs) {
  covs <- covs[!grepl("^cov_adj", covs, ignore.case = TRUE)]
  covs <- covs[!tolower(covs) %in% c("cv_adj_l0", "cv_outcome", "cov_outcome")]

  # Surgical criteria/indication are part of the selection mechanism of a
  # surgical cohort and are therefore NOT treated as ordinary bias covariates.
  surgical_pattern <- paste(
    c(
      "surgical_criteria", "surgery_criteria", "surgical_indication",
      "surgery_indication", "resection_criteria", "resection_indication",
      "absolute_indication", "operation_criteria", "operative_indication"
    ),
    collapse = "|"
  )
  covs <- covs[!grepl(surgical_pattern, covs, ignore.case = TRUE)]
  unique(covs)
}

is_mpd_indicator <- function(indicator_name) {
  # MPD-related analyses include both duct-diameter indicators and IPMN
  # duct-subtype indicators.  The latter do not necessarily contain the
  # literal string "mpd" in their analysis names, but imaging- versus
  # pathology-defined subtype is still a prespecified source of bias.
  grepl("mpd", indicator_name, ignore.case = TRUE) ||
    tolower(indicator_name) %in% c(
      "md_ipmn",
      "mixed_ipmn",
      "md_or_mixed_ipmn"
    )
}

is_imaging_indicator <- function(indicator_name) {
  grepl(paste(imaging_indicator_patterns, collapse = "|"),
        indicator_name, ignore.case = TRUE)
}

candidate_covariates_for_indicator <- function(dt, indicator_name) {
  allowed <- c(core_bias_covs, disease_scope_covs)
  if (is_imaging_indicator(indicator_name)) {
    allowed <- c(
      allowed,
      imaging_quality_covs,
      mutually_exclusive_imaging_covs
    )
  }
  if (is_mpd_indicator(indicator_name)) {
    allowed <- c(allowed, mpd_subgroup_covs)
  }
  intersect(unique(allowed), names(dt))
}

screen_covariates <- function(dt, covs, min_group_n = 2L) {
  rbindlist(lapply(covs, function(cv) {
    x <- suppressWarnings(as.numeric(dt[[cv]]))
    valid <- x[is.finite(x)]
    values <- sort(unique(valid))
    n0 <- sum(valid == 0)
    n1 <- sum(valid == 1)
    reason <- if (length(valid) == 0) {
      "all_missing"
    } else if (length(valid) < nrow(dt)) {
      "contains_missing_or_nonfinite"
    } else if (!all(values %in% c(0, 1))) {
      "not_binary_0_1"
    } else if (length(values) < 2) {
      "single_value_or_zero_variance"
    } else if (n0 < min_group_n) {
      "reference_group_too_small"
    } else if (n1 < min_group_n) {
      "exposed_group_too_small"
    } else {
      "eligible"
    }
    data.table(
      covariate = cv, n_nonmissing = length(valid), n_0 = n0, n_1 = n1,
      eligible = reason == "eligible", exclusion_reason = reason
    )
  }), fill = TRUE)
}

design_is_full_rank <- function(dt, covs) {
  if (length(covs) == 0) return(TRUE)
  x <- cbind(intercept = 1, as.matrix(dt[, ..covs]))
  all(is.finite(x)) && qr(x)$rank == ncol(x)
}

cov_strength_score <- function(dt, cov) {
  x <- dt[[cov]]
  y <- dt$ln_effect
  se <- dt$ln_se

  ok <- is.finite(x) & is.finite(y) & is.finite(se) & se > 0
  x <- x[ok]
  y <- y[ok]
  w <- 1 / se[ok]^2

  if (length(unique(x[!is.na(x)])) < 2) return(0)
  abs(weighted.mean(y[x == 1], w[x == 1], na.rm = TRUE) -
        weighted.mean(y[x == 0], w[x == 0], na.rm = TRUE))
}

enforce_mutually_exclusive_covariates <- function(selected_covs, dt) {
  selected_covs <- unique(selected_covs)
  imaging_selected <- intersect(selected_covs, mutually_exclusive_imaging_covs)

  if (length(imaging_selected) <= 1) return(selected_covs)

  strengths <- sapply(imaging_selected, function(cov) cov_strength_score(dt, cov))
  keep <- names(which.max(strengths))
  drop <- setdiff(imaging_selected, keep)

  message("  Mutually exclusive imaging covariates selected: ",
          paste(imaging_selected, collapse = ", "),
          ". Keeping ", keep, " and dropping ", paste(drop, collapse = ", "), ".")

  unique(c(setdiff(selected_covs, drop), keep))
}

enforce_reference_dummy_pair <- function(selected_covs, dt, pair) {
  selected_pair <- intersect(pair, selected_covs)

  if (length(selected_pair) < 2) {
    return(selected_covs)
  }

  pair_sum <- rowSums(
    as.data.frame(dt[, ..pair]),
    na.rm = FALSE
  )

  if (all(pair_sum == 1)) {
    strengths <- sapply(
      pair,
      function(x) cov_strength_score(dt, x)
    )

    keep <- names(which.max(strengths))
    drop <- setdiff(pair, keep)

    message(
      "Complementary dummy variables detected: ",
      paste(pair, collapse = ", "),
      ". Keeping ", keep,
      " and dropping ", paste(drop, collapse = ", ")
    )

    selected_covs <- setdiff(selected_covs, drop)
  }

  selected_covs
}

select_covariates <- function(dt, trm) {
  candidate_covs <- candidate_covariates_for_indicator(dt, indicator)
  audit <- screen_covariates(dt, candidate_covs)
  eligible_covs_initial <- audit[eligible == TRUE, covariate]
  k_effective <- max(0L, floor(nrow(dt) * trm))

  if (length(eligible_covs_initial) == 0) {
    audit[, `:=`(
      selected = FALSE,
      maximum_allowed_covariates = max(0L, k_effective - 2L),
      k_effective = k_effective,
      covfinder_candidate_final = FALSE,
      covfinder_retry_dropped = FALSE,
      covfinder_retry_count = 0L,
      covfinder_last_error = NA_character_
    )]
    return(list(
      covariates = character(0),
      audit = audit,
      covfinder_status = "no_eligible_covariates",
      covfinder_retry_count = 0L,
      covfinder_dropped = character(0)
    ))
  }

  # -------------------------------------------------------------------------
  # CovFinder estimability safeguard
  #
  # CovFinder itself can fail with numpy.linalg.LinAlgError: Singular matrix
  # before the final MR-BRT fit is reached. We therefore first try the full
  # eligible candidate set. Only when CovFinder itself is numerically
  # non-estimable do we remove the weakest candidate group and retry.
  #
  # This is NOT p-value screening. Candidate strength is used only to choose
  # which redundant/weakest candidate is removed after a numerical failure.
  # Reference-coded disease/MPD dummy sets are treated as indivisible groups.
  # -------------------------------------------------------------------------

  make_candidate_groups <- function(covs) {
    covs <- unique(covs)
    groups <- list()
    used <- character(0)

    for (pair in list(disease_scope_covs, mpd_subgroup_covs)) {
      present <- intersect(pair, covs)
      if (length(present) > 0) {
        groups[[paste0("PAIR:", paste(pair, collapse = "+"))]] <- present
        used <- c(used, present)
      }
    }

    remaining <- setdiff(covs, used)
    for (cv in remaining) groups[[cv]] <- cv
    groups
  }

  drop_weakest_candidate_group <- function(covs) {
    groups <- make_candidate_groups(covs)
    if (length(groups) == 0) {
      return(list(covariates = character(0), dropped = character(0)))
    }

    strengths <- vapply(groups, function(g) {
      vals <- vapply(g, function(cv) cov_strength_score(dt, cv), numeric(1))
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) return(NA_real_)
      max(vals)
    }, numeric(1))

    if (all(!is.finite(strengths))) {
      victim <- tail(names(groups), 1)
    } else {
      tmp <- strengths
      tmp[!is.finite(tmp)] <- Inf
      victim <- names(which.min(tmp))[1]
    }

    dropped <- groups[[victim]]
    list(
      covariates = setdiff(covs, dropped),
      dropped = dropped
    )
  }

  run_covfinder_once <- function(covs) {
    data_covfinder <- mrtool$MRData()
    data_covfinder$load_df(
      data = dt,
      col_obs = "ln_effect",
      col_obs_se = "ln_se",
      col_covs = as.list(covs),
      col_study_id = "pmid"
    )

    covfinder <- mrtool$CovFinder(
      data = data_covfinder,
      covs = as.list(covs),
      pre_selected_covs = list("intercept"),
      normalized_covs = FALSE,
      num_samples = 1000L,
      power_range = list(-4, 4),
      power_step_size = 0.05,
      inlier_pct = trm,
      laplace_threshold = 1e-5,
      # CovFinder tests whether the central (1-alpha) interval excludes zero.
      # alpha=0.05 therefore uses the 2.5th and 97.5th percentiles.
      alpha = 0.05
    )
    covfinder$select_covs(verbose = FALSE)
    covfinder
  }

  eligible_covs <- unique(eligible_covs_initial)
  covfinder_dropped <- character(0)
  covfinder_retry_count <- 0L
  covfinder_last_error <- NA_character_
  covfinder <- NULL

  repeat {
    # With zero candidate covariates there is nothing left for CovFinder to
    # select; return an intercept-only estimability status rather than crash.
    if (length(eligible_covs) == 0) break

    cf_try <- tryCatch(
      run_covfinder_once(eligible_covs),
      error = function(e) e
    )

    if (!inherits(cf_try, "error")) {
      covfinder <- cf_try
      break
    }

    covfinder_last_error <- conditionMessage(cf_try)

    singular_like <- grepl(
      "singular|LinAlgError|not positive definite|ill-conditioned|rank",
      covfinder_last_error,
      ignore.case = TRUE
    )

    if (!singular_like) {
      stop(cf_try)
    }

    prune <- drop_weakest_candidate_group(eligible_covs)
    if (length(prune$dropped) == 0 ||
        identical(prune$covariates, eligible_covs)) {
      break
    }

    covfinder_retry_count <- covfinder_retry_count + 1L
    covfinder_dropped <- unique(c(covfinder_dropped, prune$dropped))

    message(
      "  CovFinder singular candidate set; dropping for estimability: ",
      paste(prune$dropped, collapse = ", "),
      " | previous error: ", covfinder_last_error
    )

    eligible_covs <- prune$covariates
  }

  if (is.null(covfinder)) {
    selected_covs <- character(0)
    covfinder_status <- if (length(eligible_covs_initial) > 0) {
      "covfinder_not_estimable_after_candidate_pruning"
    } else {
      "no_eligible_covariates"
    }
  } else {
    selected_covs <- as.character(covfinder$selected_covs)
    selected_covs <- selected_covs[!grepl("intercept", selected_covs)]
    selected_covs <- exclude_unwanted_covariates(selected_covs)
    selected_covs <- enforce_mutually_exclusive_covariates(
      selected_covs,
      dt
    )
    covfinder_status <- if (covfinder_retry_count > 0) {
      "ok_after_candidate_pruning"
    } else {
      "ok"
    }
  }

  # Keep reference-coded three-level variables as indivisible pairs, but only
  # if every member survived the CovFinder estimability candidate screen.
  for (pair in list(disease_scope_covs, mpd_subgroup_covs)) {
    if (any(pair %in% selected_covs)) {
      if (all(pair %in% eligible_covs)) {
        selected_covs <- unique(c(selected_covs, pair))
      } else {
        selected_covs <- setdiff(selected_covs, pair)
      }
    }
  }

  # Strategy B: only mathematical identifiability and full-rank constraints
  # are applied after CovFinder.
  minimal_limit <- max(0L, k_effective - 2L)

  drop_lowest_priority <- function(x) {
    if (length(x) == 0) return(x)

    # Drop the weakest selected group, not simply the last column.
    groups <- make_candidate_groups(x)
    strengths <- vapply(groups, function(g) {
      vals <- vapply(g, function(cv) cov_strength_score(dt, cv), numeric(1))
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) return(NA_real_)
      max(vals)
    }, numeric(1))

    if (all(!is.finite(strengths))) {
      victim <- tail(names(groups), 1)
    } else {
      tmp <- strengths
      tmp[!is.finite(tmp)] <- Inf
      victim <- names(which.min(tmp))[1]
    }

    setdiff(x, groups[[victim]])
  }

  while (length(selected_covs) > minimal_limit ||
         !design_is_full_rank(dt, selected_covs)) {
    old <- selected_covs
    selected_covs <- drop_lowest_priority(selected_covs)
    if (identical(old, selected_covs) || length(selected_covs) == 0) break
  }

  audit[, `:=`(
    selected = covariate %in% selected_covs,
    maximum_allowed_covariates = minimal_limit,
    k_effective = k_effective,
    covfinder_candidate_final = covariate %in% eligible_covs,
    covfinder_retry_dropped = covariate %in% covfinder_dropped,
    covfinder_retry_count = covfinder_retry_count,
    covfinder_last_error = covfinder_last_error
  )]

  audit[
    covariate %in% covfinder_dropped,
    exclusion_reason := "dropped_before_CovFinder_retry_for_estimability"
  ]

  audit[
    eligible == TRUE &
      selected == FALSE &
      !(covariate %in% covfinder_dropped),
    exclusion_reason := "not_selected_or_minimal_identifiability_rank"
  ]

  list(
    covariates = selected_covs,
    audit = audit,
    covfinder_status = covfinder_status,
    covfinder_retry_count = covfinder_retry_count,
    covfinder_dropped = covfinder_dropped,
    covfinder_last_error = covfinder_last_error
  )
}

fit_mrbrt_model <- function(dt, selected_covs, trm) {
  data <- mrtool$MRData()
  if (length(selected_covs) == 0) {
    data$load_df(
      data = dt,
      col_obs = "ln_effect",
      col_obs_se = "ln_se",
      col_study_id = "pmid"
    )
  } else {
    data$load_df(
      data = dt,
      col_obs = "ln_effect",
      col_obs_se = "ln_se",
      col_covs = as.list(selected_covs),
      col_study_id = "pmid"
    )
  }

  cov_models <- list(mrtool$LinearCovModel("intercept", use_re = TRUE))
  for (cov in selected_covs) {
    cov_models <- append(
      cov_models,
      list(do.call(mrtool$LinearCovModel, c(list(alt_cov = cov, use_re = FALSE))))
    )
  }

  model <- mrtool$MRBRT(
    data = data,
    cov_models = cov_models,
    inlier_pct = trm
  )
  model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)
  model
}

plot_forest <- function(forest_dt, summary_row, score_row, fig_path) {
  dt <- copy(forest_dt)[data == 1]
  dt <- dt[order(val)]
  dt[, row := seq_len(.N)]

  indicator_name <- if ("indicator" %in% names(summary_row)) {
    as.character(summary_row$indicator)
  } else {
    indicator
  }

  included_plot <- is.na(dt$included) | dt$included > 0
  point_size <- pmin(1.7, pmax(0.45, 0.15 / dt$se))

  draw_scaled <- function(scale = c("RR", "lnRR")) {
    scale <- match.arg(scale)
    bprf_curve_log <- if ("score_curve_log" %in% names(score_row) &&
                          is.finite(score_row$score_curve_log)) {
      score_row$score_curve_log
    } else {
      log(score_row$BPRF_OR)
    }

    if (scale == "RR") {
      x_val <- dt$val
      x_lo <- dt$lower
      x_hi <- dt$upper
      pooled <- exp(score_row$beta)
      pooled_lo <- exp(score_row$inner_lo)
      pooled_hi <- exp(score_row$inner_hi)
      bprf <- exp(bprf_curve_log)
      outer <- exp(c(score_row$outer_lo, score_row$outer_hi))
      inner <- exp(c(score_row$inner_lo, score_row$inner_hi))
      null_x <- 1
      xlab <- "Effect size (95% UI)"
      right_header <- "RR [95% CI]"
      value_label <- sprintf("%.2f [%.2f, %.2f]", x_val, x_lo, x_hi)
      log_axis <- "x"
    } else {
      x_val <- log(dt$val)
      x_lo <- log(dt$lower)
      x_hi <- log(dt$upper)
      pooled <- score_row$beta
      pooled_lo <- score_row$inner_lo
      pooled_hi <- score_row$inner_hi
      bprf <- bprf_curve_log
      outer <- c(score_row$outer_lo, score_row$outer_hi)
      inner <- c(score_row$inner_lo, score_row$inner_hi)
      null_x <- 0
      xlab <- "lnRR (95% UI)"
      right_header <- "lnRR [95% CI]"
      value_label <- sprintf("%.2f [%.2f, %.2f]", x_val, x_lo, x_hi)
      log_axis <- ""
    }

    if (any(!is.finite(outer))) outer <- c(bprf, pooled_hi)
    if (any(!is.finite(inner))) inner <- c(pooled_lo, pooled_hi)

    all_x <- c(x_lo, x_hi, pooled_lo, pooled_hi, bprf, outer, inner, null_x)
    data_min <- min(all_x, na.rm = TRUE)
    data_max <- max(all_x, na.rm = TRUE)

    if (scale == "RR") {
      plot_min <- max(0.03, data_min * 0.75)
      plot_max <- min(300, data_max * 1.10)
      xlim <- c(plot_min, plot_max * 2.15)
      label_x <- plot_max * 1.18
      tick_at <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)
      tick_at <- tick_at[tick_at >= xlim[1] & tick_at <= plot_max]
    } else {
      pad <- max(0.25, 0.12 * diff(range(all_x, na.rm = TRUE)))
      plot_min <- data_min - pad
      plot_max <- data_max + pad
      xlim <- c(plot_min, plot_max + 0.42 * (plot_max - plot_min))
      label_x <- plot_max + 0.06 * (plot_max - plot_min)
      tick_at <- pretty(c(plot_min, plot_max), n = 6)
    }

    ylim <- c(-2.1, nrow(dt) + 1.9)
    seg_lo <- pmax(x_lo, xlim[1])
    seg_hi <- pmin(x_hi, plot_max)
    # EXACT previous figure annotation:
    #   ROS_rounded + legacy star rating.
    # Do not substitute beta, BPRF_OR, or a newly derived score here.
    title_text <- paste0(
      indicator_name,
      ", ROS=",
      ifelse(
        is.finite(score_row$ROS_rounded),
        sprintf("%.4f", score_row$ROS_rounded),
        "NA"
      ),
      ", ",
      ifelse(
        !is.na(score_row$stars),
        star_label(score_row$stars),
        "Not estimable"
      )
    )

    draw_one <- function() {
      par(mar = c(6.2, 10.8, 3.0, 6.8), xpd = FALSE)
      plot(NA, xlim = xlim, ylim = ylim, log = log_axis, yaxt = "n",
           xaxt = "n", xlab = "", ylab = "", main = title_text,
           bty = "o")

      rect(max(min(outer), xlim[1]), ylim[1],
           min(max(outer), plot_max), ylim[2],
           col = adjustcolor("#c8e7ed", 0.65), border = NA)
      rect(max(min(inner), xlim[1]), ylim[1],
           min(max(inner), plot_max), ylim[2],
           col = adjustcolor("#8ecbd3", 0.70), border = NA)

      abline(v = null_x, lty = "dashed", lwd = 1.3, col = "gray25")
      abline(v = pooled, lwd = 2.2, col = "#008080")
      abline(v = bprf, lwd = 2.2, col = "#AA0C2D")

      segments(seg_lo, dt$row, seg_hi, dt$row,
               col = ifelse(included_plot, "black", "red"), lwd = 1.35)
      points(x_val[included_plot], dt$row[included_plot],
             pch = 19, cex = point_size[included_plot], col = "gray25")
      points(x_val[!included_plot], dt$row[!included_plot],
             pch = 4, cex = point_size[!included_plot] * 1.2,
             col = "red", lwd = 2.2)

      axis(1, at = tick_at, labels = if (scale == "RR") tick_at else round(tick_at, 2),
           cex.axis = 0.82)
      axis(2, at = dt$row,
           labels = dt$study_label,
           las = 2, cex.axis = 0.70, tick = FALSE)
      mtext(xlab, side = 1, line = 4.0, cex = 0.95)
      mtext("Study ID", side = 2, line = 7.5, cex = 0.95)

      text(label_x, nrow(dt) + 1.0, right_header, adj = 0,
           cex = 0.74, font = 2)
      text(label_x, dt$row, value_label, adj = 0, cex = 0.64)
    }

    suffix <- if (scale == "RR") "_forest_RR" else "_forest_lnRR"
    tiff(paste0(fig_path, suffix, ".tiff"), width = 1200, height = 900,
         units = "px", pointsize = 12, res = 100)
    draw_one()
    dev.off()

    pdf(paste0(fig_path, suffix, ".pdf"), width = 12, height = 9)
    draw_one()
    dev.off()
  }

  draw_scaled("RR")
  draw_scaled("lnRR")
}


# ---------------------------------------------------------------------------
# Define diagnostic membership without changing the scoring model
# ---------------------------------------------------------------------------
#
# Every stage is scored directly from its own MR-BRT fit. For Stage 4 this
# means beta, gamma, Fisher information, BPRF and ROS all come from the model
# fitted with inlier_pct=0.9. The repository weight threshold defines the
# diagnostic inlier set for Egger post and the post-trimming funnel plot;
# it never triggers a second model fit.
# ---------------------------------------------------------------------------

build_score_model <- function(model, original_dt, selected_covs, trm) {

  if (trm >= 0.999999) {
    return(list(
      score_model = model,
      score_dt = copy(original_dt),
      score_model_type = "stage_specific_MRBRT_direct_no_trimming",
      target_inliers = nrow(original_dt),
      requested_inliers = nrow(original_dt),
      trimming_count_matches_request = TRUE,
      selected_row_indices = seq_len(nrow(original_dt)),
      score_selected_covs = selected_covs,
      score_dropped_for_estimability = character(0)
    ))
  }

  w <- as.numeric(model$w_soln)
  if (length(w) != nrow(original_dt)) {
    stop(
      "Trimming weights are not aligned with the original analysis rows: ",
      "length(w_soln)=", length(w),
      ", nrow(data)=", nrow(original_dt)
    )
  }

  if (any(!is.finite(w))) {
    stop("Non-finite trimming weights were returned by MR-BRT.")
  }

  # IHME dichotomous diagnostic rule: the trimmed diagnostic set contains
  # observations with MR-BRT robust weight >= 0.1. Scoring still uses the
  # fitted stage model directly; this membership is for Egger/funnel/QC.
  requested_inliers <- max(
    2L,
    min(nrow(original_dt), as.integer(ceiling(trm * nrow(original_dt))))
  )

  keep_idx <- which(w >= 0.1)
  if (length(keep_idx) < 2L) {
    stop("Fewer than two observations retained by w_soln >= 0.1")
  }
  keep_idx <- sort(keep_idx)
  dt_inlier <- copy(original_dt[keep_idx])

  list(
    score_model = model,
    score_dt = dt_inlier,
    score_model_type = "stage_specific_MRBRT_direct_trim09",
    target_inliers = length(keep_idx),
    requested_inliers = requested_inliers,
    trimming_count_matches_request = length(keep_idx) == requested_inliers,
    selected_row_indices = keep_idx,
    score_selected_covs = selected_covs,
    score_dropped_for_estimability = character(0)
  )
}



# ---------------------------------------------------------------------------
# Safe MR-BRT fitting for numerically singular adjusted models
# ---------------------------------------------------------------------------
#
# CovFinder selection and ordinary QR full-rank screening are retained.
# If the actual MR-BRT optimizer still encounters a singular information
# matrix, the adjusted model is not numerically identifiable in that
# covariate configuration. We therefore remove the weakest selected
# covariate (or an indivisible reference-coded pair) and refit.
#
# This is an estimability safeguard, not significance-based model selection.
# Every removal is written to the covariate QC table.
# ---------------------------------------------------------------------------

drop_one_for_estimability <- function(selected_covs, dt) {
  selected_covs <- unique(selected_covs)
  if (length(selected_covs) == 0) {
    return(list(covariates = character(0), dropped = character(0)))
  }

  # Preserve reference-coded disease / MPD dummy sets as units when present.
  groups <- list()
  used <- character(0)

  for (pair in list(disease_scope_covs, mpd_subgroup_covs)) {
    present <- intersect(pair, selected_covs)
    if (length(present) > 0) {
      groups[[paste0("PAIR:", paste(pair, collapse = "+"))]] <- present
      used <- c(used, present)
    }
  }

  remaining <- setdiff(selected_covs, used)
  for (cv in remaining) groups[[cv]] <- cv

  # Group strength = maximum univariable weighted absolute association among
  # members. Drop the weakest group first.
  group_strength <- vapply(groups, function(g) {
    vals <- vapply(g, function(cv) cov_strength_score(dt, cv), numeric(1))
    max(vals, na.rm = TRUE)
  }, numeric(1))

  if (all(!is.finite(group_strength))) {
    victim_name <- tail(names(groups), 1)
  } else {
    group_strength[!is.finite(group_strength)] <- Inf
    victim_name <- names(which.min(group_strength))[1]
  }

  dropped <- groups[[victim_name]]
  kept <- setdiff(selected_covs, dropped)

  list(covariates = kept, dropped = dropped)
}


fit_mrbrt_model_with_estimability <- function(
    dt,
    selected_covs,
    trm,
    max_retries = NULL
) {
  selected_covs <- unique(selected_covs)
  if (is.null(max_retries)) max_retries <- length(selected_covs) + 1L

  current_covs <- selected_covs
  dropped_log <- list()
  attempt <- 0L
  last_error <- NA_character_

  repeat {
    attempt <- attempt + 1L

    fit <- tryCatch(
      fit_mrbrt_model(dt, current_covs, trm),
      error = function(e) e
    )

    if (!inherits(fit, "error")) {
      return(list(
        model = fit,
        selected_covs = current_covs,
        dropped_for_estimability = unique(unlist(dropped_log)),
        n_fit_attempts = attempt,
        last_fit_error = last_error
      ))
    }

    last_error <- conditionMessage(fit)

    singular_like <- grepl(
      "singular|LinAlgError|not positive definite|ill-conditioned|rank",
      last_error,
      ignore.case = TRUE
    )

    if (!singular_like) {
      stop(fit)
    }

    if (length(current_covs) == 0 || attempt >= max_retries) {
      stop(
        "MR-BRT remained non-estimable after covariate pruning. Last error: ",
        last_error
      )
    }

    pruned <- drop_one_for_estimability(current_covs, dt)

    if (identical(pruned$covariates, current_covs) ||
        length(pruned$dropped) == 0) {
      stop(
        "Unable to resolve singular MR-BRT design. Last error: ",
        last_error
      )
    }

    message(
      "  MR-BRT singular adjusted model; dropping for estimability: ",
      paste(pruned$dropped, collapse = ", "),
      " | previous error: ", last_error
    )

    dropped_log[[length(dropped_log) + 1L]] <- pruned$dropped
    current_covs <- pruned$covariates
  }
}


# ---------------------------------------------------------------------------
# MR-BRT-adjusted residual diagnostics and Egger-type intercept
# ---------------------------------------------------------------------------

make_adjusted_diagnostics <- function(
    model, dt, selected_covs, scenario_id, retained_row_indices
) {

  # Fixed-effect prediction is essential here. We remove the fitted signal and
  # selected study-level bias covariates, but DO NOT absorb study random effects
  # into each study's own fitted value.
  fixed_pred <- as.numeric(model$predict(
    data = model$data,
    predict_for_study = FALSE,
    sort_by_data_id = TRUE
  ))

  obs <- as.numeric(model$data$obs)
  obs_se <- as.numeric(model$data$obs_se)
  gamma <- max(0, as.numeric(model$gamma_soln)[1])
  weights <- as.numeric(model$w_soln)

  if (length(obs) != nrow(dt) || length(fixed_pred) != nrow(dt)) {
    stop("MR-BRT diagnostics row alignment failed.")
  }

  out <- copy(dt)
  out[, scenario_id := scenario_id]
  out[, model_obs := obs]
  out[, model_obs_se := obs_se]
  out[, fixed_effect_prediction := fixed_pred]
  out[, residual := model_obs - fixed_effect_prediction]
  # Match the public IHME dichotomous implementation: residual uncertainty is
  # the overall residual SD sqrt(reported sampling SE^2 + fitted gamma).
  out[, reported_se := model_obs_se]
  out[, residual_sd := sqrt(model_obs_se^2 + gamma)]
  out[, standardized_residual := residual / residual_sd]
  out[, model_weight := weights]
  out[, original_row_index := .I]
  out[, retained_for_final := original_row_index %in% retained_row_indices]
  out[, trimmed_outlier := !retained_for_final]
  out[, selected_covariates := paste(selected_covs, collapse = ";")]

  out
}


adjusted_egger_intercept <- function(diag_dt, dataset_label) {

  d <- copy(diag_dt)
  if (dataset_label == "post_trimming") {
    d <- d[trimmed_outlier == FALSE]
  }

  d <- d[
    is.finite(residual) &
      is.finite(residual_sd) &
      residual_sd > 0
  ]

  k <- uniqueN(d$pmid)
  n <- nrow(d)

  ans <- data.table(
    dataset = dataset_label,
    n_observations = n,
    k_unique_studies = k,
    egger_intercept = NA_real_,
    egger_intercept_se = NA_real_,
    test_statistic = NA_real_,
    p_value = NA_real_,
    formally_interpretable_k_ge_10 = k >= MIN_K_EGGER,
    small_study_effect_signal = NA
  )

  if (n < 3) return(ans)

  # Egger-type test after removing MR-BRT signal/bias covariates.
  # standardized residual ~ precision
  z <- d$residual / d$residual_sd
  precision <- 1 / d$residual_sd

  fit <- lm(z ~ precision)
  cc <- summary(fit)$coefficients

  ans[, `:=`(
    egger_intercept = unname(cc["(Intercept)", "Estimate"]),
    egger_intercept_se = unname(cc["(Intercept)", "Std. Error"]),
    test_statistic = unname(cc["(Intercept)", "t value"]),
    p_value = unname(cc["(Intercept)", "Pr(>|t|)"])
  )]

  if (k >= MIN_K_EGGER) {
    ans[, small_study_effect_signal := p_value < EGGER_ALPHA]
  }

  ans
}


ihme_repository_residual_test <- function(diag_dt, dataset_label) {
  d <- copy(diag_dt)
  if (dataset_label == "post_trimming") {
    d <- d[trimmed_outlier == FALSE]
  }

  x <- d$standardized_residual[
    is.finite(d$standardized_residual)
  ]

  n <- length(x)
  k <- uniqueN(d$pmid)

  if (n < 2) {
    return(data.table(
      dataset = dataset_label,
      n_observations = n,
      k_unique_studies = k,
      mean_standardized_residual = NA_real_,
      se_mean = NA_real_,
      z_abs = NA_real_,
      p_one_sided = NA_real_,
      formally_interpretable_k_ge_10 = k >= MIN_K_EGGER
    ))
  }

  r_mean <- mean(x)
  r_sd <- 1 / sqrt(n)
  z_abs <- abs(r_mean / r_sd)
  p <- 1 - pnorm(z_abs)

  data.table(
    dataset = dataset_label,
    n_observations = n,
    k_unique_studies = k,
    mean_standardized_residual = r_mean,
    se_mean = r_sd,
    z_abs = z_abs,
    p_one_sided = p,
    formally_interpretable_k_ge_10 = k >= MIN_K_EGGER
  )
}


# ---------------------------------------------------------------------------
# Consecutive sampling vs study-design dependence check
# ---------------------------------------------------------------------------

assess_consecutive_vs_design <- function(dt) {
  cons_names <- c("cv_consecutive_sampling", "cv_consecutive")
  design_names <- c("cv_studydesign", "cv_study_design")

  cons <- intersect(cons_names, names(dt))[1]
  design <- intersect(design_names, names(dt))[1]

  if (is.na(cons) || is.na(design)) {
    return(data.table(
      consecutive_covariate = ifelse(is.na(cons), NA_character_, cons),
      study_design_covariate = ifelse(is.na(design), NA_character_, design),
      n_complete = NA_integer_,
      cramers_v = NA_real_,
      fisher_p = NA_real_,
      interpretation = "One or both covariates were not available."
    ))
  }

  x <- dt[[cons]]
  y <- dt[[design]]
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 4 || uniqueN(x) < 2 || uniqueN(y) < 2) {
    return(data.table(
      consecutive_covariate = cons,
      study_design_covariate = design,
      n_complete = length(x),
      cramers_v = NA_real_,
      fisher_p = NA_real_,
      interpretation = "Insufficient variation for dependence assessment."
    ))
  }

  tab <- table(x, y)
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  denom <- min(nrow(tab) - 1, ncol(tab) - 1)
  v <- if (denom > 0) sqrt(as.numeric(chi$statistic) / sum(tab) / denom) else NA_real_
  fp <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)

  interpretation <- if (!is.finite(v)) {
    "Not estimable"
  } else if (v >= 0.70) {
    "Strong dependence; avoid automatic simultaneous inclusion."
  } else if (v >= 0.50) {
    "Moderate-to-strong dependence; inspect sparsity/collinearity."
  } else {
    "Low-to-moderate dependence; retain as conceptually distinct candidates."
  }

  data.table(
    consecutive_covariate = cons,
    study_design_covariate = design,
    n_complete = sum(tab),
    cramers_v = v,
    fisher_p = fp,
    interpretation = interpretation
  )
}


plot_funnel <- function(diag_dt, fig_path, dataset_label) {
  dt <- copy(diag_dt)
  if (identical(dataset_label, "post_trimming")) {
    dt <- dt[retained_for_final == TRUE]
  }
  max_se <- quantile(dt$residual_sd, 0.99, na.rm = TRUE)
  if (!is.finite(max_se) || max_se <= 0) max_se <- max(dt$residual_sd, na.rm = TRUE)

  draw_one <- function() {
    par(mar = c(5, 5, 4, 2))
    plot(dt$residual, dt$residual_sd,
         pch = ifelse(dt$retained_for_final, 19, 4),
         col = ifelse(dt$retained_for_final, "gray25", "red"),
         cex = pmin(1.8, pmax(0.45, 0.12 / dt$residual_sd)),
         lwd = ifelse(dt$retained_for_final, 1, 2.5),
         xlim = c(-2.5 * max_se, 2.5 * max_se),
         ylim = c(max_se, 0),
         xlab = "Fixed-effect residual log effect", ylab = "Overall residual SD",
         main = paste(unique(dt$scenario_id), dataset_label))
    polygon(c(0, -1.96 * max_se, 1.96 * max_se),
            c(0, max_se, max_se),
            col = adjustcolor("#c8e7ed", 0.45), border = NA)
    abline(v = 0, lty = "dashed", lwd = 1.4, col = "black")
    lines(c(0, -1.96 * max_se), c(0, max_se), col = "#008080", lwd = 1.5)
    lines(c(0, 1.96 * max_se), c(0, max_se), col = "#008080", lwd = 1.5)
    points(dt$residual, dt$residual_sd,
           pch = ifelse(dt$retained_for_final, 19, 4),
           col = ifelse(dt$retained_for_final, "gray25", "red"),
           cex = pmin(1.8, pmax(0.45, 0.12 / dt$residual_sd)),
           lwd = ifelse(dt$retained_for_final, 1, 2.5))
  }

  tiff(paste0(fig_path, "_funnel_", dataset_label, ".tiff"), width = 1000, height = 850,
       units = "px", pointsize = 12, res = 100)
  draw_one()
  dev.off()

  pdf(paste0(fig_path, "_funnel_", dataset_label, ".pdf"), width = 10, height = 8.5)
  draw_one()
  dev.off()
}

run_scenario <- function(dt, scenario, locked_covariates = NULL) {
  # -------------------------------------------------------------------------
  # MODEL-STAGE DISPATCH
  # The same execution pipeline is used for both adjusted stages. The scenario row
  # determines whether CovFinder is used and whether trimming is applied.
  # -------------------------------------------------------------------------
  scenario_id <- scenario$scenario_id
  model_stage <- scenario$model_stage
  model_stage_label <- scenario$model_stage_label
  covariate_mode <- scenario$covariate_mode
  scenario_dir <- file.path(save_dir, indicator, scenario_id)
  scenario_fig_dir <- file.path(fig_dir, indicator, scenario_id)
  dir.create(scenario_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(scenario_fig_dir, showWarnings = FALSE, recursive = TRUE)

  out_base <- file.path(scenario_dir, paste0(indicator, "_", scenario_id))
  fig_base <- file.path(scenario_fig_dir, paste0(indicator, "_", scenario_id))

  message("\nRunning Stage ", model_stage, " / ", scenario_id,
          ": ", scenario$description)

  selection <- if (identical(covariate_mode, "discover_by_stage")) {
    select_covariates(dt, scenario$trm)
  } else {
    audit <- screen_covariates(
      dt,
      candidate_covariates_for_indicator(dt, indicator)
    )
    audit[, `:=`(
      selected = FALSE,
      maximum_allowed_covariates = 0L,
      k_effective = max(0L, floor(nrow(dt) * scenario$trm))
    )]
    audit[eligible == TRUE, exclusion_reason := "intercept_only_scenario"]
    list(
      covariates = character(0),
      audit = audit,
      covfinder_status = "not_applicable_intercept_only_scenario",
      covfinder_retry_count = 0L,
      covfinder_dropped = character(0),
      covfinder_last_error = NA_character_
    )
  }
  selected_covs <- selection$covariates
  k_effective <- max(0L, floor(nrow(dt) * scenario$trm))

  covariate_qc <- copy(selection$audit)
  covariate_qc[, `:=`(
    indicator = indicator,
    model_stage = model_stage,
    model_stage_label = model_stage_label,
    scenario_id = scenario_id,
    n_rows_original = nrow(dt),
    candidate_covariate_count = nrow(selection$audit),
    eligible_covariate_count = sum(selection$audit$eligible),
    selected_covariate_count = length(selected_covs),
    covariate_mode = covariate_mode,
    covariate_source = fifelse(
      covariate_mode == "discover_by_stage",
      paste0("CovFinder_same_stage_inlier_pct_", scenario$trm),
      "none"
    ),
    covfinder_status = selection$covfinder_status,
    covfinder_retry_count = selection$covfinder_retry_count,
    covfinder_dropped = paste(selection$covfinder_dropped, collapse = ";")
  )]
  fwrite(covariate_qc, paste0(out_base, "_covariate_QC.csv"))

  message("  Selected covariates: ",
          ifelse(length(selected_covs) == 0, "none", paste(selected_covs, collapse = ", ")))

  if (identical(
    selection$covfinder_status,
    "covfinder_not_estimable_after_candidate_pruning"
  )) {
    stop(
      "CovFinder was not estimable for ", scenario_id, ": ",
      selection$covfinder_last_error
    )
  }

  # FITTING STEP: estimate beta, selected fixed effects, gamma and robust
  # observation weights under the current stage specification.
  fit_bundle <- fit_mrbrt_model_with_estimability(
    dt = dt,
    selected_covs = selected_covs,
    trm = scenario$trm
  )
  model <- fit_bundle$model

  # The final estimable covariate set may be smaller than the raw CovFinder
  # selection if the latter produces a singular MR-BRT information matrix.
  initial_selected_covs <- selected_covs
  selected_covs <- fit_bundle$selected_covs

  if (length(fit_bundle$dropped_for_estimability) > 0) {
    covariate_qc[
      covariate %in% fit_bundle$dropped_for_estimability,
      `:=`(
        selected = FALSE,
        exclusion_reason = "dropped_after_CovFinder_for_MRBRT_estimability"
      )
    ]
  }

  covariate_qc[, `:=`(
    selected_after_estimability_check = covariate %in% selected_covs,
    n_fit_attempts = fit_bundle$n_fit_attempts,
    fit_error_resolved = ifelse(
      is.na(fit_bundle$last_fit_error),
      NA_character_,
      fit_bundle$last_fit_error
    )
  )]

  fwrite(covariate_qc, paste0(out_base, "_covariate_QC.csv"))
  design_matrix <- if (length(selected_covs) == 0) {
    matrix(1, nrow = nrow(dt), ncol = 1)
  } else {
    cbind(intercept = 1, as.matrix(dt[, ..selected_covs]))
  }
  covariate_qc[, `:=`(
    model_rank = qr(design_matrix)$rank,
    model_columns = ncol(design_matrix),
    residual_df_after_beta_and_gamma =
      max(0L, k_effective - ncol(design_matrix) - 1L)
  )]
  fwrite(covariate_qc, paste0(out_base, "_covariate_QC.csv"))

  # MRTool binary documentation states that `intercept` is added
  # automatically by MRData. For the adjusted reference study (all selected
  # bias covariates = 0), the target log effect is therefore the intercept.
  # We take it from the BPRF/Scorelator result below rather than constructing a
  # separate prediction MRData object that redundantly passes `intercept`.

  data_model <- model$data
  obs_pred <- as.numeric(model$predict(
    data = data_model,
    predict_for_study = TRUE,
    sort_by_data_id = TRUE
  ))

  analysis_dt <- copy(dt)
  analysis_dt[, scenario_id := scenario_id]
  analysis_dt[, model_obs := as.numeric(data_model$obs)]
  analysis_dt[, model_obs_se := as.numeric(data_model$obs_se)]
  analysis_dt[, model_pred := obs_pred]
  analysis_dt[, model_residual := model_obs - model_pred]
  analysis_dt[, model_weight := as.numeric(model$w_soln)]
  analysis_dt[, included := model_weight > 0]
  analysis_dt[, selected_covs := paste(selected_covs, collapse = ";")]

  # Construct the single authoritative diagnostic membership. The same rows
  # drive Egger post and the post-trimming funnel plot. ROS/BPRF are calculated
  # directly from the stage-specific MR-BRT model.
  score_bundle <- build_score_model(
    model = model,
    original_dt = dt,
    selected_covs = selected_covs,
    trm = scenario$trm
  )
  score_model <- score_bundle$score_model

  # Small-study-effect / reporting-bias diagnostics.
  # We deliberately calculate BOTH pre-trimming and post-trimming results.
  diag_dt <- make_adjusted_diagnostics(
    model = model,
    dt = dt,
    selected_covs = selected_covs,
    scenario_id = scenario_id,
    retained_row_indices = score_bundle$selected_row_indices
  )

  adjusted_egger_dt <- rbindlist(list(
    adjusted_egger_intercept(diag_dt, "pre_trimming"),
    adjusted_egger_intercept(diag_dt, "post_trimming")
  ), fill = TRUE)
  adjusted_egger_dt[, `:=`(
    indicator = indicator,
    scenario_id = scenario_id
  )]

  ihme_sse_dt <- rbindlist(list(
    ihme_repository_residual_test(diag_dt, "pre_trimming"),
    ihme_repository_residual_test(diag_dt, "post_trimming")
  ), fill = TRUE)
  ihme_sse_dt[, `:=`(
    indicator = indicator,
    scenario_id = scenario_id
  )]

  fwrite(diag_dt, paste0(out_base, "_small_study_diagnostics.csv"))
  fwrite(adjusted_egger_dt, paste0(out_base, "_adjusted_Egger_pre_post.csv"))
  fwrite(ihme_sse_dt, paste0(out_base, "_IHME_residual_small_study_test.csv"))

  forest_dt <- data.table(
    scenario_id = scenario_id,
    val = exp(analysis_dt$model_obs),
    se = analysis_dt$model_obs_se,
    study = as.character(analysis_dt$pmid),
    included = analysis_dt$model_weight,
    study_label = analysis_dt$study_label,
    data = 1L
  )
  for (cov in selected_covs) forest_dt[, (cov) := analysis_dt[[cov]]]
  forest_dt[, lower := exp(log(val) - 1.96 * se)]
  forest_dt[, upper := exp(log(val) + 1.96 * se)]
  forest_dt <- forest_dt[order(val)]
  forest_dt[included > 0 & included < 1, included := 0.5]

  # Score the stage-specific model directly. Stage 4 therefore retains beta,
  # gamma and Fisher information from MR-BRT at inlier_pct=0.9. If the official
  # Scorelator cannot process a trimmed object, the Fisher fallbacks operate on
  # this same model; no post-trimming refit is permitted.
  score_dt <- ihme_bprf_score(
    score_model,
    indicator,
    selected_covs = score_bundle$score_selected_covs
  )

  # Hard stop if the values used in tables and figures diverge from the
  # previous validated scoring convention.
  if (grepl("^ok", score_dt$score_status[1])) {
    validate_legacy_score_fields(score_dt)
  }

  score_dt[, scenario_id := scenario_id]
  score_dt[, `:=`(
    n_score_rows = nrow(score_model$data$to_df()),
    score_model_type = score_bundle$score_model_type,
    score_target_inliers = score_bundle$target_inliers,
    score_requested_inliers = score_bundle$requested_inliers,
    trimming_count_matches_request =
      score_bundle$trimming_count_matches_request,
    score_selected_covs = paste(score_bundle$score_selected_covs, collapse = ";"),
    score_dropped_for_estimability = paste(
      score_bundle$score_dropped_for_estimability,
      collapse = ";"
    )
  )]

  # Save the repository-threshold membership used for post-trimming diagnostics.
  score_membership <- copy(dt)
  score_membership[, original_row_index := .I]
  score_membership[, robust_weight := as.numeric(model$w_soln)]
  score_membership[, retained_for_posttrim_diagnostics :=
    original_row_index %in% score_bundle$selected_row_indices]
  fwrite(
    score_membership,
    paste0(out_base, "_posttrim_diagnostic_membership.csv")
  )

  covariate_qc[, `:=`(
    gamma = score_dt$gamma[1],
    gamma_var = score_dt$gamma_var[1],
    gamma_status = score_dt$gamma_status[1],
    gamma_raw = score_dt$gamma_raw[1],
    gamma_fisher_boundary_used = score_dt$gamma_fisher_boundary_used[1],
    gamma_fisher_epsilon = score_dt$gamma_fisher_epsilon[1],
    scorelator_trigger_error = score_dt$scorelator_trigger_error[1],
    formula_fallback_error = score_dt$formula_fallback_error[1],
    bprf_score_status = score_dt$score_status[1]
  )]
  fwrite(covariate_qc, paste0(out_base, "_covariate_QC.csv"))

  if (!grepl("^ok", score_dt$score_status[1])) {
    fwrite(score_dt, paste0(out_base, "_BPRF_NOT_ESTIMABLE_QC.csv"))
    warning(
      "IHME BPRF/ROS not estimable for ", indicator, " / ", scenario_id,
      ": ", score_dt$score_error[1],
      ". gamma uncertainty was not replaced by zero.",
      call. = FALSE
    )
    return(cbind(
      data.table(
        indicator = indicator,
        model_stage = model_stage,
        model_stage_label = model_stage_label,
        scenario_id = scenario_id,
        scenario_description = scenario$description,
        trm = scenario$trm,
        use_covfinder = scenario$use_covfinder,
        covariate_mode = covariate_mode,
        n_rows = nrow(dt),
        n_included_weight_gt0 = sum(as.numeric(model$w_soln) > 0),
        selected_covs = paste(selected_covs, collapse = ";"),
        selected_cov_count = length(selected_covs)
      ),
      score_dt
    ))
  }
  score_row <- score_dt[1]

  result_dt <- data.table(
    scenario_id = scenario_id,
    val = exp(score_row$beta),
    se = NA_real_,
    study = "IHME BPRF score model",
    included = 1,
    study_label = "IHME BPRF score model",
    data = 2L,
    lower = score_row$outer_OR_lo,
    upper = score_row$outer_OR_hi,
    selected_covs = paste(selected_covs, collapse = ";"),
    prediction_covariate_values = paste(
      if (length(selected_covs) == 0) "none" else paste0(selected_covs, "=0"),
      collapse = ";"
    )
  )
  for (cov in selected_covs) result_dt[, (cov) := 0]
  forest_dt <- rbindlist(list(result_dt, forest_dt), fill = TRUE)
  forest_dt[, row := seq_len(.N)]

  summary_dt <- data.table(
    indicator = indicator,
    model_stage = model_stage,
    model_stage_label = model_stage_label,
    scenario_id = scenario_id,
    scenario_description = scenario$description,
    trm = scenario$trm,
    use_covfinder = scenario$use_covfinder,
    covariate_mode = covariate_mode,
    n_rows = nrow(dt),
    n_included_weight_gt0 = sum(as.numeric(model$w_soln) > 0),
    selected_covs_initial = paste(initial_selected_covs, collapse = ";"),
    selected_covs = paste(selected_covs, collapse = ";"),
    selected_cov_count = length(selected_covs),
    dropped_for_estimability = paste(
      fit_bundle$dropped_for_estimability,
      collapse = ";"
    ),
    effect_model = ifelse(
      length(selected_covs) == 0,
      "intercept_random_effect_only",
      "intercept_random_effect_plus_selected_cv_fixed_effects"
    ),
    prediction_covariate_values = paste(
      if (length(selected_covs) == 0) "none" else paste0(selected_covs, "=0"),
      collapse = ";"
    ),
    adjusted_model_log_prediction = score_row$beta,
    gamma_main = max(0, as.numeric(model$gamma_soln)[1]),
    score_model_type = score_bundle$score_model_type,
    score_target_inliers = score_bundle$target_inliers,
    score_requested_inliers = score_bundle$requested_inliers,
    trimming_count_matches_request =
      score_bundle$trimming_count_matches_request,
    score_model_n_rows = nrow(score_model$data$to_df()),
    retained_row_indices = paste(
      sort(score_bundle$selected_row_indices), collapse = ";"
    ),
    trimmed_row_indices = paste(
      setdiff(seq_len(nrow(dt)), score_bundle$selected_row_indices),
      collapse = ";"
    ),
    covariate_lock_preserved = NA,
    score_covariate_lock_preserved = NA,

    egger_pre_intercept = adjusted_egger_dt[
      dataset == "pre_trimming", egger_intercept
    ][1],
    egger_pre_p = adjusted_egger_dt[
      dataset == "pre_trimming", p_value
    ][1],
    egger_post_intercept = adjusted_egger_dt[
      dataset == "post_trimming", egger_intercept
    ][1],
    egger_post_p = adjusted_egger_dt[
      dataset == "post_trimming", p_value
    ][1],

    reporting_interval_type = paste(
      "no_gamma_95pct_CI_beta_only;",
      "gamma_point_only_95pct_UI_diagnostic_not_for_manuscript;",
      "with_gamma_95pct_UI_beta_plus_gamma_upper"
    ),
    score_interval_type = "IHME_BPRF_tail_5pct_or_95pct_gamma_plus_2SD"
  )

  summary_with_score <- cbind(
    summary_dt[rep(1, nrow(score_dt))],
    score_dt
  )

  # Explicit reporting aliases.
  #
  # Published BPRF terminology:
  #   - "without gamma" = beta uncertainty only;
  #   - "with gamma" = beta uncertainty plus the upper uncertainty bound of
  #     gamma (gamma_ub), i.e. the conservative outer interval.
  #
  # The point-gamma/inner interval is diagnostic only. It is NOT the published
  # "with gamma" interval because it ignores uncertainty in gamma. If the
  # fitted gamma is zero, this diagnostic interval correctly equals the
  # no-gamma interval, while the published with-gamma interval remains wider
  # whenever gamma_var > 0.
  summary_with_score[, `:=`(
    MRBRT_no_gamma_CI95_lo = no_gamma_OR_lo,
    MRBRT_no_gamma_CI95_hi = no_gamma_OR_hi,
    MRBRT_gamma_point_only_UI95_lo = inner_OR_lo,
    MRBRT_gamma_point_only_UI95_hi = inner_OR_hi,
    MRBRT_with_gamma_UI95_lo = outer_OR_lo,
    MRBRT_with_gamma_UI95_hi = outer_OR_hi,
    BPRF_conservative_UI95_lo = outer_OR_lo,
    BPRF_conservative_UI95_hi = outer_OR_hi,
    legacy_score_log = score_log,
    legacy_BPRF_OR_for_scoring = BPRF_OR,
    legacy_score_curve_log_for_plot = score_curve_log,
    legacy_BPRF_curve_OR_for_plot = BPRF_curve_OR,
    legacy_ROS = ROS,
    legacy_ROS_rounded_figure_label = ROS_rounded,
    legacy_stars = stars
  )]

  # Interval and gamma provenance audit. This prevents the point-gamma
  # diagnostic interval from being silently reported as the published
  # with-gamma UI and documents that scoring used the same stage model.
  interval_tol <- 1e-10
  summary_with_score[, `:=`(
    no_gamma_width_log = no_gamma_hi - no_gamma_lo,
    gamma_point_only_width_log = inner_hi - inner_lo,
    with_gamma_width_log = outer_hi - outer_lo,
    with_gamma_equal_no_gamma =
      abs(outer_lo - no_gamma_lo) <= interval_tol &
      abs(outer_hi - no_gamma_hi) <= interval_tol,
    point_gamma_equal_no_gamma =
      abs(inner_lo - no_gamma_lo) <= interval_tol &
      abs(inner_hi - no_gamma_hi) <= interval_tol,
    scenario_model_gamma = gamma_main,
    score_model_gamma = gamma,
    score_refit_changed_gamma =
      is.finite(gamma_main) & is.finite(gamma) &
      abs(gamma - gamma_main) > 1e-8
  )]

  summary_with_score[, with_gamma_interval_status := fifelse(
    !is.finite(gamma_var) | gamma_var < 0,
    "ERROR_gamma_uncertainty_not_estimable",
    fifelse(
      with_gamma_equal_no_gamma & gamma_ub > interval_tol,
      "ERROR_with_gamma_alias_or_formula_did_not_use_gamma_ub",
      fifelse(
        with_gamma_equal_no_gamma,
        "equal_only_because_gamma_ub_is_numerically_zero",
        "ok_with_gamma_wider_than_without_gamma"
      )
    )
  )]

  if (any(summary_with_score$with_gamma_interval_status ==
          "ERROR_with_gamma_alias_or_formula_did_not_use_gamma_ub")) {
    stop(
      "With-gamma UI equals no-gamma UI despite positive gamma_ub. ",
      "The reporting alias or interval formula is incorrect."
    )
  }

  fwrite(analysis_dt, paste0(out_base, "_analysis_data.csv"))
  write.csv(forest_dt, paste0(out_base, "_forest_data.csv"), row.names = FALSE)
  fwrite(summary_with_score, paste0(out_base, "_bprf_summary.csv"))
  save_dill(model, paste0(out_base, "_model.pkl"))
  save_dill(score_model, paste0(out_base, "_score_model.pkl"))

  plot_forest(forest_dt, summary_dt[1], score_row, fig_base)
  plot_funnel(diag_dt, fig_base, "pre_trimming")
  plot_funnel(diag_dt, fig_base, "post_trimming")

  summary_with_score
}


# ---------------------------------------------------------------------------
# BATCH EXECUTION
# ---------------------------------------------------------------------------

run_one_indicator_file <- function(input_file) {

  indicator_name <- sub("_all_extraction\\.csv$", "", basename(input_file))
  indicator <<- indicator_name

  message("\n============================================================")
  message("RUNNING INDICATOR: ", indicator)
  message("Input: ", input_file)
  message("============================================================")

  dt_sub <- fread(input_file)

  required_cols <- c("ln_effect", "ln_se", "pmid")
  missing_cols <- setdiff(required_cols, names(dt_sub))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  dt_sub <- dt_sub[
    is.finite(ln_effect) &
      is.finite(ln_se) &
      ln_se > 0 &
      !is.na(pmid)
  ]

  if (nrow(dt_sub) < 2) {
    indicator_out_dir <- file.path(save_dir, indicator)
    dir.create(indicator_out_dir, recursive = TRUE, showWarnings = FALSE)

    insufficient_dt <- data.table(
      indicator = indicator,
      scenario_id = NA_character_,
      n_usable_effect_rows = nrow(dt_sub),
      score_status = "insufficient_studies_for_MRBRT",
      score_error = "Fewer than two usable study-level effect estimates"
    )

    fwrite(
      insufficient_dt,
      file.path(
        indicator_out_dir,
        paste0(indicator, "_INSUFFICIENT_STUDIES.csv")
      )
    )

    warning(
      "Skipping ", indicator,
      ": fewer than two usable study-level effect rows for MR-BRT.",
      call. = FALSE
    )

    return(insufficient_dt)
  }

  dt_sub[, study_label := make_study_label(dt_sub)]

  # Drop zero-variance candidate covariates before CovFinder.
  cv_cols_all <- names(dt_sub)[grepl("^cv_", names(dt_sub))]
  zero_var_cv <- cv_cols_all[vapply(cv_cols_all, function(col) {
    vals <- dt_sub[[col]]
    length(unique(vals[!is.na(vals)])) <= 1
  }, logical(1))]

  if (length(zero_var_cv) > 0) {
    message(
      "Dropping zero-variance cv_ columns: ",
      paste(zero_var_cv, collapse = ", ")
    )
    dt_sub[, (zero_var_cv) := NULL]
  }

  # Master covariate audit
  all_cv_now <- names(dt_sub)[grepl("^cv_", names(dt_sub))]
  candidate_before_exclusion <- candidate_covariates_for_indicator(dt_sub, indicator)
  candidate_after_exclusion <- exclude_unwanted_covariates(candidate_before_exclusion)

  covariate_master_audit <- data.table(
    covariate = sort(unique(all_cv_now))
  )

  covariate_master_audit[, in_candidate_family :=
    covariate %in% candidate_before_exclusion]

  covariate_master_audit[, allowed_after_exclusion :=
    covariate %in% candidate_after_exclusion]

  covariate_master_audit[, exclusion_reason := fifelse(
    covariate %in% imaging_quality_covs & !is_imaging_indicator(indicator),
    "excluded_reader_quality_for_nonimaging_indicator",
    fifelse(
    grepl(
      "surgical_criteria|surgery_criteria|surgical_indication|surgery_indication|resection_criteria|resection_indication|absolute_indication|operation_criteria|operative_indication",
      covariate,
      ignore.case = TRUE
    ),
    "excluded_surgical_selection_mechanism",
    fifelse(
      grepl("^cov_adj", covariate, ignore.case = TRUE) |
        tolower(covariate) %in% c("cv_adj_l0", "cv_outcome", "cov_outcome"),
      "excluded_by_prespecified_rule",
      fifelse(
        allowed_after_exclusion,
        "eligible_candidate_family",
        "not_in_prespecified_candidate_family"
      )
    )
    )
  )]

  indicator_out_dir <- file.path(save_dir, indicator)
  dir.create(indicator_out_dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(
    covariate_master_audit,
    file.path(
      indicator_out_dir,
      paste0(indicator, "_MASTER_covariate_audit.csv")
    )
  )

  # Consecutive vs study-design dependence
  dependence_dt <- assess_consecutive_vs_design(dt_sub)
  fwrite(
    dependence_dt,
    file.path(
      indicator_out_dir,
      paste0(indicator, "_consecutive_vs_studydesign.csv")
    )
  )

  # Run the two adjusted model stages in order:
  #   Stage 3 = CovFinder(1.0) + adjusted MR-BRT(1.0)
  #   Stage 4 = CovFinder(0.9) + adjusted MR-BRT(0.9), primary
  # This mirrors the public IHME dichotomous script: CovFinder and the final
  # model use the same inlier_pct within each independently fitted pipeline.
  scenario_results <- vector("list", nrow(scenarios))

  for (i in seq_len(nrow(scenarios))) {

    scenario_id_i <- scenarios$scenario_id[i]
    model_stage_i <- scenarios$model_stage[i]
    model_stage_label_i <- scenarios$model_stage_label[i]

    message("\n>>> Starting Stage ", model_stage_i, ": ", model_stage_label_i)

    scenario_results[[i]] <- tryCatch(
      {
        result_i <- run_scenario(
          dt_sub,
          scenarios[i],
          locked_covariates = NULL
        )

        result_i
      },
      error = function(e) {

        err_dt <- data.table(
          indicator = indicator,
          model_stage = model_stage_i,
          model_stage_label = model_stage_label_i,
          scenario_id = scenario_id_i,
          scenario_description = scenarios$description[i],
          trm = scenarios$trm[i],
          use_covfinder = scenarios$use_covfinder[i],
          score_status = "scenario_failed",
          score_error = conditionMessage(e)
        )

        fwrite(
          err_dt,
          file.path(
            indicator_out_dir,
            paste0(
              indicator, "_", scenario_id_i,
              "_SCENARIO_ERROR.csv"
            )
          )
        )

        warning(
          "Scenario failed for ",
          indicator, " / ", scenario_id_i,
          ": ", conditionMessage(e),
          call. = FALSE
        )

        err_dt
      }
    )
  }

  all_summary <- rbindlist(scenario_results, fill = TRUE)

  # -----------------------------------------------------------------------
  # Prespecified trimming-effect audit
  #   S4 - S3: complete adjusted-pipeline difference after enabling trimming
  # CovFinder is rerun at each stage, so S3/S4 may have different covariates.
  # The comparison is not a same-covariate pure trimming contrast.
  # -----------------------------------------------------------------------
  parse_index_set <- function(x) {
    if (length(x) == 0 || is.na(x) || !nzchar(x)) return(integer(0))
    as.integer(strsplit(x, ";", fixed = TRUE)[[1]])
  }

  get_stage_row <- function(stage_id) {
    z <- all_summary[scenario_id == stage_id]
    if (nrow(z) == 0 || !("score_status" %in% names(z)) ||
        !grepl("^ok", z$score_status[1])) return(NULL)
    z[1]
  }

  make_trim_comparison <- function(reference_id, trimmed_id, comparison_name) {
    ref <- get_stage_row(reference_id)
    trm <- get_stage_row(trimmed_id)
    if (is.null(ref) || is.null(trm)) return(NULL)

    ref_set <- parse_index_set(ref$retained_row_indices)
    trm_set <- parse_index_set(trm$retained_row_indices)
    union_set <- union(ref_set, trm_set)

    data.table(
      indicator = indicator,
      comparison = comparison_name,
      reference_scenario = reference_id,
      trimmed_scenario = trimmed_id,
      reference_covariates = ref$selected_covs,
      trimmed_covariates = trm$selected_covs,
      same_covariates = identical(
        sort(strsplit(ifelse(is.na(ref$selected_covs), "", ref$selected_covs),
                      ";", fixed = TRUE)[[1]]),
        sort(strsplit(ifelse(is.na(trm$selected_covs), "", trm$selected_covs),
                      ";", fixed = TRUE)[[1]])
      ),
      n_rows = trm$n_rows,
      n_retained_reference = length(ref_set),
      n_retained_trimmed = length(trm_set),
      n_trimmed = trm$n_rows - length(trm_set),
      trimmed_proportion = (trm$n_rows - length(trm_set)) / trm$n_rows,
      retained_set_jaccard = ifelse(
        length(union_set) == 0,
        NA_real_,
        length(intersect(ref_set, trm_set)) / length(union_set)
      ),
      delta_beta = trm$beta - ref$beta,
      beta_ratio_OR_scale = exp(trm$beta - ref$beta),
      delta_gamma = trm$gamma - ref$gamma,
      delta_gamma_ub = trm$gamma_ub - ref$gamma_ub,
      delta_BPRF = trm$BPRF_OR - ref$BPRF_OR,
      BPRF_ratio = trm$BPRF_OR / ref$BPRF_OR,
      delta_ROS = trm$ROS - ref$ROS,
      delta_stars = trm$stars - ref$stars
    )
  }

  trim_comparison_list <- Filter(
    Negate(is.null),
    list(
      make_trim_comparison(
        "S3_cvFE_trim1", "S4_cvFE_trim09",
        "adjusted_trim_effect_S4_minus_S3"
      )
    )
  )

  trim_effect_audit <- if (length(trim_comparison_list) > 0) {
    rbindlist(trim_comparison_list, fill = TRUE)
  } else {
    data.table(
      indicator = character(0),
      comparison = character(0),
      reference_scenario = character(0),
      trimmed_scenario = character(0)
    )
  }

  fwrite(
    trim_effect_audit,
    file.path(indicator_out_dir, paste0(indicator, "_TRIM_EFFECT_AUDIT.csv"))
  )
  all_summary[, `:=`(
    trim_reference_scenario = NA_character_,
    trim_effect_delta_beta = NA_real_,
    trim_effect_beta_ratio_OR_scale = NA_real_,
    trim_effect_delta_gamma = NA_real_,
    trim_effect_delta_gamma_ub = NA_real_,
    trim_effect_delta_BPRF = NA_real_,
    trim_effect_BPRF_ratio = NA_real_,
    trim_effect_delta_ROS = NA_real_,
    trim_effect_delta_stars = NA_real_,
    trim_effect_same_covariates = NA
  )]

  if (nrow(trim_effect_audit) > 0) {
    for (j in seq_len(nrow(trim_effect_audit))) {
      z <- trim_effect_audit[j]
      all_summary[scenario_id == z$trimmed_scenario, `:=`(
        trim_reference_scenario = z$reference_scenario,
        trim_effect_delta_beta = z$delta_beta,
        trim_effect_beta_ratio_OR_scale = z$beta_ratio_OR_scale,
        trim_effect_delta_gamma = z$delta_gamma,
        trim_effect_delta_gamma_ub = z$delta_gamma_ub,
        trim_effect_delta_BPRF = z$delta_BPRF,
        trim_effect_BPRF_ratio = z$BPRF_ratio,
        trim_effect_delta_ROS = z$delta_ROS,
        trim_effect_delta_stars = z$delta_stars,
        trim_effect_same_covariates = z$same_covariates
      )]
    }
  }

  # Match the repository dichotomous summary: report the all-observation
  # residual test from the final trimmed adjusted model (S4), while the S4
  # post result excludes observations with w_soln < 0.1.
  s4_egger_path <- file.path(
    indicator_out_dir, "S4_cvFE_trim09",
    paste0(indicator, "_S4_cvFE_trim09_IHME_residual_small_study_test.csv")
  )
  manuscript_egger <- if (file.exists(s4_egger_path)) {
    fread(s4_egger_path)[dataset == "pre_trimming"]
  } else data.table()

  manuscript_egger_status <- "ok_S4_all_observations_pre_exclusion"
  if (nrow(manuscript_egger) != 1L) {
    manuscript_egger_status <-
      "not_estimable_or_missing_S4_all_observation_result_indicator_results_preserved"
    manuscript_egger <- data.table(
      dataset = "pre_trimming",
      n_observations = NA_integer_,
      k_unique_studies = NA_integer_,
      mean_standardized_residual = NA_real_,
      se_mean = NA_real_,
      z_abs = NA_real_,
      p_one_sided = NA_real_,
      formally_interpretable_k_ge_10 = NA
    )
  }

  all_summary[, `:=`(
    manuscript_egger_source =
      "IHME_residual_test_S4_all_observations_fixed_effect_residual_sqrt_obs_se2_plus_gamma_one_sided",
    manuscript_egger_status = manuscript_egger_status,
    manuscript_egger_mean_standardized_residual =
      manuscript_egger$mean_standardized_residual[1],
    manuscript_egger_se_mean = manuscript_egger$se_mean[1],
    manuscript_egger_p = manuscript_egger$p_one_sided[1],
    manuscript_egger_k = manuscript_egger$k_unique_studies[1],
    manuscript_egger_formally_interpretable =
      manuscript_egger$formally_interpretable_k_ge_10[1]
  )]
  fwrite(
    manuscript_egger,
    file.path(indicator_out_dir, paste0(indicator, "_MANUSCRIPT_EGGER_S4_PRE.csv"))
  )

  summary_path <- file.path(
    indicator_out_dir,
    paste0(indicator, "_ADJUSTED_TWO_STAGE_COMPLETE.csv")
  )
  fwrite(all_summary, summary_path)

  interval_audit_cols <- intersect(
    c(
      "indicator", "scenario_id", "score_method", "score_model_type",
      "scenario_model_gamma", "score_model_gamma", "gamma_var", "gamma_sd",
      "gamma_ub", "score_refit_changed_gamma",
      "MRBRT_no_gamma_CI95_lo", "MRBRT_no_gamma_CI95_hi",
      "MRBRT_gamma_point_only_UI95_lo",
      "MRBRT_gamma_point_only_UI95_hi",
      "MRBRT_with_gamma_UI95_lo", "MRBRT_with_gamma_UI95_hi",
      "point_gamma_equal_no_gamma", "with_gamma_equal_no_gamma",
      "no_gamma_width_log", "gamma_point_only_width_log",
      "with_gamma_width_log", "with_gamma_interval_status",
      "gamma_fisher_boundary_used", "gamma_fisher_epsilon"
    ),
    names(all_summary)
  )
  fwrite(
    all_summary[, ..interval_audit_cols],
    file.path(indicator_out_dir, paste0(indicator, "_WITH_GAMMA_INTERVAL_AUDIT.csv"))
  )

  # Manuscript-focused adjusted pre/post trimming results
  main_table <- all_summary[
    scenario_id %in% c("S3_cvFE_trim1", "S4_cvFE_trim09")
  ]

  main_path <- file.path(
    indicator_out_dir,
    paste0(indicator, "_MANUSCRIPT_PRE_POST_TRIM.csv")
  )
  fwrite(main_table, main_path)

  # Score / plot provenance
  provenance_cols <- intersect(
    c(
      "indicator",
      "scenario_id",
      "score_method",
      "score_log",
      "BPRF_OR",
      "score_curve_log",
      "BPRF_curve_OR",
      "no_gamma_lo",
      "no_gamma_hi",
      "no_gamma_OR_lo",
      "no_gamma_OR_hi",
      "MRBRT_no_gamma_CI95_lo",
      "MRBRT_no_gamma_CI95_hi",
      "MRBRT_gamma_point_only_UI95_lo",
      "MRBRT_gamma_point_only_UI95_hi",
      "MRBRT_with_gamma_UI95_lo",
      "MRBRT_with_gamma_UI95_hi",
      "BPRF_conservative_UI95_lo",
      "BPRF_conservative_UI95_hi",
      "with_gamma_interval_status",
      "point_gamma_equal_no_gamma",
      "with_gamma_equal_no_gamma",
      "inner_lo",
      "inner_hi",
      "inner_OR_lo",
      "inner_OR_hi",
      "outer_lo",
      "outer_hi",
      "outer_OR_lo",
      "outer_OR_hi",
      "bprf_q05",
      "bprf_q95",
      "bprf_OR_q05",
      "bprf_OR_q95",
      "ROS",
      "ROS_rounded",
      "stars",
      "star_label",
      "score_direction",
      "gamma",
      "gamma_var",
      "gamma_ub",
      "score_model_type",
      "score_target_inliers",
      "score_requested_inliers",
      "trimming_count_matches_request",
      "score_model_n_rows",
      "score_status",
      "score_error"
    ),
    names(all_summary)
  )

  score_provenance <- all_summary[, ..provenance_cols]

  score_provenance_path <- file.path(
    indicator_out_dir,
    paste0(indicator, "_SCORE_PLOT_PROVENANCE.csv")
  )
  fwrite(score_provenance, score_provenance_path)

  message("\nCompleted indicator: ", indicator)

  all_summary
}


# ---------------------------------------------------------------------------
# Discover all indicator files
# ---------------------------------------------------------------------------

indicator_files <- list.files(
  in_dir,
  pattern = "_all_extraction\\.csv$",
  full.names = TRUE
)

if (length(indicator_files) == 0) {
  stop(
    "No input files matching *_all_extraction.csv were found in: ",
    in_dir
  )
}

indicator_files <- sort(indicator_files)

message("\n============================================================")
message("BATCH BPRF ANALYSIS")
message("Input directory: ", in_dir)
message("Number of indicators: ", length(indicator_files))
message("============================================================")


# ---------------------------------------------------------------------------
# Batch loop
# ---------------------------------------------------------------------------

batch_results <- vector("list", length(indicator_files))
batch_status <- vector("list", length(indicator_files))

for (i in seq_along(indicator_files)) {

  f <- indicator_files[i]
  indicator_name <- sub("_all_extraction\\.csv$", "", basename(f))

  message("\n\n############################################################")
  message("[", i, "/", length(indicator_files), "] ", indicator_name)
  message("############################################################")

  result <- tryCatch(
    {
      res <- run_one_indicator_file(f)

      batch_status[[i]] <- data.table(
        indicator = indicator_name,
        input_file = f,
        batch_status = "completed",
        batch_error = NA_character_
      )

      res
    },
    error = function(e) {

      batch_status[[i]] <- data.table(
        indicator = indicator_name,
        input_file = f,
        batch_status = "failed",
        batch_error = conditionMessage(e)
      )

      warning(
        "Indicator failed: ",
        indicator_name,
        " | ",
        conditionMessage(e),
        call. = FALSE
      )

      if (!CONTINUE_ON_ERROR) {
        stop(e)
      }

      data.table(
        indicator = indicator_name,
        scenario_id = NA_character_,
        score_status = "indicator_failed",
        score_error = conditionMessage(e)
      )
    }
  )

  batch_results[[i]] <- result

  # Progressive checkpoint after every indicator
  status_so_far <- rbindlist(
    batch_status[seq_len(i)],
    fill = TRUE
  )

  fwrite(
    status_so_far,
    file.path(
      save_dir,
      "CHECKPOINT_BPRF_binary_batch_status.csv"
    )
  )

  results_so_far <- rbindlist(
    batch_results[seq_len(i)],
    fill = TRUE
  )

  fwrite(
    results_so_far,
    file.path(
      save_dir,
      "CHECKPOINT_BPRF_binary_results_so_far.csv"
    )
  )
}


# ---------------------------------------------------------------------------
# Final combined outputs
# ---------------------------------------------------------------------------

all_results <- rbindlist(batch_results, fill = TRUE)
status_dt <- rbindlist(batch_status, fill = TRUE)

fwrite(
  status_dt,
  file.path(
    save_dir,
    "00_BPRF_binary_BATCH_STATUS.csv"
  )
)

fwrite(
  all_results,
  file.path(
    save_dir,
    "01_BPRF_binary_ALL_INDICATORS_ALL_SCENARIOS.csv"
  )
)

# Successful BPRF rows only
successful_results <- all_results[
  !is.na(scenario_id) &
    (
      is.na(score_status) |
      grepl("^ok", score_status)
    )
]

fwrite(
  successful_results,
  file.path(
    save_dir,
    "02_BPRF_binary_SUCCESSFUL_SCENARIOS.csv"
  )
)

# S3-S4 wide summary for quick QC
if (nrow(successful_results) > 0) {

  summary_cols <- intersect(
    c(
      "indicator",
      "scenario_id",
      "n_rows",
      "n_included_weight_gt0",
      "selected_covs",
      "selected_cov_count",
      "covariate_mode",
      "covariate_lock_preserved",
      "score_covariate_lock_preserved",
      "score_model_type",
      "score_target_inliers",
      "score_model_n_rows",
      "beta",
      "no_gamma_OR_lo",
      "no_gamma_OR_hi",
      "MRBRT_no_gamma_CI95_lo",
      "MRBRT_no_gamma_CI95_hi",
      "MRBRT_gamma_point_only_UI95_lo",
      "MRBRT_gamma_point_only_UI95_hi",
      "MRBRT_with_gamma_UI95_lo",
      "MRBRT_with_gamma_UI95_hi",
      "BPRF_conservative_UI95_lo",
      "BPRF_conservative_UI95_hi",
      "point_gamma_equal_no_gamma",
      "with_gamma_equal_no_gamma",
      "with_gamma_interval_status",
      "scenario_model_gamma",
      "score_model_gamma",
      "score_refit_changed_gamma",
      "inner_OR_lo",
      "inner_OR_hi",
      "outer_OR_lo",
      "outer_OR_hi",
      "bprf_OR_q05",
      "bprf_OR_q95",
      "BPRF_OR",
      "BPRF_curve_OR",
      "ROS",
      "ROS_rounded",
      "stars",
      "star_label",
      "gamma",
      "gamma_var",
      "gamma_ub",
      "egger_pre_intercept",
      "egger_pre_p",
      "egger_post_intercept",
      "egger_post_p",
      "manuscript_egger_source",
      "manuscript_egger_status",
      "manuscript_egger_mean_standardized_residual",
      "manuscript_egger_se_mean",
      "manuscript_egger_p",
      "manuscript_egger_k",
      "manuscript_egger_formally_interpretable",
      "trim_reference_scenario",
      "trim_effect_delta_beta",
      "trim_effect_beta_ratio_OR_scale",
      "trim_effect_delta_gamma",
      "trim_effect_delta_gamma_ub",
      "trim_effect_delta_BPRF",
      "trim_effect_BPRF_ratio",
      "trim_effect_delta_ROS",
      "trim_effect_delta_stars",
      "trim_effect_same_covariates"
    ),
    names(successful_results)
  )

  compact_summary <- successful_results[, ..summary_cols]

  fwrite(
    compact_summary,
    file.path(
      save_dir,
      "03_BPRF_binary_COMPACT_SUMMARY.csv"
    )
  )

  # Primary manuscript table:
  # adjusted no-trim + adjusted post-trim.
  manuscript_main <- successful_results[
    scenario_id %in% c("S3_cvFE_trim1", "S4_cvFE_trim09")
  ]

  fwrite(
    manuscript_main,
    file.path(
      save_dir,
      "04_MANUSCRIPT_MAIN_BPRF_PRE_POST_TRIM.csv"
    )
  )

  # Primary post-trimming S4 only
  primary_s4 <- successful_results[
    scenario_id == "S4_cvFE_trim09"
  ]

  fwrite(
    primary_s4,
    file.path(
      save_dir,
      "05_PRIMARY_S4_ADJUSTED_TRIM09.csv"
    )
  )

  # Score/plot provenance across all indicators
  provenance_cols <- intersect(
    c(
      "indicator",
      "scenario_id",
      "score_method",
      "score_log",
      "BPRF_OR",
      "score_curve_log",
      "BPRF_curve_OR",
      "no_gamma_lo",
      "no_gamma_hi",
      "no_gamma_OR_lo",
      "no_gamma_OR_hi",
      "MRBRT_no_gamma_CI95_lo",
      "MRBRT_no_gamma_CI95_hi",
      "MRBRT_gamma_point_only_UI95_lo",
      "MRBRT_gamma_point_only_UI95_hi",
      "MRBRT_with_gamma_UI95_lo",
      "MRBRT_with_gamma_UI95_hi",
      "BPRF_conservative_UI95_lo",
      "BPRF_conservative_UI95_hi",
      "point_gamma_equal_no_gamma",
      "with_gamma_equal_no_gamma",
      "with_gamma_interval_status",
      "scenario_model_gamma",
      "score_model_gamma",
      "score_refit_changed_gamma",
      "inner_lo",
      "inner_hi",
      "inner_OR_lo",
      "inner_OR_hi",
      "outer_lo",
      "outer_hi",
      "outer_OR_lo",
      "outer_OR_hi",
      "bprf_q05",
      "bprf_q95",
      "bprf_OR_q05",
      "bprf_OR_q95",
      "ROS",
      "ROS_rounded",
      "stars",
      "star_label",
      "score_direction",
      "gamma",
      "gamma_var",
      "gamma_ub",
      "score_model_type",
      "score_target_inliers",
      "score_model_n_rows"
    ),
    names(successful_results)
  )

  fwrite(
    successful_results[, ..provenance_cols],
    file.path(
      save_dir,
      "06_SCORE_PLOT_PROVENANCE_ALL_INDICATORS.csv"
    )
  )

  combined_interval_audit_cols <- intersect(
    c(
      "indicator", "scenario_id", "score_method", "score_model_type",
      "scenario_model_gamma", "score_model_gamma", "gamma_var", "gamma_sd",
      "gamma_ub", "score_refit_changed_gamma",
      "MRBRT_no_gamma_CI95_lo", "MRBRT_no_gamma_CI95_hi",
      "MRBRT_gamma_point_only_UI95_lo",
      "MRBRT_gamma_point_only_UI95_hi",
      "MRBRT_with_gamma_UI95_lo", "MRBRT_with_gamma_UI95_hi",
      "point_gamma_equal_no_gamma", "with_gamma_equal_no_gamma",
      "no_gamma_width_log", "gamma_point_only_width_log",
      "with_gamma_width_log", "with_gamma_interval_status",
      "gamma_fisher_boundary_used", "gamma_fisher_epsilon"
    ),
    names(successful_results)
  )
  fwrite(
    successful_results[, ..combined_interval_audit_cols],
    file.path(save_dir, "08_WITH_GAMMA_INTERVAL_AUDIT_ALL.csv")
  )
}

# Failure summary
failure_rows <- all_results[
  grepl(
    "failed|not_estimable",
    ifelse(is.na(score_status), "", score_status),
    ignore.case = TRUE
  ) |
    (!is.na(score_error) & nzchar(score_error))
]

fwrite(
  failure_rows,
  file.path(
    save_dir,
    "07_BPRF_FAILURE_LOG.csv"
  )
)

message("\n============================================================")
message("BATCH ANALYSIS COMPLETE")
message("Output directory: ", save_dir)
message("Indicators completed: ", sum(status_dt$batch_status == "completed"))
message("Indicators failed:    ", sum(status_dt$batch_status == "failed"))
message("============================================================")

print(status_dt)

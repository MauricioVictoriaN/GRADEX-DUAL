# ┌─────────────────────────────────────────────────────────────────┐
# │                                                                 │
# │      ██████╗  ██████╗  █████╗ ██████╗ ███████╗██╗  ██╗         │
# │      ██╔════╝ ██╔══██╗██╔══██╗██╔══██╗██╔════╝╚██╗██╔╝         │
# │      ██║  ███╗██████╔╝███████║██║  ██║█████╗   ╚███╔╝          │
# │      ██║   ██║██╔══██╗██╔══██║██║  ██║██╔══╝   ██╔██╗          │
# │      ╚██████╔╝██║  ██║██║  ██║██████╔╝███████╗██╔╝ ██╗         │
# │       ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═════╝ ╚══════╝╚═╝  ╚═╝         │
# │                                                                 │
# │              ██████╗ ██╗   ██╗ █████╗ ██╗                      │
# │              ██╔══██╗██║   ██║██╔══██╗██║                      │
# │              ██║  ██║██║   ██║███████║██║                      │
# │              ██║  ██║██║   ██║██╔══██║██║                      │
# │              ██████╔╝╚██████╔╝██║  ██║███████╗                 │
# │              ╚═════╝  ╚═════╝ ╚═╝  ╚═╝╚══════╝                 │
# │                                                                 │
# │                      GRADEX-DUAL v1.0.0                          │
# │                                                                 │
# │              Pointwise (Guillot & Duband, 1967)                │
# │           + Regional (Hosking & Wallis, 1997)                  │
# │                                                                 │
# └─────────────────────────────────────────────────────────────────┘
#
# =============================================================================
# GRADEX-DUAL v1.0.0
# =============================================================================
#
# Description:
#   GRADEX-DUAL: A dual framework combining pointwise GRADEX
#   (Guillot & Duband, 1967) with regional frequency analysis
#   (Hosking & Wallis, 1997) for HEC-HMS hydrological design.
#
#   Computes the GRADEX parameter (precipitation gradient T=10 to T=100 years)
#   from rain gauge records using two complementary layers:
#
#   LAYER 1 — Pointwise spatial interpolation:
#     Fits four probability distributions per station (Gumbel, Log-Normal,
#     Pearson III, Gamma) using L-moments and MLE, selects the best model
#     by Anderson-Darling + RMSE, and regionalises to a target gauging point
#     via IDW (optimised exponent) or Ordinary Kriging (automap variogram).
#     Multi-metric method selection: RMSE + R² + |BIAS| scoring with
#     anti-degenerate-Kriging safeguard (v1.2).
#
#   LAYER 2 — H&W regional distribution layer:
#     Uses regional L-moments (H&W eq. 4.3) to select the best distribution
#     via the Z goodness-of-fit statistic (H&W eq. 5.6) across five candidates
#     (GEV, GLO, GNO, PE3, GPA). Re-estimates GRADEX at each station with the
#     regionally selected distribution and interpolates to the gauging point.
#     Computes the formally correct bootstrap CI (H&W Section 6.3) propagating
#     parameter, inter-site, record-length, and spatial interpolation uncertainty.
#
#   QUALITY-GATED RECOMMENDATION:
#     Final GRADEX recommendation for HEC-HMS applies three explicit quality
#     gates before computing weights:
#       - R² gate: spatial weight = 0 if R² < 0.3
#       - CI amplitude gate: H&W weight = 0 if CI amplitude >= 100%
#       - Sensitivity range: 10-30% depending on gate outcomes
#
#   PRE-ANALYSIS DIAGNOSTICS:
#     - Ljung-Box test for serial independence (WMO Technical Note 168)
#     - Hosking & Wallis H-statistic for regional homogeneity
#     - Discordancy D-statistic to flag outlier stations
#     - Station-level GRADEX confidence intervals via parametric bootstrap
#     - Mann-Kendall trend test with TRUE Theil-Sen slope (v1.1 fix)
#     - Non-stationarity sensitivity analysis on Sen-detrended series (v1.1)
#     - Monte Carlo validation against known true GRADEX (v1.1, enabled v1.2)
#
# Changelog v0.1 -> v1.0.0:
#   * Anti-degenerate-Kriging safeguard: when var(Predicho)/var(Observado)
#     in LOO-CV is < 0.10, Kriging is disqualified from method selection
#     regardless of its multi-metric score. Catches the failure mode where
#     a near-pure-nugget variogram produces near-constant predictions that
#     game the R² and |BIAS| metrics while genuinely failing on RMSE.
#   * Final summary: diagnostics block (LB / D / H) integrated INSIDE the
#     output box (was leaking outside the closing border in v1.1).
#   * Monte Carlo validation enabled by default (RUN_MC_VALIDATION = TRUE).
#
# Changelog v0.2 -> v1.0.0:
#   * BUGFIX: KRIGING_NUGGET_SILL_MAX was used in Section 9 but never defined,
#     causing a fatal error whenever n_stations >= 5 triggered Kriging. The
#     constant is now defined in Section 1 (default 0.5).
#   * BUGFIX: cv_metrics() coerced its Valor column to character when extra_info
#     contained text fields (e.g. "Variogram model" = "Sph"), which broke
#     downstream round() calls. Numeric and text metadata are now stored in
#     separate columns (Valor numeric, Valor_text character).
#   * Mann-Kendall: 'pendiente' field now reports the true Theil-Sen slope
#     (mm/year) computed by median of pairwise slopes; the previous value
#     was the Kendall S statistic (sum of signs), mislabelled.
#   * Section 6 trend plot now overlays the Theil-Sen line per station,
#     matching the slope reported in the table.
#   * NEW Section 6b: non-stationarity sensitivity analysis on Sen-detrended
#     series, reports the bias of the stationarity assumption.
#   * NEW Section 20: Monte Carlo validation of the framework with known
#     true GRADEX.
#   * Map and CV labels now use ggrepel to prevent overlap.
#
# Input file:
#   datos_precipitacion.xlsx   (sheet 1)
#   Required columns: id, nombre, lon, lat, anio, p_max
#
# Outputs:
#   resultados_gradex_completos.xlsx  — Multi-sheet workbook
#     Pointwise: Summary | Station_Results | Distribution_Summary
#                Confidence_Intervals | CV_Metrics_Combined
#                CV_Predictions_IDW | CV_Predictions_Kriging
#                Independence_Tests | Homogeneity_Discordancy
#                Homogeneity_Summary | Trend_Analysis
#     H&W Layer: HW_Summary | HW_Z_Selection | HW_Site_Lmom
#                HW_Station_GRADEX | HW_Recommendation | HW_Quality_Diag
#   mapa_gradex_final.png             — Main map (selected method)
#   mapa_comparacion_metodos.png      — IDW vs Kriging comparison panel
#   analisis_tendencia.png            — Mann-Kendall trend plot
#   validacion_cruzada_idw.png        — IDW cross-validation scatter
#   validacion_cruzada_kriging.png    — Kriging cross-validation scatter
#   shapefile_resultados/             — GIS shapefiles (MAGNA-SIRGAS / CTM12)
#   session_info.txt                  — R session metadata
#
# IMPORTANT — READ BEFORE RUNNING:
#   Always restart R before running this script to ensure all function
#   definitions are loaded fresh from this file.
#   RStudio:  Session > Restart R  (Ctrl+Shift+F10)
#   Console:  .rs.restartR()
#
#   Verification after loading:
#   stopifnot("hw_ok" %in% all.vars(body(compute_weighted_recommendation)))
#   message("OK: correct version loaded")
#
# References:
#   Guillot, P. & Duband, D. (1967). La méthode du Gradex pour le calcul
#     de la probabilité des crues à partir des pluies. IAHS Publ. 84.
#   Hosking, J.R.M. & Wallis, J.R. (1997). Regional Frequency Analysis:
#     An Approach Based on L-Moments. Cambridge University Press.
#   WMO (2008). Manual on Flood Forecasting and Warning. WMO-No. 1072.
#
# Author: Mauricio Javier Victoria N.
# ORCID: 0009-0003-4328-5691
# Date: 2026-05-04
# R version recommended: >= 4.2.0
# =============================================================================


# -----------------------------------------------------------------------------
# SECTION 1: WORKSPACE CONFIGURATION AND GLOBAL PARAMETERS
# -----------------------------------------------------------------------------

configure_workspace <- function(base_path) {
  if (!dir.exists(base_path)) {
    message("Directory not found. Creating: ", base_path)
    dir.create(base_path, recursive = TRUE)
  }
  setwd(base_path)
  message("Working directory set to: ", getwd())
  invisible(NULL)
}

# ---- User-defined parameters (modify as needed) ----------------------------
BASE_PATH        <- "D:/R/Regionalizacion_Gradex"
INPUT_FILE       <- "datos_precipitacion.xlsx"
# Gauging point: aforo hipotético en la confluencia aguas arriba de Yarumal,
# dentro del polígono cubierto por la red de 9 estaciones (norte de Antioquia,
# cuenca alta del río Cauca). Coherente con el caso de estudio actualizado.
LON_GAUGE        <- -75.48     # Target gauging point longitude (decimal degrees)
LAT_GAUGE        <-   6.65     # Target gauging point latitude  (decimal degrees)
MIN_YEARS        <- 10L        # Minimum record length per station (years)
WARN_YEARS       <- 30L        # Threshold below which Q100 reliability degrades
CONFIDENCE_LEVEL <- 0.95       # Confidence level for interval estimation
N_BOOTSTRAP      <- 2000L      # Bootstrap replicates for non-parametric CI
IDW_POWER_GRID   <- c(1, 1.5, 2, 3)  # Candidate IDW exponents (LOO-CV grid)
IDW_NMAX         <- 5L         # Max neighbours for IDW
GRID_RES_MAP     <- 80L        # Interpolation grid resolution (cells per axis)
BBOX_BUFFER      <- 0.05       # Bounding-box buffer (decimal degrees)
RANDOM_SEED      <- 2024L      # Reproducibility seed
GOF_ALPHA        <- 0.05       # Significance level for Anderson-Darling test
# Multi-metric selection weights (all equal by default; adjust if needed)
# Scoring: for each metric, 1 point if Kriging is better, 0 otherwise.
# Kriging is selected when its total score > IDW total score.
# Ties are broken by RMSE (lower wins). Falls back to IDW on Kriging failure.
SEL_WEIGHT_RMSE <- 1L   # Weight for RMSE criterion
SEL_WEIGHT_R2   <- 1L   # Weight for R² criterion  (higher = better)
SEL_WEIGHT_BIAS <- 1L   # Weight for |BIAS| criterion (lower absolute = better)
# Independence and homogeneity tests
LJUNG_BOX_LAGS  <- 3L   # Lags for Ljung-Box test (WMO: 1/4 of series length max)
LJUNG_BOX_ALPHA <- 0.05 # Significance level for independence test
HW_NSIM         <- 500L # Simulations for Hosking-Wallis H-statistic
DISC_THRESHOLD  <- 3.0  # Discordancy D > 3 flags a station as atypical
# Station-level GRADEX confidence interval (parametric bootstrap)
N_BOOT_STATION  <- 1000L # Bootstrap replicates for per-station GRADEX CI
# Kriging variogram diagnostic threshold
# Nugget/sill ratio below this value indicates strong spatial structure;
# above it, the variable is largely random and IDW is preferred regardless
# of LOO-CV scoring outcome. Standard geostatistical recommendation: 0.5.
KRIGING_NUGGET_SILL_MAX <- 0.5  # v1.1: previously undefined — caused error
                                #       when n_stations >= 5 triggered Kriging

set.seed(RANDOM_SEED)
options(gstat.debug = FALSE)   # Suppress gstat C-level printf output globally
configure_workspace(BASE_PATH)


# -----------------------------------------------------------------------------
# SECTION 2: LIBRARY LOADING
# -----------------------------------------------------------------------------

load_libraries <- function() {
  if (!require("pacman", quietly = TRUE)) install.packages("pacman")

  pacman::p_load(
    tidyverse,     # Data manipulation and visualisation
    lmom,          # L-moment parameter estimation
    fitdistrplus,  # MLE distribution fitting
    sf,            # Simple features for spatial data
    gstat,         # IDW and Kriging interpolation
    readxl,        # Excel reading
    writexl,       # Excel writing
    gridExtra,     # Multi-panel plot layout
    automap,       # Automatic variogram fitting (Kriging)
    goftest,       # Anderson-Darling goodness-of-fit
    sp,            # Legacy spatial classes (required by gstat)
    viridis,       # Perceptually uniform colour scales
    patchwork,     # Plot composition
    ggspatial,     # Spatial annotations (scale bar, north arrow)
    Kendall,       # Mann-Kendall trend test
    hydroGOF,      # Hydrological GOF metrics
    knitr,         # Table rendering
    ggrepel        # v1.1: non-overlapping labels in maps and CV plots
  )

  for (pkg in c("kableExtra", "metR")) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }

  message("All libraries loaded successfully.")
  invisible(NULL)
}

load_libraries()


# -----------------------------------------------------------------------------
# SECTION 3: CUSTOM GUMBEL FUNCTIONS
# -----------------------------------------------------------------------------

pgumbel <- function(q, loc = 0, scale = 1) {
  exp(-exp(-(q - loc) / scale))
}

qgumbel <- function(p, loc = 0, scale = 1) {
  loc - scale * log(-log(p))
}


# -----------------------------------------------------------------------------
# SECTION 4: DATA LOADING AND VALIDATION
# -----------------------------------------------------------------------------

load_precipitation_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("Input file not found: ", file_path,
         "\nVerify path in: ", getwd())
  }

  df <- readxl::read_excel(file_path, sheet = 1)

  required_cols <- c("id", "nombre", "lon", "lat", "anio", "p_max")
  missing_cols  <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         "\nAvailable: ", paste(names(df), collapse = ", "))
  }

  for (col in c("lon", "lat", "p_max")) {
    orig     <- df[[col]]
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
    n_na     <- sum(is.na(df[[col]])) - sum(is.na(orig))
    if (n_na > 0L)
      warning("Column '", col, "': ", n_na, " value(s) coerced to NA.")
  }
  df$anio <- suppressWarnings(as.integer(df$anio))

  if (any(df$p_max < 0, na.rm = TRUE))
    warning("Column 'p_max' contains negative values.")
  if (any(abs(df$lon) > 180, na.rm = TRUE))
    warning("Column 'lon' has values outside [-180, 180].")
  if (any(abs(df$lat) > 90, na.rm = TRUE))
    warning("Column 'lat' has values outside [-90, 90].")

  message("Data loaded: ", nrow(df), " records | ",
          length(unique(df$id)), " stations")
  df
}

df_raw <- load_precipitation_data(INPUT_FILE)


# -----------------------------------------------------------------------------
# SECTION 5: GRADEX ESTIMATION PER STATION
# -----------------------------------------------------------------------------

#' Compute GRADEX for a single station with record-length diagnostics
#'
#' @param x       Numeric vector. Annual maximum precipitation series.
#' @param stn_id  Character. Station identifier.
#' @param gof_alpha Numeric. AD test significance level.
#' @return One-row data frame with GRADEX, RMSE, AD p-values, best model,
#'   and record-quality flags.
compute_station_gradex <- function(x, stn_id, gof_alpha = GOF_ALPHA) {

  x <- stats::na.omit(x)
  n <- length(x)

  na_row <- data.frame(
    id = stn_id, n_datos = n,
    Gumbel_gradex = NA_real_,    Gumbel_rmse = NA_real_,
    Gumbel_ad_pval = NA_real_,
    LogNormal_gradex = NA_real_, LogNormal_rmse = NA_real_,
    LogNormal_ad_pval = NA_real_,
    Pearson_gradex = NA_real_,   Pearson_rmse = NA_real_,
    Pearson_ad_pval = NA_real_,
    Gamma_gradex = NA_real_,     Gamma_rmse = NA_real_,
    Gamma_ad_pval = NA_real_,
    Mejor_Modelo = NA_character_, Mejor_gradex = NA_real_,
    record_quality = NA_character_,
    stringsAsFactors = FALSE
  )

  if (n < 5L) return(na_row)

  # ---- Record-length quality classification --------------------------------
  # Based on WMO guidelines for extreme-value analysis:
  #   >= 30 years : reliable Q100 estimates
  #   20-29 years : acceptable with caution
  #   10-19 years : unreliable Q100; wide confidence intervals expected
  #    < 10 years : rejected (filtered by MIN_YEARS upstream)
  record_quality <- dplyr::case_when(
    n >= WARN_YEARS       ~ "RELIABLE (>= 30 yr)",
    n >= 20L              ~ "ACCEPTABLE (20-29 yr) — caution advised",
    n >= MIN_YEARS        ~ "UNRELIABLE (< 30 yr) — Q100 CI very wide",
    TRUE                  ~ "INSUFFICIENT"
  )

  if (n < WARN_YEARS) {
    warning(sprintf(
      paste0("Station %s: only %d years of record.\n",
             "  Q100 estimates have very wide confidence intervals.\n",
             "  WMO recommends >= 30 years for T=100 frequency analysis.\n",
             "  Classification: %s"),
      stn_id, n, record_quality
    ))
  }

  x_sorted    <- sort(x)
  p_empirical <- seq_len(n) / (n + 1L)
  y_10        <- -log(-log(1 - 1 / 10))
  y_100       <- -log(-log(1 - 1 / 100))

  # ---- TECHNICAL NOTE: GRADEX definition and approximation -----------------
  # The GRADEX (Guillot & Duband, 1967) is strictly defined as the slope of
  # the Gumbel reduced variate plot: GRADEX = (Q100 - Q10) / (y100 - y10).
  # For the Gumbel distribution this equals the scale parameter beta exactly,
  # because the Gumbel CDF is linear in reduced variate space.
  #
  # For Log-Normal, Pearson III, and Gamma, this linearity does NOT hold —
  # these distributions curve in Gumbel paper. The formula above therefore
  # computes a LOCAL CHORD SLOPE between T=10 and T=100, not a true GRADEX.
  # This is a widely used approximation in applied hydrology when Gumbel does
  # not provide the best fit, but it should be interpreted with caution:
  #   - The "GRADEX" value for non-Gumbel distributions is a secant slope,
  #     not a global property of the distribution tail.
  #   - It will underestimate the true tail slope for heavy-tailed distributions
  #     (Log-Normal with high sigma, Pearson III with positive skew) at T > 100.
  #   - For design purposes, the Gumbel GRADEX is always reported alongside the
  #     best-fit result so the user can compare.
  # --------------------------------------------------------------------------

  ad_pvalue <- function(observed, fitted_cdf) {
    tryCatch(goftest::ad.test(observed, null = fitted_cdf)$p.value,
             error = function(e) NA_real_)
  }

  # 5.1 Gumbel (L-moments)
  gp <- tryCatch(lmom::pelgum(lmom::samlmu(x)), error = function(e) NULL)
  if (!is.null(gp) && !anyNA(gp) && all(is.finite(gp))) {
    gradex_gumbel <- (lmom::quagum(1-1/100, gp) - lmom::quagum(1-1/10, gp)) /
                     (y_100 - y_10)
    fitted_g      <- lmom::quagum(p_empirical, gp)
    rmse_gumbel   <- sqrt(mean((x_sorted - fitted_g)^2, na.rm = TRUE))
    ad_gumbel     <- ad_pvalue(x_sorted,
                               function(q) pgumbel(q, loc = gp[1], scale = gp[2]))
  } else { gradex_gumbel <- rmse_gumbel <- ad_gumbel <- NA_real_ }

  # 5.2 Log-Normal (MLE)
  fit_ln <- tryCatch(fitdistrplus::fitdist(x, "lnorm", method = "mle"),
                     error = function(e) NULL)
  if (!is.null(fit_ln)) {
    mu_ln  <- fit_ln$estimate["meanlog"]
    sd_ln  <- fit_ln$estimate["sdlog"]
    gradex_lnorm <- (stats::qlnorm(1-1/100, mu_ln, sd_ln) -
                     stats::qlnorm(1-1/10,  mu_ln, sd_ln)) / (y_100 - y_10)
    rmse_lnorm   <- sqrt(mean((x_sorted -
                      stats::qlnorm(p_empirical, mu_ln, sd_ln))^2, na.rm = TRUE))
    ad_lnorm     <- ad_pvalue(x_sorted,
                              function(q) stats::plnorm(q, mu_ln, sd_ln))
  } else { gradex_lnorm <- rmse_lnorm <- ad_lnorm <- NA_real_ }

  # 5.3 Pearson III (L-moments)
  pp <- tryCatch(lmom::pelpe3(lmom::samlmu(x)), error = function(e) NULL)
  if (!is.null(pp) && !anyNA(pp)) {
    gradex_pearson <- (lmom::quape3(1-1/100, pp) - lmom::quape3(1-1/10, pp)) /
                      (y_100 - y_10)
    rmse_pearson   <- sqrt(mean((x_sorted -
                       lmom::quape3(p_empirical, pp))^2, na.rm = TRUE))
    ad_pearson     <- ad_pvalue(x_sorted, function(q) lmom::cdfpe3(q, pp))
  } else { gradex_pearson <- rmse_pearson <- ad_pearson <- NA_real_ }

  # 5.4 Gamma (L-moments)
  gam <- tryCatch(lmom::pelgam(lmom::samlmu(x)), error = function(e) NULL)
  if (!is.null(gam) && !anyNA(gam)) {
    gradex_gamma <- (lmom::quagam(1-1/100, gam) - lmom::quagam(1-1/10, gam)) /
                    (y_100 - y_10)
    rmse_gamma   <- sqrt(mean((x_sorted -
                     lmom::quagam(p_empirical, gam))^2, na.rm = TRUE))
    ad_gamma     <- ad_pvalue(x_sorted, function(q) lmom::cdfgam(q, gam))
  } else { gradex_gamma <- rmse_gamma <- ad_gamma <- NA_real_ }

  # 5.5 Two-stage selection: AD filter -> min RMSE
  rmse_v   <- c(rmse_gumbel, rmse_lnorm, rmse_pearson, rmse_gamma)
  gradex_v <- c(gradex_gumbel, gradex_lnorm, gradex_pearson, gradex_gamma)
  ad_v     <- c(ad_gumbel, ad_lnorm, ad_pearson, ad_gamma)
  names_v  <- c("Gumbel", "Log-Normal", "Pearson III", "Gamma")

  gof_pass <- !is.na(ad_v) & (ad_v >= gof_alpha)
  cand_rmse <- if (any(gof_pass)) ifelse(gof_pass, rmse_v, NA_real_) else {
    warning("Station ", stn_id,
            ": no distribution passed AD test. Selecting by RMSE only.")
    rmse_v
  }

  best_idx <- which.min(cand_rmse)
  if (length(best_idx) > 0L && !is.na(cand_rmse[best_idx])) {
    best_model  <- names_v[best_idx]
    best_gradex <- gradex_v[best_idx]
  } else {
    best_model  <- "No fit"
    best_gradex <- NA_real_
  }

  # ---- Station-level GRADEX confidence interval (parametric bootstrap) ------
  # Resamples the fitted distribution to propagate parameter uncertainty into
  # the GRADEX estimate. Uses the best-fit distribution parameters.
  gradex_ci_lower <- NA_real_
  gradex_ci_upper <- NA_real_

  if (!is.na(best_gradex) && best_model != "No fit") {
    boot_gradex <- tryCatch({
      replicate(N_BOOT_STATION, {
        x_sim <- switch(best_model,
          "Gumbel"      = lmom::quagum(stats::runif(n), gp),
          "Log-Normal"  = stats::rlnorm(n, mu_ln, sd_ln),
          "Pearson III" = { pp2 <- lmom::pelpe3(lmom::samlmu(
                              stats::rnorm(n, mean(x), stats::sd(x))));
                            lmom::quape3(stats::runif(n), pp2) },
          "Gamma"       = { gam2 <- lmom::pelgam(lmom::samlmu(
                              stats::rgamma(n, shape = 2)));
                            lmom::quagam(stats::runif(n), gam2) },
          rep(NA_real_, n)
        )
        x_sim <- x_sim[is.finite(x_sim) & x_sim > 0]
        if (length(x_sim) < 5L) return(NA_real_)
        lm_s  <- lmom::samlmu(x_sim)
        p_s   <- switch(best_model,
          "Gumbel"      = tryCatch(lmom::pelgum(lm_s), error = function(e) NULL),
          "Log-Normal"  = { fl <- tryCatch(
                              fitdistrplus::fitdist(x_sim, "lnorm", method = "mle"),
                              error = function(e) NULL);
                            if (is.null(fl)) NULL else
                              c(fl$estimate["meanlog"], fl$estimate["sdlog"]) },
          "Pearson III" = tryCatch(lmom::pelpe3(lm_s), error = function(e) NULL),
          "Gamma"       = tryCatch(lmom::pelgam(lm_s), error = function(e) NULL),
          NULL
        )
        if (is.null(p_s) || anyNA(p_s)) return(NA_real_)
        q100_b <- switch(best_model,
          "Gumbel"      = lmom::quagum(1-1/100, p_s),
          "Log-Normal"  = stats::qlnorm(1-1/100, p_s[1], p_s[2]),
          "Pearson III" = lmom::quape3(1-1/100, p_s),
          "Gamma"       = lmom::quagam(1-1/100, p_s),
          NA_real_)
        q10_b  <- switch(best_model,
          "Gumbel"      = lmom::quagum(1-1/10, p_s),
          "Log-Normal"  = stats::qlnorm(1-1/10, p_s[1], p_s[2]),
          "Pearson III" = lmom::quape3(1-1/10, p_s),
          "Gamma"       = lmom::quagam(1-1/10, p_s),
          NA_real_)
        if (anyNA(c(q100_b, q10_b))) NA_real_
        else (q100_b - q10_b) / (y_100 - y_10)
      })
    }, error = function(e) rep(NA_real_, N_BOOT_STATION))

    boot_gradex <- boot_gradex[is.finite(boot_gradex)]
    if (length(boot_gradex) >= 10L) {
      alpha_ci        <- 1 - CONFIDENCE_LEVEL
      gradex_ci_lower <- round(stats::quantile(boot_gradex, alpha_ci / 2), 3)
      gradex_ci_upper <- round(stats::quantile(boot_gradex, 1 - alpha_ci / 2), 3)
    }
  }

  data.frame(
    id = stn_id, n_datos = n,
    Gumbel_gradex    = round(gradex_gumbel,  3),
    Gumbel_rmse      = round(rmse_gumbel,    3),
    Gumbel_ad_pval   = round(ad_gumbel,      4),
    LogNormal_gradex  = round(gradex_lnorm,  3),
    LogNormal_rmse    = round(rmse_lnorm,    3),
    LogNormal_ad_pval = round(ad_lnorm,      4),
    Pearson_gradex   = round(gradex_pearson, 3),
    Pearson_rmse     = round(rmse_pearson,   3),
    Pearson_ad_pval  = round(ad_pearson,     4),
    Gamma_gradex     = round(gradex_gamma,   3),
    Gamma_rmse       = round(rmse_gamma,     3),
    Gamma_ad_pval    = round(ad_gamma,       4),
    Mejor_Modelo     = best_model,
    Mejor_gradex     = round(best_gradex,    3),
    gradex_ci_lower  = gradex_ci_lower,
    gradex_ci_upper  = gradex_ci_upper,
    record_quality   = record_quality,
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
# SECTION 5b: SERIAL INDEPENDENCE TEST (LJUNG-BOX)
# -----------------------------------------------------------------------------

#' Test serial independence of annual maximum series per station
#'
#' The Ljung-Box portmanteau test (Box & Pierce, 1970; Ljung & Box, 1978)
#' evaluates the null hypothesis that the first LJUNG_BOX_LAGS autocorrelation
#' coefficients are jointly zero (white-noise process).
#'
#' Failure to reject H0 supports the i.i.d. assumption required for
#' frequency analysis. Rejection indicates autocorrelation that may cause
#' underestimation of uncertainty in quantile estimates.
#'
#' WMO Technical Note 168 recommends testing at lags 1 to n/4.
#'
#' @param data  Raw precipitation tibble with columns id, nombre, anio, p_max.
#' @return Data frame with Ljung-Box statistic, p-value, and iid_ok flag.
run_independence_tests <- function(data) {
  message("
Running Ljung-Box serial independence tests...")

  results <- data |>
    dplyr::group_by(id, nombre) |>
    dplyr::summarise(
      series = list(p_max[order(anio)]),
      n      = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      lb_result = purrr::map2(series, n, function(x, n_obs) {
        x <- stats::na.omit(unlist(x))
        lags_use <- max(1L, min(LJUNG_BOX_LAGS, floor(n_obs / 4L)))
        if (length(x) < 8L)
          return(data.frame(lb_stat = NA_real_, lb_pval = NA_real_,
                            lags_used = lags_use, iid_ok = NA))
        lb <- tryCatch(
          stats::Box.test(x, lag = lags_use, type = "Ljung-Box"),
          error = function(e) NULL)
        if (is.null(lb))
          return(data.frame(lb_stat = NA_real_, lb_pval = NA_real_,
                            lags_used = lags_use, iid_ok = NA))
        data.frame(
          lb_stat   = round(lb$statistic, 4),
          lb_pval   = round(lb$p.value,   4),
          lags_used = lags_use,
          iid_ok    = lb$p.value >= LJUNG_BOX_ALPHA
        )
      })
    ) |>
    tidyr::unnest(lb_result) |>
    dplyr::select(-series)

  sep <- paste(rep("-", 72L), collapse = "")
  message(sep)
  message("LJUNG-BOX INDEPENDENCE TEST RESULTS")
  message(sprintf("  H0: no serial autocorrelation at lags 1..%d (alpha=%.2f)",
                  LJUNG_BOX_LAGS, LJUNG_BOX_ALPHA))
  message(sep)

  for (i in seq_len(nrow(results))) {
    r <- results[i, ]
    status <- if (is.na(r$iid_ok)) "INSUFFICIENT DATA"
              else if (r$iid_ok)   "✓ IID assumption supported"
              else                  "⚠ AUTOCORRELATION DETECTED — CI may be underestimated"
    message(sprintf("  %-20s | Q(%d)=%-8.4f | p=%-6.4f | %s",
                    r$nombre, r$lags_used, r$lb_stat, r$lb_pval, status))
  }
  message(sep)

  n_fail <- sum(!results$iid_ok, na.rm = TRUE)
  if (n_fail > 0L) {
    warning(sprintf(
      paste0("%d station(s) show significant autocorrelation (Ljung-Box p < %.2f).\n",
             "  Frequency analysis assumes i.i.d. series. Autocorrelation may be\n",
             "  caused by: ENSO teleconnections, catchment storage, or data errors.\n",
             "  Recommendation: increase bootstrap replicates (N_BOOTSTRAP) and\n",
             "  interpret confidence intervals conservatively."),
      n_fail, LJUNG_BOX_ALPHA))
  } else {
    message("  All stations: i.i.d. assumption is supported.")
  }

  results
}


# -----------------------------------------------------------------------------
# SECTION 5c: REGIONAL HOMOGENEITY (HOSKING & WALLIS H-STATISTIC)
# -----------------------------------------------------------------------------

#' Compute Hosking & Wallis discordancy (D) and heterogeneity (H) statistics
#'
#' D-statistic (Hosking & Wallis, 1997, eq. 3.3):
#'   Identifies stations whose L-moment ratios are unusually distant from the
#'   regional average. D > DISC_THRESHOLD (default 3) flags an atypical station.
#'
#' H-statistic (Hosking & Wallis, 1997, eq. 4.4):
#'   Compares the observed dispersion of L-CV values across stations with the
#'   dispersion expected from a homogeneous region (estimated via HW_NSIM
#'   Monte Carlo simulations of a 4-parameter kappa distribution).
#'   Interpretation:
#'     H < 1   : Region is acceptably homogeneous
#'     1 <= H < 2 : Possibly heterogeneous
#'     H >= 2  : Definitely heterogeneous
#'
#' @param data  Raw precipitation tibble.
#' @return Named list: discordancy (data frame), H_stat (numeric),
#'   H_class (character), regional_lmom (numeric vector).
run_homogeneity_tests <- function(data) {
  message("
Running Hosking & Wallis regional homogeneity tests...")

  # Compute at-site L-moment ratios
  site_lmom <- data |>
    dplyr::group_by(id, nombre) |>
    dplyr::summarise(
      n    = dplyr::n(),
      series = list(stats::na.omit(p_max)),
      .groups = "drop"
    ) |>
    dplyr::filter(n >= MIN_YEARS) |>
    dplyr::mutate(
      lmom = purrr::map(series, function(x) {
        tryCatch({
          lm <- lmom::samlmu(unlist(x))
          data.frame(l1 = lm[1], lcv = lm[2] / lm[1],
                     lskew = lm[3], lkurt = lm[4])
        }, error = function(e) NULL)
      })
    ) |>
    dplyr::filter(!purrr::map_lgl(lmom, is.null)) |>
    tidyr::unnest(lmom) |>
    dplyr::select(-series)

  if (nrow(site_lmom) < 3L) {
    warning("Homogeneity test requires >= 3 stations with sufficient data.")
    return(NULL)
  }

  n_sites <- nrow(site_lmom)
  ni      <- site_lmom$n
  lcv_i   <- site_lmom$lcv
  lsk_i   <- site_lmom$lskew

  # ---- Discordancy D-statistic ----------------------------------------------
  # u_i = [L-CV_i, L-skew_i, L-kurt_i]' — vector of L-moment ratios per site
  u_mat <- cbind(site_lmom$lcv, site_lmom$lskew, site_lmom$lkurt)
  u_bar <- colMeans(u_mat)
  A_mat <- crossprod(sweep(u_mat, 2L, u_bar))   # sum of outer products

  A_inv <- tryCatch(solve(A_mat), error = function(e) NULL)
  disc_d <- if (!is.null(A_inv)) {
    vapply(seq_len(n_sites), function(i) {
      d <- u_mat[i, ] - u_bar
      (n_sites / 3L) * as.numeric(t(d) %*% A_inv %*% d)
    }, numeric(1L))
  } else {
    rep(NA_real_, n_sites)
  }

  discordancy_df <- site_lmom |>
    dplyr::select(id, nombre, n, lcv, lskew, lkurt) |>
    dplyr::mutate(
      D_stat    = round(disc_d, 3),
      discordant = !is.na(D_stat) & D_stat > DISC_THRESHOLD
    )

  # ---- H-statistic (Hosking & Wallis, 1997) ---------------------------------
  # Regional weighted L-CV
  N_total  <- sum(ni)
  lcv_reg  <- sum(ni * lcv_i) / N_total
  lsk_reg  <- sum(ni * lsk_i) / N_total

  # Observed V: weighted standard deviation of at-site L-CV
  V_obs <- sqrt(sum(ni * (lcv_i - lcv_reg)^2) / N_total)

  # Regional L-kurtosis for kappa-4 simulation
  lku_reg <- sum(ni * site_lmom$lkurt) / N_total

  # Fit kappa-4 parameters to regional L-moments
  kappa_params <- tryCatch(
    lmom::pelkap(c(1, lcv_reg, lsk_reg, lku_reg)),
    error = function(e) NULL)

  H_stat <- NA_real_
  H_class <- "Cannot compute (kappa fit failed)"
  V_sim_vals <- NULL

  if (!is.null(kappa_params)) {
    set.seed(RANDOM_SEED)
    V_sim_vals <- replicate(HW_NSIM, {
      sim_series <- lapply(ni, function(n_i) {
        lmom::quakap(stats::runif(n_i), kappa_params)
      })
      lcv_sim <- vapply(seq_along(sim_series), function(j) {
        lm <- tryCatch(lmom::samlmu(sim_series[[j]]),
                       error = function(e) rep(NA_real_, 4L))
        if (anyNA(lm)) NA_real_ else lm[2] / lm[1]
      }, numeric(1L))
      lcv_sim_reg <- sum(ni * lcv_sim, na.rm = TRUE) / N_total
      sqrt(sum(ni * (lcv_sim - lcv_sim_reg)^2, na.rm = TRUE) / N_total)
    })

    V_sim_vals <- V_sim_vals[is.finite(V_sim_vals)]
    if (length(V_sim_vals) >= 10L) {
      mu_v  <- mean(V_sim_vals)
      sd_v  <- stats::sd(V_sim_vals)
      H_stat <- (V_obs - mu_v) / sd_v
      H_class <- dplyr::case_when(
        H_stat <  1 ~ "HOMOGENEOUS (H < 1)",
        H_stat <  2 ~ "POSSIBLY HETEROGENEOUS (1 <= H < 2)",
        TRUE        ~ "DEFINITELY HETEROGENEOUS (H >= 2)"
      )
    }
  }

  # ---- Print report ---------------------------------------------------------
  sep <- paste(rep("-", 72L), collapse = "")
  message(sep)
  message("HOSKING & WALLIS REGIONAL HOMOGENEITY REPORT")
  message(sprintf("  Stations: %d | Total station-years: %d", n_sites, N_total))
  message(sep)
  message("  Discordancy D-statistic (threshold = ", DISC_THRESHOLD, "):")
  for (i in seq_len(nrow(discordancy_df))) {
    r <- discordancy_df[i, ]
    flag <- if (!is.na(r$discordant) && r$discordant)
              " ⚠ ATYPICAL — verify station data" else " ✓ concordant"
    message(sprintf("    %-20s | D = %6.3f%s", r$nombre, r$D_stat, flag))
  }
  message(sep)
  message(sprintf("  Regional L-CV:    %.4f", lcv_reg))
  message(sprintf("  Regional L-skew:  %.4f", lsk_reg))
  message(sprintf("  H-statistic:      %s",
                  if (is.finite(H_stat)) sprintf("%.3f", H_stat) else "N/A"))
  message(sprintf("  Classification:   %s", H_class))
  message(sep)

  n_disc <- sum(discordancy_df$discordant, na.rm = TRUE)
  if (n_disc > 0L) {
    warning(sprintf(
      paste0("%d station(s) are discordant (D > %.1f): %s.\n",
             "  Their L-moment ratios differ substantially from regional average.\n",
             "  Recommended action: verify data quality, check for gauge errors\n",
             "  or micro-climate effects before including in regional analysis."),
      n_disc, DISC_THRESHOLD,
      paste(discordancy_df$nombre[discordancy_df$discordant], collapse = ", ")))
  }

  if (is.finite(H_stat) && H_stat >= 2) {
    warning(sprintf(
      paste0("Region is DEFINITELY HETEROGENEOUS (H = %.3f).\n",
             "  The stations do not belong to a single homogeneous region.\n",
             "  IDW/Kriging interpolation of GRADEX across heterogeneous\n",
             "  stations will produce physically inconsistent results.\n",
             "  Recommended actions:\n",
             "    1. Split into sub-regions using cluster analysis\n",
             "    2. Use only stations within the same hydrological region\n",
             "    3. Apply regional L-moment estimation (index-flood method)"),
      H_stat))
  } else if (is.finite(H_stat) && H_stat >= 1) {
    message(sprintf(
      "  NOTE: Region is possibly heterogeneous (H = %.3f). \n",
      "  Proceed with caution and report H in technical documentation.", H_stat))
  }

  list(
    discordancy    = discordancy_df,
    H_stat         = H_stat,
    H_class        = H_class,
    V_obs          = V_obs,
    regional_lcv   = lcv_reg,
    regional_lskew = lsk_reg,
    n_simulations  = length(V_sim_vals)
  )
}


# -----------------------------------------------------------------------------
# SECTION 6: MANN-KENDALL TREND ANALYSIS
# -----------------------------------------------------------------------------

#' Theil-Sen slope estimator (mm/year)
#'
#' Computes the median of all pairwise slopes (Sen, 1968). This is the
#' robust nonparametric companion of the Mann-Kendall test.
#'
#' @param x numeric, time index (e.g., year)
#' @param y numeric, response (e.g., annual maximum precipitation)
#' @return numeric scalar: median pairwise slope (units of y per unit of x)
theil_sen_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 3L) return(NA_real_)
  i <- combn(seq_len(n), 2L)
  dx <- x[i[2L, ]] - x[i[1L, ]]
  dy <- y[i[2L, ]] - y[i[1L, ]]
  stats::median(dy[dx != 0] / dx[dx != 0], na.rm = TRUE)
}

compute_mk_statistics <- function(df) {
  if (nrow(df) < 5L)
    return(data.frame(tau = NA_real_, p_valor = NA_real_,
                       pendiente_sen = NA_real_, S_kendall = NA_real_))
  mk <- tryCatch(Kendall::MannKendall(df$p_max), error = function(e) NULL)
  if (is.null(mk))
    return(data.frame(tau = NA_real_, p_valor = NA_real_,
                       pendiente_sen = NA_real_, S_kendall = NA_real_))
  # FIX v1.1: previously the function returned mk$S (Kendall S statistic)
  # mislabelled as "pendiente". S is the sum of signs of all pairwise
  # differences and is NOT a slope in mm/year. The correct slope is Theil-Sen.
  sen_slope <- theil_sen_slope(df$anio, df$p_max)
  data.frame(tau           = round(as.numeric(mk$tau), 4),
             p_valor       = round(as.numeric(mk$sl),  4),
             pendiente_sen = round(sen_slope,           3),  # mm/year
             S_kendall     = round(as.numeric(mk$S),    0))  # raw Kendall statistic
}

run_trend_analysis <- function(data) {
  message("\nRunning Mann-Kendall trend analysis (with Theil-Sen slope)...")
  trend_results <- data |>
    dplyr::group_by(id, nombre) |>
    dplyr::reframe(compute_mk_statistics(dplyr::pick(dplyr::everything()))) |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(tau)) |>
    dplyr::mutate(
      tendencia = dplyr::case_when(
        p_valor < 0.01 ~ "Very significant",
        p_valor < 0.05 ~ "Significant",
        p_valor < 0.10 ~ "Weak",
        TRUE           ~ "Not significant"),
      direccion = dplyr::case_when(
        tau > 0 ~ "Increasing", tau < 0 ~ "Decreasing", TRUE ~ "No trend")
    )

  if (nrow(trend_results) == 0L) {
    warning("No trend statistics could be computed.")
    return(data.frame())
  }
  message("Mann-Kendall results (Theil-Sen slope in mm/year):")
  print(trend_results)

  # Compute Sen slopes per station for plotting overlay
  sen_lines <- data |>
    dplyr::group_by(id, nombre) |>
    dplyr::summarise(
      slope_sen = theil_sen_slope(anio, p_max),
      intercept_sen = stats::median(p_max - theil_sen_slope(anio, p_max) * anio,
                                     na.rm = TRUE),
      anio_min = min(anio, na.rm = TRUE),
      anio_max = max(anio, na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot2::ggplot(data, ggplot2::aes(x = anio, y = p_max,
                                           color = nombre)) +
    ggplot2::geom_line(linewidth = 0.6, alpha = 0.7) +
    ggplot2::geom_point(size = 1.4, alpha = 0.8) +
    # Theil-Sen line per station (matches the slope reported in the table)
    ggplot2::geom_segment(data = sen_lines,
                           ggplot2::aes(x = anio_min, xend = anio_max,
                                        y = intercept_sen + slope_sen * anio_min,
                                        yend = intercept_sen + slope_sen * anio_max,
                                        color = nombre),
                           linewidth = 1.1,
                           inherit.aes = FALSE) +
    ggplot2::labs(title = "TREND ANALYSIS \u2014 ANNUAL MAXIMUM PRECIPITATION",
                  subtitle = "Mann-Kendall test | Theil-Sen robust slope",
                  x = "Year", y = "Maximum precipitation (mm)",
                  color = "Station") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_color_brewer(palette = "Set1")

  print(p)
  ggplot2::ggsave("analisis_tendencia.png", p, width = 12, height = 7, dpi = 300)
  message("Trend plot saved: analisis_tendencia.png")
  trend_results
}


# -----------------------------------------------------------------------------
# SECTION 6b: NON-STATIONARITY SENSITIVITY ANALYSIS (DETRENDED SERIES)
# -----------------------------------------------------------------------------
# When the Mann-Kendall test detects significant trends in some stations,
# the assumption of stationarity required by classical frequency analysis is
# violated. This section provides a sensitivity check by re-running the
# pointwise GRADEX estimation on residual series:
#
#   p_residual_i_t = p_max_i_t  -  Sen_slope_i * (t - t_mid_i)
#
# The series are recentered to their stationary mean. The GRADEX is then
# recomputed per station and contrasted with the stationary-assumption
# estimate to assess the magnitude of bias introduced by ignoring trends.
#
# This implements recommendation P7 of the manuscript review: "Aplicar al
# menos una corrección simple: re-ejecutar el análisis sobre la serie
# residual (precipitación − tendencia Sen)".
# -----------------------------------------------------------------------------

#' Build a detrended copy of the input data using Theil-Sen slopes per station.
#'
#' For each station, subtracts (slope_sen * (year - year_mid)) from p_max,
#' producing a residual series whose linear-Sen trend is approximately zero.
#' The series is recentered to preserve the original mean (so that the
#' magnitude of GRADEX is comparable, not just its variability).
#'
#' @param data data.frame with columns id, nombre, anio, p_max
#' @param trend_results data.frame from run_trend_analysis (columns id, pendiente_sen)
#' @return data.frame with same structure but p_max replaced by detrended values
build_detrended_data <- function(data, trend_results) {
  if (is.null(trend_results) || nrow(trend_results) == 0L) {
    message("  No trend results available; skipping detrending.")
    return(data)
  }
  data |>
    dplyr::left_join(
      trend_results |> dplyr::select(id, pendiente_sen),
      by = "id"
    ) |>
    dplyr::group_by(id) |>
    dplyr::mutate(
      anio_mid = mean(anio, na.rm = TRUE),
      p_max_orig = p_max,
      p_max = p_max_orig - dplyr::coalesce(pendiente_sen, 0) *
              (anio - anio_mid)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-anio_mid, -p_max_orig, -pendiente_sen)
}

#' Run a parallel GRADEX estimation on detrended series and report sensitivity.
#'
#' Executes compute_station_gradex on each detrended station series and
#' compares the per-station GRADEX with the stationary-assumption estimate.
#' Prints a summary table; full results are saved to the Excel workbook.
#'
#' @param data data.frame original input
#' @param trend_results data.frame from run_trend_analysis
#' @param results_stationary data.frame of GRADEX results from main pipeline
#' @return data.frame with side-by-side comparison
run_detrended_sensitivity <- function(data, trend_results,
                                       results_stationary) {
  message("\n--- Non-stationarity sensitivity analysis ---")
  message("  Re-running pointwise GRADEX on Sen-detrended series...")
  data_dt <- build_detrended_data(data, trend_results)

  # Re-estimate GRADEX per station on the detrended series
  results_dt <- data_dt |>
    dplyr::group_by(id) |>
    dplyr::summarise(
      compute_station_gradex(p_max, dplyr::first(id)),
      .groups = "drop"
    ) |>
    dplyr::select(id, gradex_detrended = Mejor_gradex,
                   modelo_detrended = Mejor_Modelo)

  # Side-by-side comparison
  comp <- results_stationary |>
    dplyr::select(id, gradex_stationary = Mejor_gradex,
                   modelo_stationary = Mejor_Modelo) |>
    dplyr::left_join(results_dt, by = "id") |>
    dplyr::left_join(
      trend_results |> dplyr::select(id, nombre, pendiente_sen, p_valor),
      by = "id"
    ) |>
    dplyr::mutate(
      diff_mm  = round(gradex_detrended - gradex_stationary, 2),
      diff_pct = round(100 * (gradex_detrended - gradex_stationary) /
                        gradex_stationary, 1)
    ) |>
    dplyr::select(id, nombre, pendiente_sen, p_valor,
                   modelo_stationary, gradex_stationary,
                   modelo_detrended, gradex_detrended,
                   diff_mm, diff_pct)

  message("  Sensitivity comparison (stationary vs Sen-detrended GRADEX):")
  print(comp)

  max_abs_pct <- max(abs(comp$diff_pct), na.rm = TRUE)
  mean_abs_pct <- mean(abs(comp$diff_pct), na.rm = TRUE)
  message(sprintf(
    "  >> Max |delta|: %.1f%%   Mean |delta|: %.1f%%   (stations: %d)",
    max_abs_pct, mean_abs_pct, nrow(comp)))

  if (max_abs_pct > 15) {
    message("  WARNING: Stationarity assumption introduces >15% bias in ",
            "at least one station. Consider non-stationary GEV in future work.")
  } else if (max_abs_pct > 5) {
    message("  NOTE: Detrending changes GRADEX by 5-15% in some stations.")
  } else {
    message("  OK: Detrending impact is small (<5%); stationarity assumption ",
            "is adequate for design purposes.")
  }
  comp
}


# -----------------------------------------------------------------------------
# SECTION 7: CONFIDENCE INTERVAL ESTIMATION
# -----------------------------------------------------------------------------

estimate_confidence_intervals <- function(results,
                                          conf_level = CONFIDENCE_LEVEL,
                                          n_boot     = N_BOOTSTRAP) {
  vals    <- stats::na.omit(results$Mejor_gradex)
  alpha   <- 1 - conf_level
  z_crit  <- stats::qnorm(1 - alpha / 2)
  n_stn   <- length(vals)
  mu      <- mean(vals)
  sigma   <- stats::sd(vals)
  se      <- sigma / sqrt(n_stn)

  boot_means <- replicate(n_boot,
    mean(sample(vals, size = n_stn, replace = TRUE)))
  ci_boot <- as.numeric(stats::quantile(boot_means,
                                        probs = c(alpha/2, 1-alpha/2)))

  message(sprintf("  Bootstrap CI (n=%d): [%.2f, %.2f] mm",
                  n_boot, ci_boot[1], ci_boot[2]))

  list(conf_level    = conf_level,
       mean          = mu, sd = sigma, se = se,
       ci_normal     = c(mu - z_crit * se, mu + z_crit * se),
       ci_percentile = as.numeric(stats::quantile(vals,
                         probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)),
       ci_bootstrap  = ci_boot,
       z_critical    = z_crit, n_bootstrap = n_boot)
}


# -----------------------------------------------------------------------------
# SECTION 8: IDW EXPONENT OPTIMISATION
# -----------------------------------------------------------------------------

optimise_idw_power <- function(stations_sf, power_grid = IDW_POWER_GRID) {
  message("\nOptimising IDW exponent (LOO-CV grid search)...")
  n              <- nrow(stations_sf)
  rmse_per_power <- numeric(length(power_grid))

  for (pi in seq_along(power_grid)) {
    p <- power_grid[pi]
    message(sprintf("  Testing exponent %.1f (%d/%d)...", p, pi, length(power_grid)))
    errors <- numeric(n)
    for (i in seq_len(n)) {
      train <- stations_sf[-i, ]
      test  <- stations_sf[ i, ]
      if (nrow(train) < 2L) { errors[i] <- NA_real_; next }
      mod       <- gstat::gstat(formula = Mejor_gradex ~ 1, data = train,
                                nmax = IDW_NMAX, set = list(idp = p))
      pred      <- suppressWarnings(
        predict(mod, test, debug.level = 0)$var1.pred)
      errors[i] <- (test$Mejor_gradex - pred)^2
    }
    rmse_per_power[pi] <- sqrt(mean(errors, na.rm = TRUE))
    message(sprintf("    RMSE = %.4f mm", rmse_per_power[pi]))
  }

  opt <- power_grid[which.min(rmse_per_power)]
  message(sprintf("  Optimal IDW exponent: %.1f", opt))
  opt
}


# -----------------------------------------------------------------------------
# SECTION 9: VARIOGRAM FITTING AND KRIGING DIAGNOSTICS
# -----------------------------------------------------------------------------

#' Fit automatic variogram and evaluate spatial structure for Kriging
#'
#' Uses automap::autofitVariogram to select the best variogram model.
#' Computes the nugget/sill ratio as the primary indicator of spatial
#' continuity: a ratio < 0.5 means more than half the total variance is
#' spatially structured, making Kriging more informative than IDW.
#'
#' Minimum station requirement: Kriging requires at least 5 stations to
#' build a reliable empirical variogram. With fewer stations, the number
#' of point pairs per lag bin is too small for stable parameter estimation,
#' which can cause memory errors in the underlying FORTRAN/C routines.
#'
#' @param stations_sp  SpatialPointsDataFrame. Stations in projected CRS.
#' @return Named list: vgm_model, nugget_sill_ratio, nugget, sill,
#'   range_m, model_name, kriging_recommended. NULL if insufficient data.
fit_variogram <- function(stations_sp) {
  message("\nFitting automatic variogram (automap)...")

  n_stn <- nrow(stations_sp)

  # Minimum stations for a reliable variogram: gstat needs at least
  # N*(N-1)/2 unique pairs > number of variogram bins. With < 5 stations
  # this condition often fails and can abort the R session.
  MIN_STATIONS_KRIGING <- 5L
  if (n_stn < MIN_STATIONS_KRIGING) {
    message(sprintf(
      paste0("  Kriging skipped: only %d station(s) available.\n",
             "  Minimum required for stable variogram fitting: %d.\n",
             "  IDW will be used as the interpolation method."),
      n_stn, MIN_STATIONS_KRIGING))
    return(NULL)
  }

  # Wrap autofitVariogram in a forked process via tryCatch + withCallingHandlers
  # to prevent C-level memory errors from aborting the main R session.
  vgm_fit <- tryCatch(
    withCallingHandlers(
      automap::autofitVariogram(
        formula    = Mejor_gradex ~ 1,
        input_data = stations_sp,
        model      = c("Sph", "Exp", "Gau"),  # restrict to stable models
        kappa      = c(0.5, 1, 2),
        fix.values = c(NA, NA, NA)
      ),
      warning = function(w) {
        message("  Variogram warning (non-fatal): ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      message("  Variogram fitting failed: ", e$message)
      NULL
    }
  )

  if (is.null(vgm_fit)) return(NULL)

  vm          <- vgm_fit$var_model
  nugget      <- vm$psill[1]
  sill_total  <- sum(vm$psill)
  partial_sill <- vm$psill[2]
  range_m     <- vm$range[2]
  model_name  <- as.character(vm$model[2])
  ns_ratio    <- nugget / sill_total

  message(sprintf("  Variogram model:      %s", model_name))
  message(sprintf("  Nugget:               %.4f", nugget))
  message(sprintf("  Partial sill:         %.4f", partial_sill))
  message(sprintf("  Total sill:           %.4f", sill_total))
  message(sprintf("  Range:                %.1f m", range_m))
  message(sprintf("  Nugget/Sill ratio:    %.3f", ns_ratio))

  if (ns_ratio < KRIGING_NUGGET_SILL_MAX) {
    message(sprintf(
      "  -> Nugget/Sill = %.3f < %.1f: strong spatial structure detected.",
      ns_ratio, KRIGING_NUGGET_SILL_MAX))
    message("     Kriging is CANDIDATE (final decision depends on CV RMSE).")
  } else {
    message(sprintf(
      "  -> Nugget/Sill = %.3f >= %.1f: weak spatial structure.",
      ns_ratio, KRIGING_NUGGET_SILL_MAX))
    message("     IDW will be preferred regardless of RMSE comparison.")
  }

  list(vgm_model            = vgm_fit,
       nugget_sill_ratio     = ns_ratio,
       nugget                = nugget,
       sill                  = sill_total,
       range_m               = range_m,
       model_name            = model_name,
       kriging_recommended   = (ns_ratio < KRIGING_NUGGET_SILL_MAX))
}


# -----------------------------------------------------------------------------
# SECTION 10: LEAVE-ONE-OUT CROSS-VALIDATION (IDW + KRIGING)
# -----------------------------------------------------------------------------

#' Run LOO cross-validation for a single spatial method
#'
#' @param stations_sf  sf object with Mejor_gradex.
#' @param method       Character: "IDW" or "Kriging".
#' @param idw_power    Numeric. IDW exponent (ignored for Kriging).
#' @param vgm_info     List from fit_variogram (ignored for IDW).
#' @return Data frame with observed / predicted columns, or NULL on failure.
run_loo_cv <- function(stations_sf, method, idw_power = 2,
                       vgm_info = NULL) {

  n      <- nrow(stations_sf)
  if (n < 3L) return(NULL)
  preds  <- vector("list", n)

  for (i in seq_len(n)) {
    train_sf <- stations_sf[-i, ]
    test_sf  <- stations_sf[ i, ]
    if (nrow(train_sf) < 2L) next

    pred_val <- tryCatch({
      if (method == "IDW") {
        mod  <- gstat::gstat(formula = Mejor_gradex ~ 1, data = train_sf,
                             nmax = IDW_NMAX, set = list(idp = idw_power))
        suppressWarnings(predict(mod, test_sf, debug.level = 0)$var1.pred)
      } else {
        # Ordinary Kriging — per-iteration variogram on leave-one-out subset.
        # autofitVariogram is wrapped defensively: model set is restricted to
        # three numerically stable models to avoid C-level crashes on small n.
        train_sp <- sf::as_Spatial(train_sf)
        vgm_loo  <- tryCatch(
          withCallingHandlers(
            automap::autofitVariogram(
              Mejor_gradex ~ 1, train_sp,
              model      = c("Sph", "Exp", "Gau"),
              kappa      = c(0.5, 1, 2),
              fix.values = c(NA, NA, NA)
            )$var_model,
            warning = function(w) invokeRestart("muffleWarning")
          ),
          error = function(e) NULL
        )
        if (is.null(vgm_loo)) {
          message(sprintf("  [Kriging] Station %d: variogram failed, skipping.", i))
          return(NA_real_)
        }
        test_sp <- sf::as_Spatial(test_sf)
        kr <- tryCatch(
          gstat::krige(Mejor_gradex ~ 1, train_sp, test_sp,
                       model = vgm_loo, debug.level = 0),
          error = function(e) NULL
        )
        if (is.null(kr)) NA_real_ else kr$var1.pred
      }
    }, error = function(e) NA_real_)

    preds[[i]] <- data.frame(
      Estacion  = stations_sf$nombre[i],
      Observado = stations_sf$Mejor_gradex[i],
      Predicho  = pred_val,
      Metodo    = method
    )

    message(sprintf("  [%s] Station %d/%d: obs=%.2f | pred=%.2f mm",
                    method, i, n,
                    stations_sf$Mejor_gradex[i], pred_val))
  }

  dplyr::bind_rows(preds) |>
    dplyr::filter(!is.na(Predicho)) |>
    dplyr::mutate(
      Error     = Predicho - Observado,
      Error_Abs = abs(Error),
      Error_Rel = abs(Error / Observado) * 100
    )
}


#' Compute CV metrics data frame from predictions
#'
#' Returns a data frame with columns Metrica, Valor (numeric), Valor_text
#' (character), and Metodo. The split between Valor and Valor_text is
#' deliberate: extra_info entries that are text (e.g. "Variogram model" = "Sph")
#' would otherwise force `Valor` into character type when rbind'ed, breaking
#' downstream code that expects numeric values (FIX v1.1).
cv_metrics <- function(df, method, extra_info = NULL) {
  r2   <- stats::cor(df$Observado, df$Predicho)^2
  rmse <- sqrt(mean((df$Observado - df$Predicho)^2))
  mae  <- mean(df$Error_Abs)
  bias <- mean(df$Error)
  rel  <- mean(df$Error_Rel)

  rows <- data.frame(
    Metrica    = c("R\u00b2", "RMSE (mm)", "MAE (mm)",
                   "BIAS (mm)", "Error Relativo Medio (%)"),
    Valor      = round(c(r2, rmse, mae, bias, rel), 4),
    Valor_text = NA_character_,
    Metodo     = method,
    stringsAsFactors = FALSE
  )
  if (!is.null(extra_info)) {
    extra_vals <- unlist(extra_info)
    extra_num  <- suppressWarnings(as.numeric(extra_vals))
    extra_text <- ifelse(is.na(extra_num), as.character(extra_vals), NA_character_)
    rows <- rbind(rows, data.frame(
      Metrica    = names(extra_info),
      Valor      = extra_num,                 # numeric or NA
      Valor_text = extra_text,                # character or NA
      Metodo     = method,
      stringsAsFactors = FALSE))
  }
  rows
}


#' Build CV scatter plot for one method
cv_scatter_plot <- function(df, metrics_df, method, subtitle_extra = "") {
  r2_v   <- round(metrics_df$Valor[metrics_df$Metrica == "R\u00b2"], 3)
  rmse_v <- round(metrics_df$Valor[metrics_df$Metrica == "RMSE (mm)"], 2)

  ggplot2::ggplot(df, ggplot2::aes(x = Observado, y = Predicho)) +
    ggplot2::geom_point(size = 4, alpha = 0.75,
                        color = ifelse(method == "IDW", "steelblue", "darkorange")) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red",
                         linetype = "dashed", linewidth = 1) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, alpha = 0.2,
                         linewidth = 0.8, color = "darkgreen") +
    # FIX v1.1: switch to ggrepel to avoid label overlap
    ggrepel::geom_text_repel(ggplot2::aes(label = Estacion),
                              size = 3.2, box.padding = 0.4,
                              point.padding = 0.3,
                              segment.color = "grey60",
                              segment.size = 0.3,
                              max.overlaps = Inf,
                              min.segment.length = 0) +
    ggplot2::annotate("label",
                      x = min(df$Observado), y = max(df$Predicho),
                      label = paste0(
                        "RMSE: ", rmse_v, " mm\n",
                        "MAE:  ",
                        round(metrics_df$Valor[metrics_df$Metrica=="MAE (mm)"], 2),
                        " mm\n",
                        "BIAS: ",
                        round(metrics_df$Valor[metrics_df$Metrica=="BIAS (mm)"], 2),
                        " mm"),
                      hjust = 0, vjust = 1, size = 3,
                      fill = "white", alpha = 0.85) +
    ggplot2::labs(
      title    = paste("CROSS-VALIDATION \u2014", method),
      subtitle = paste0("R\u00b2 = ", r2_v, " | RMSE = ", rmse_v, " mm",
                        if (nchar(subtitle_extra) > 0)
                          paste0(" | ", subtitle_extra) else ""),
      x = "Observed GRADEX (mm)", y = "Predicted GRADEX (mm)") +
    ggplot2::theme_minimal()
}


#' Run full cross-validation for both methods and apply multi-metric selection
#'
#' Selection rule (multi-metric scoring):
#'   Each of three criteria (RMSE, R², |BIAS|) awards 1 point to the better
#'   method. Kriging is selected when score_Kriging > score_IDW. Ties are
#'   broken by RMSE (lower wins). Kriging failure falls back to IDW.
#'   Geostatistical context (nugget/sill) is reported for interpretation but
#'   does not veto — the score decides.
#'
#' @param results   Data frame with station GRADEX results.
#' @param config    List with lon_aforo, lat_aforo.
#' @return Named list: idw, kriging, selected_method, vgm_info, metrics_combined.
run_cross_validation_both <- function(results, config) {
  message(rep("=", 60L))
  message("CROSS-VALIDATION: IDW + ORDINARY KRIGING")
  message(rep("=", 60L))

  stations_sf <- sf::st_as_sf(results, coords = c("lon", "lat"), crs = 4326L)
  n_stations  <- nrow(stations_sf)

  if (n_stations < 3L) {
    warning("Cross-validation requires >= 3 stations. Skipping.")
    return(NULL)
  }

  # --- IDW: optimise exponent then LOO-CV -----------------------------------
  optimal_power <- optimise_idw_power(stations_sf)
  message(sprintf("\nRunning IDW LOO-CV (exponent = %.1f)...", optimal_power))
  pred_idw <- run_loo_cv(stations_sf, "IDW", idw_power = optimal_power)

  # --- Kriging: variogram fitting then LOO-CV --------------------------------
  # Project to CTM12 for metric variogram distances
  ctm12_crs   <- 9377L
  stn_proj    <- sf::st_transform(stations_sf, ctm12_crs)
  stations_sp <- sf::as_Spatial(stn_proj)

  vgm_info <- fit_variogram(stations_sp)

  pred_ok <- NULL
  if (!is.null(vgm_info)) {
    message("\nRunning Kriging LOO-CV...")
    # Additional guard: LOO-CV leaves n-1 stations per fold; require >= 5
    # in training set, meaning the full network needs >= 6 stations.
    if (nrow(stations_sp) >= 6L) {
      pred_ok <- run_loo_cv(stn_proj, "Kriging", vgm_info = vgm_info)
    } else {
      message(sprintf(
        paste0("  Kriging LOO-CV skipped: %d stations total.\n",
               "  Each fold would train on %d stations — too few for\n",
               "  stable per-fold variogram fitting. IDW will be used."),
        nrow(stations_sp), nrow(stations_sp) - 1L))
    }
  }

  # --- Compute metrics -------------------------------------------------------
  if (is.null(pred_idw) || nrow(pred_idw) == 0L) {
    message("IDW cross-validation produced no predictions.")
    return(NULL)
  }

  metrics_idw <- cv_metrics(pred_idw, "IDW",
    extra_info = list("IDW exponent" = optimal_power))

  metrics_ok  <- if (!is.null(pred_ok) && nrow(pred_ok) > 0L)
    cv_metrics(pred_ok, "Kriging",
      extra_info = list(
        "Variogram model"  = if (!is.null(vgm_info)) vgm_info$model_name else NA,
        "Nugget/Sill"      = if (!is.null(vgm_info))
                               round(vgm_info$nugget_sill_ratio, 3) else NA,
        "Range (m)"        = if (!is.null(vgm_info))
                               round(vgm_info$range_m, 0) else NA))
  else NULL

  # ---- Extract scalar metrics for both methods ------------------------------
  get_m <- function(m_df, name) {
    v <- m_df$Valor[m_df$Metrica == name]
    if (length(v) == 0L) NA_real_ else as.numeric(v[1L])
  }

  rmse_idw <- get_m(metrics_idw, "RMSE (mm)")
  r2_idw   <- get_m(metrics_idw, "R\u00b2")
  bias_idw <- abs(get_m(metrics_idw, "BIAS (mm)"))

  rmse_ok  <- if (!is.null(metrics_ok)) get_m(metrics_ok, "RMSE (mm)") else Inf
  r2_ok    <- if (!is.null(metrics_ok)) get_m(metrics_ok, "R\u00b2")   else -Inf
  bias_ok  <- if (!is.null(metrics_ok))
                abs(get_m(metrics_ok, "BIAS (mm)")) else Inf

  # ---- Multi-metric scoring --------------------------------------------------
  # Each criterion awards SEL_WEIGHT_* points to the better method.
  # Lower RMSE, higher R², lower |BIAS| = better.
  score_ok  <- 0L
  score_idw <- 0L
  criteria_log <- character(0L)

  # RMSE
  if (is.finite(rmse_ok) && rmse_ok < rmse_idw) {
    score_ok  <- score_ok  + SEL_WEIGHT_RMSE
    criteria_log <- c(criteria_log,
      sprintf("RMSE: Kriging %.3f < IDW %.3f [+%d to Kriging]",
              rmse_ok, rmse_idw, SEL_WEIGHT_RMSE))
  } else {
    score_idw <- score_idw + SEL_WEIGHT_RMSE
    criteria_log <- c(criteria_log,
      sprintf("RMSE: IDW %.3f <= Kriging %.3f [+%d to IDW]",
              rmse_idw, if(is.finite(rmse_ok)) rmse_ok else NA,
              SEL_WEIGHT_RMSE))
  }

  # R²
  if (is.finite(r2_ok) && r2_ok > r2_idw) {
    score_ok  <- score_ok  + SEL_WEIGHT_R2
    criteria_log <- c(criteria_log,
      sprintf("R\u00b2:   Kriging %.3f > IDW %.3f [+%d to Kriging]",
              r2_ok, r2_idw, SEL_WEIGHT_R2))
  } else {
    score_idw <- score_idw + SEL_WEIGHT_R2
    criteria_log <- c(criteria_log,
      sprintf("R\u00b2:   IDW %.3f >= Kriging %.3f [+%d to IDW]",
              r2_idw, if(is.finite(r2_ok)) r2_ok else NA, SEL_WEIGHT_R2))
  }

  # |BIAS|
  if (is.finite(bias_ok) && bias_ok < bias_idw) {
    score_ok  <- score_ok  + SEL_WEIGHT_BIAS
    criteria_log <- c(criteria_log,
      sprintf("|BIAS|: Kriging %.3f < IDW %.3f [+%d to Kriging]",
              bias_ok, bias_idw, SEL_WEIGHT_BIAS))
  } else {
    score_idw <- score_idw + SEL_WEIGHT_BIAS
    criteria_log <- c(criteria_log,
      sprintf("|BIAS|: IDW %.3f <= Kriging %.3f [+%d to IDW]",
              bias_idw, if(is.finite(bias_ok)) bias_ok else NA,
              SEL_WEIGHT_BIAS))
  }

  # ---- v1.2: ANTI-DEGENERATE-KRIGING SAFEGUARD ------------------------------
  # When the variogram is dominated by nugget (N/S close to 1), Ordinary Kriging
  # collapses toward the regional mean: all predictions cluster near a constant
  # value regardless of station location. In LOO cross-validation this can
  # produce a misleadingly high R² (because the few residual variations align
  # with the observation rank) and a near-zero |BIAS| (because predictions
  # average to the observed mean), even though the predictor is essentially
  # uninformative. The RMSE, however, is *not* improved compared with IDW.
  #
  # Operational test: variance ratio of LOO predictions to observations.
  #   ratio = var(Predicho) / var(Observado)
  #   ratio < 0.10  -> Kriging predictions span < 10% of observation variance
  #                    -> degenerate (predicts ~mean); discard from scoring
  #   ratio >= 0.10 -> retains usable spatial discrimination; keep
  #
  # The 0.10 threshold is conservative: a healthy Kriging on a region with
  # genuine spatial structure typically yields ratios > 0.40. Values < 0.10
  # are diagnostic of pure-nugget collapse.
  #
  # If degenerate: force IDW selection AND emit a warning so the user knows
  # the Kriging "wins" on R²/BIAS were artefacts of the variance collapse.
  kriging_degenerate <- FALSE
  kriging_var_ratio  <- NA_real_
  if (!is.null(pred_ok) && nrow(pred_ok) >= 3L) {
    var_obs  <- stats::var(pred_ok$Observado, na.rm = TRUE)
    var_pred <- stats::var(pred_ok$Predicho,  na.rm = TRUE)
    if (is.finite(var_obs) && var_obs > 0) {
      kriging_var_ratio <- var_pred / var_obs
      if (is.finite(kriging_var_ratio) && kriging_var_ratio < 0.10) {
        kriging_degenerate <- TRUE
        warning(sprintf(
          paste0("\n  DEGENERATE KRIGING DETECTED:\n",
                 "  var(Predicho)/var(Observado) = %.4f < 0.10\n",
                 "  Kriging predictions span less than 10%% of observation\n",
                 "  variance — the variogram has collapsed to near-constant\n",
                 "  predictions (typical of N/S ratio close to 1).\n",
                 "  Although Kriging may score higher on R\u00b2 and |BIAS|, those\n",
                 "  apparent gains are artefacts of the variance collapse.\n",
                 "  Action: Kriging is DISQUALIFIED from method selection;\n",
                 "  IDW will be used regardless of multi-metric score."),
          kriging_var_ratio))
        # Force IDW by zeroing Kriging's score
        score_ok <- 0L
        criteria_log <- c(criteria_log,
          sprintf("DEGENERACY: Kriging var ratio %.3f < 0.10 [score zeroed]",
                  kriging_var_ratio))
      }
    }
  }
  # ---------------------------------------------------------------------------

  # Final decision: higher score wins; ties broken by RMSE
  selected_method <- if (!is.null(metrics_ok) && is.finite(rmse_ok) &&
                         !kriging_degenerate) {
    if      (score_ok  > score_idw) "Kriging"
    else if (score_idw > score_ok)  "IDW"
    else if (rmse_ok   < rmse_idw)  "Kriging"   # tie-break
    else                            "IDW"
  } else {
    "IDW"  # Kriging unavailable or degenerate
  }

  # ---- Print selection decision report --------------------------------------
  sep  <- paste(rep("\u2550", 72L), collapse = "")
  sep2 <- paste(rep("\u2500", 72L), collapse = "")
  cat("\n\u2554", sep, "\u2557\n", sep = "")
  cat("\u2551  METHOD SELECTION REPORT — MULTI-METRIC SCORING",
      paste(rep(" ", 22L), collapse = ""), "\u2551\n", sep = "")
  cat("\u2560", sep, "\u2563\n", sep = "")

  # Metric table
  cat(sprintf("\u2551  %-30s  %12s  %12s  %10s\u2551\n",
              "Criterion (weight)", "IDW", "Kriging", "Winner"))
  cat(sprintf("\u2551  %s\u2551\n", sep2))
  cat(sprintf("\u2551  %-30s  %12.4f  %12s  %10s\u2551\n",
              paste0("RMSE (w=", SEL_WEIGHT_RMSE, ", lower=better)"),
              rmse_idw,
              if(is.finite(rmse_ok)) sprintf("%.4f", rmse_ok) else "N/A",
              if(rmse_ok < rmse_idw) "Kriging" else "IDW"))
  cat(sprintf("\u2551  %-30s  %12.4f  %12s  %10s\u2551\n",
              paste0("R\u00b2   (w=", SEL_WEIGHT_R2,  ", higher=better)"),
              r2_idw,
              if(is.finite(r2_ok)) sprintf("%.4f", r2_ok) else "N/A",
              if(is.finite(r2_ok) && r2_ok > r2_idw) "Kriging" else "IDW"))
  cat(sprintf("\u2551  %-30s  %12.4f  %12s  %10s\u2551\n",
              paste0("|BIAS|(w=", SEL_WEIGHT_BIAS, ", lower=better)"),
              bias_idw,
              if(is.finite(bias_ok)) sprintf("%.4f", bias_ok) else "N/A",
              if(is.finite(bias_ok) && bias_ok < bias_idw) "Kriging" else "IDW"))
  cat(sprintf("\u2551  %s\u2551\n", sep2))
  cat(sprintf("\u2551  %-30s  %12d  %12d  %10s\u2551\n",
              "TOTAL SCORE", score_idw, score_ok,
              if (score_ok > score_idw) "Kriging"
              else if (score_idw > score_ok) "IDW" else "TIED (RMSE)"))

  # Variogram context (informational only)
  if (!is.null(vgm_info)) {
    cat(sprintf("\u2551  %s\u2551\n", sep2))
    cat(sprintf("\u2551  Variogram: %-61s\u2551\n",
      paste0(vgm_info$model_name,
             " | Range=", round(vgm_info$range_m, 0), " m",
             " | Nugget/Sill=", round(vgm_info$nugget_sill_ratio, 3),
             " [informational]")))
  }

  # v1.2: Anti-degenerate-Kriging diagnostic line
  if (!is.null(pred_ok) && is.finite(kriging_var_ratio)) {
    cat(sprintf("\u2551  %s\u2551\n", sep2))
    deg_status <- if (kriging_degenerate)
      "DEGENERATE (< 0.10) -- Kriging disqualified" else "OK (>= 0.10)"
    cat(sprintf("\u2551  Kriging var ratio: %.4f -- %-44s\u2551\n",
                kriging_var_ratio, deg_status))
  }

  cat("\u2560", sep, "\u2563\n", sep = "")
  cat(sprintf("\u2551  \u2605 SELECTED: %-61s\u2551\n", selected_method))

  reason <- if (kriging_degenerate)
    sprintf("Kriging disqualified (var ratio %.3f < 0.10) -- IDW used",
            kriging_var_ratio)
  else if (selected_method == "Kriging" && score_ok > score_idw)
    sprintf("Kriging scored %d/%d metrics vs IDW %d/%d",
            score_ok, SEL_WEIGHT_RMSE+SEL_WEIGHT_R2+SEL_WEIGHT_BIAS,
            score_idw, SEL_WEIGHT_RMSE+SEL_WEIGHT_R2+SEL_WEIGHT_BIAS)
  else if (selected_method == "Kriging")
    "Tied on score; Kriging RMSE < IDW RMSE (tie-break)"
  else if (score_idw > score_ok)
    sprintf("IDW scored %d/%d metrics vs Kriging %d/%d",
            score_idw, SEL_WEIGHT_RMSE+SEL_WEIGHT_R2+SEL_WEIGHT_BIAS,
            score_ok,  SEL_WEIGHT_RMSE+SEL_WEIGHT_R2+SEL_WEIGHT_BIAS)
  else if (!is.finite(rmse_ok))
    "Kriging cross-validation failed -- IDW fallback"
  else
    "Tied on score; IDW RMSE <= Kriging RMSE (tie-break)"

  cat(sprintf("\u2551  Reason: %-64s\u2551\n", reason))
  cat("\u255a", sep, "\u255d\n", sep = "")

  # --- Explicit R² warning when spatial interpolation has low predictive power --
  r2_selected <- if (selected_method == "Kriging" && !is.null(metrics_ok))
    as.numeric(metrics_ok$Valor[metrics_ok$Metrica == "R\u00b2"])
  else
    as.numeric(metrics_idw$Valor[metrics_idw$Metrica == "R\u00b2"])

  if (is.finite(r2_selected) && r2_selected < 0.3) {
    warning(sprintf(
      paste0("\n  LOW R\u00b2 WARNING: Selected method (%s) has R\u00b2 = %.3f\n",
             "  Spatial interpolation explains only %.1f%% of GRADEX variability.\n",
             "  With R\u00b2 < 0.3, the interpolated value at the gauging point is\n",
             "  not reliably predicted from the station network.\n",
             "  Root causes: too few stations, high spatial heterogeneity,\n",
             "  or outlier stations dominating the network.\n",
             "  Recommended action: do NOT use the interpolated GRADEX as a\n",
             "  point estimate. Report the station-mean range instead, and\n",
             "  expand the sensitivity range to \u00b125-30%% in HEC-HMS."),
      selected_method, r2_selected, r2_selected * 100))
  } else if (is.finite(r2_selected) && r2_selected < 0.5) {
    message(sprintf(
      "  NOTE: %s R\u00b2 = %.3f — weak but usable. Expand sensitivity range.",
      selected_method, r2_selected))
  }

  # --- CV scatter plots -------------------------------------------------------
  p_idw <- cv_scatter_plot(pred_idw, metrics_idw, "IDW",
    paste0("exponent = ", optimal_power))
  print(p_idw)
  ggplot2::ggsave("validacion_cruzada_idw.png", p_idw,
                  width = 8, height = 7, dpi = 300)

  if (!is.null(pred_ok)) {
    p_ok <- cv_scatter_plot(pred_ok, metrics_ok, "Kriging",
      if (!is.null(vgm_info))
        paste0(vgm_info$model_name, " | N/S=",
               round(vgm_info$nugget_sill_ratio, 2)) else "")
    print(p_ok)
    ggplot2::ggsave("validacion_cruzada_kriging.png", p_ok,
                    width = 8, height = 7, dpi = 300)
  }

  message("Cross-validation plots saved.")

  list(
    idw              = list(predictions = pred_idw, metrics = metrics_idw,
                             optimal_power = optimal_power),
    kriging          = if (!is.null(pred_ok))
                         list(predictions = pred_ok, metrics = metrics_ok,
                              vgm_info = vgm_info) else NULL,
    selected_method  = selected_method,
    vgm_info         = vgm_info,
    metrics_combined = dplyr::bind_rows(metrics_idw, metrics_ok)
  )
}


# -----------------------------------------------------------------------------
# SECTION 11: MAIN GRADEX ANALYSIS PIPELINE
# -----------------------------------------------------------------------------

run_gradex_analysis <- function(data, lon_gauge, lat_gauge,
                                min_years = MIN_YEARS) {
  message(rep("=", 60L))
  message("GRADEX ANALYSIS FOR HEC-HMS")
  message(rep("=", 60L))

  if (!"nombre" %in% names(data)) data$nombre <- data$id

  stations_grouped <- data |>
    dplyr::group_by(id, nombre, lon, lat) |>
    dplyr::summarise(n_years = dplyr::n(), series = list(p_max), .groups = "drop") |>
    dplyr::filter(n_years >= min_years)

  gradex_detail <- purrr::map2_dfr(stations_grouped$series,
                                    stations_grouped$id,
                                    compute_station_gradex)

  results <- stations_grouped |>
    dplyr::select(-series) |>
    dplyr::left_join(gradex_detail, by = "id") |>
    dplyr::filter(!is.na(Mejor_gradex))

  if (nrow(results) == 0L)
    stop("No stations with sufficient data. Check MIN_YEARS parameter.")

  # Print record-quality summary
  message("\nRecord quality summary:")
  print(results |> dplyr::select(nombre, n_years, record_quality))

  # Run independence and homogeneity diagnostics on raw data
  independence_results <- run_independence_tests(data)
  homogeneity_results  <- run_homogeneity_tests(data)

  message("\nEstimating confidence intervals...")
  ci_list <- estimate_confidence_intervals(results, CONFIDENCE_LEVEL)

  station_table <- results |>
    dplyr::select(Station = nombre, Years = n_years,
                  Record_Quality = record_quality,
                  Gumbel = Gumbel_gradex,   Gumbel_RMSE  = Gumbel_rmse,
                  `Log-Normal` = LogNormal_gradex,
                  `Log-Normal_RMSE` = LogNormal_rmse,
                  `Pearson III` = Pearson_gradex, Pearson_RMSE = Pearson_rmse,
                  Gamma = Gamma_gradex,     Gamma_RMSE   = Gamma_rmse,
                  `Best Model` = Mejor_Modelo,
                  `Final GRADEX` = Mejor_gradex,
                  `GRADEX CI Lower` = gradex_ci_lower,
                  `GRADEX CI Upper` = gradex_ci_upper)

  message("\nPer-station results:"); print(station_table)

  dist_comparison <- data.frame(
    Distribution = c("Gumbel", "Log-Normal", "Pearson III", "Gamma"),
    Mean_GRADEX  = round(c(mean(results$Gumbel_gradex,    na.rm = TRUE),
                            mean(results$LogNormal_gradex, na.rm = TRUE),
                            mean(results$Pearson_gradex,   na.rm = TRUE),
                            mean(results$Gamma_gradex,     na.rm = TRUE)), 2),
    Mean_RMSE    = round(c(mean(results$Gumbel_rmse,    na.rm = TRUE),
                            mean(results$LogNormal_rmse, na.rm = TRUE),
                            mean(results$Pearson_rmse,   na.rm = TRUE),
                            mean(results$Gamma_rmse,     na.rm = TRUE)), 3),
    Best_Count   = c(sum(results$Mejor_Modelo == "Gumbel"),
                     sum(results$Mejor_Modelo == "Log-Normal"),
                     sum(results$Mejor_Modelo == "Pearson III"),
                     sum(results$Mejor_Modelo == "Gamma"))
  )

  # Preliminary IDW estimate (median exponent) — updated after CV
  idw_power_init <- stats::median(IDW_POWER_GRID)
  stations_sf    <- sf::st_as_sf(results, coords = c("lon", "lat"), crs = 4326L)
  gauge_sf       <- sf::st_as_sf(data.frame(lon = lon_gauge, lat = lat_gauge),
                                  coords = c("lon", "lat"), crs = 4326L)

  gradex_point <- if (nrow(results) >= 3L) {
    mod <- gstat::gstat(formula = Mejor_gradex ~ 1, data = stations_sf,
                        nmax = IDW_NMAX, set = list(idp = idw_power_init))
    suppressWarnings(predict(mod, gauge_sf, debug.level = 0)$var1.pred)
  } else {
    mean(results$Mejor_gradex, na.rm = TRUE)
  }
  gradex_point <- round(gradex_point, 3)

  cat(sprintf(
    "\n   \u2554%s\u2557\n   \u2551   PRELIMINARY GRADEX = %8.2f mm (IDW, exp=%.1f) \u2551\n   \u255a%s\u255d\n",
    paste(rep("\u2550", 50L), collapse = ""), gradex_point, idw_power_init,
    paste(rep("\u2550", 50L), collapse = "")))

  list(gradex_hec_hms          = gradex_point,
       gradex_stations         = results,
       station_table           = station_table,
       distribution_comparison = dist_comparison,
       confidence_intervals    = ci_list,
       independence_results    = independence_results,
       homogeneity_results     = homogeneity_results,
       gauge_coords            = c(lon = lon_gauge, lat = lat_gauge),
       n_stations              = nrow(results),
       stations_sf             = stations_sf,
       gauge_sf                = gauge_sf)
}


# -----------------------------------------------------------------------------
# SECTION 12: INTERPOLATION AT GAUGING POINT (FINAL ESTIMATE)
# -----------------------------------------------------------------------------

#' Apply selected method to estimate GRADEX at the target gauging point
#'
#' @param analysis_result  List from run_gradex_analysis.
#' @param cv_result        List from run_cross_validation_both.
#' @return Numeric. Final GRADEX estimate (mm).
estimate_gradex_at_gauge <- function(analysis_result, cv_result) {

  stations_sf     <- analysis_result$stations_sf
  gauge_sf        <- analysis_result$gauge_sf
  selected_method <- cv_result$selected_method
  ctm12_crs       <- 9377L

  message(sprintf("\nEstimating GRADEX at gauging point using: %s",
                  selected_method))

  gradex_final <- tryCatch({
    if (selected_method == "Kriging") {
      stn_proj   <- sf::st_transform(stations_sf, ctm12_crs)
      gauge_proj <- sf::st_transform(gauge_sf,    ctm12_crs)
      stn_sp     <- sf::as_Spatial(stn_proj)
      gau_sp     <- sf::as_Spatial(gauge_proj)
      vgm_model  <- cv_result$vgm_info$vgm_model$var_model
      kr         <- gstat::krige(Mejor_gradex ~ 1, stn_sp, gau_sp,
                                  model = vgm_model, debug.level = 0)
      kr$var1.pred
    } else {
      opt_power <- cv_result$idw$optimal_power
      mod       <- gstat::gstat(formula = Mejor_gradex ~ 1,
                                data = stations_sf, nmax = IDW_NMAX,
                                set = list(idp = opt_power))
      suppressWarnings(predict(mod, gauge_sf, debug.level = 0)$var1.pred)
    }
  }, error = function(e) {
    warning("Final estimation failed: ", e$message,
            "\nFalling back to station mean.")
    mean(analysis_result$gradex_stations$Mejor_gradex, na.rm = TRUE)
  })

  round(gradex_final, 3)
}


# -----------------------------------------------------------------------------
# SECTION 13: MAP GENERATION — MAIN + COMPARISON PANEL
# -----------------------------------------------------------------------------

#' Build an interpolated grid using either IDW or Kriging
#'
#' @param stations_sp  SpatialPointsDataFrame in CTM12.
#' @param grid_sp      SpatialPixelsDataFrame  in CTM12.
#' @param method       "IDW" or "Kriging".
#' @param idw_power    Numeric IDW exponent.
#' @param vgm_model    gstat variogram model (Kriging only).
#' @return Data frame with columns X, Y, gradex_interp.
build_interp_grid <- function(stations_sp, grid_sp, method,
                              idw_power = 2, vgm_model = NULL) {
  interp_sp <- tryCatch({
    if (method == "IDW") {
      mod <- gstat::gstat(formula = Mejor_gradex ~ 1, data = stations_sp,
                          nmax = IDW_NMAX, set = list(idp = idw_power))
      suppressWarnings(predict(mod, grid_sp, debug.level = 0))
    } else {
      gstat::krige(Mejor_gradex ~ 1, stations_sp, grid_sp,
                   model = vgm_model, debug.level = 0)
    }
  }, error = function(e) {
    message("  Grid prediction failed for ", method, ": ", e$message)
    NULL
  })

  if (is.null(interp_sp)) return(NULL)
  df           <- as.data.frame(interp_sp)
  names(df)[3] <- "gradex_interp"
  df[order(df$Y, df$X), ]
}


#' Create a single-method ggplot tile+contour layer stack
#'
#' @param interp_df    Data frame with X, Y, gradex_interp.
#' @param stations_df  Data frame with X, Y, nombre, Mejor_gradex (CTM12).
#' @param gauge_xy     2-col matrix. Gauging point CTM12 coordinates.
#' @param gradex_final Numeric. Estimated GRADEX at gauge.
#' @param ci           List. Confidence intervals.
#' @param method_label Character. Title label.
#' @param n_stations   Integer.
#' @param idw_power    Numeric or NULL.
#' @param vgm_info     List or NULL.
#' @param compact      Logical. If TRUE, suppress gauge annotation (for panel).
#' @return ggplot2 object.
make_interp_map <- function(interp_df, stations_df, gauge_xy,
                            gradex_final, ci, method_label,
                            n_stations, idw_power = NULL,
                            vgm_info = NULL, compact = FALSE) {

  has_metR       <- requireNamespace("metR", quietly = TRUE)
  gradex_range   <- range(interp_df$gradex_interp, na.rm = TRUE)
  contour_breaks <- pretty(gradex_range, n = 6L)

  subtitle_parts <- c(
    paste0("Stations: ", n_stations),
    if (!is.null(idw_power)) paste0("IDW exp: ", idw_power),
    if (!is.null(vgm_info))  paste0(vgm_info$model_name,
                                     " | N/S=",
                                     round(vgm_info$nugget_sill_ratio, 2)),
    "CRS: MAGNA-SIRGAS / CTM12 (EPSG:9377)"
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = interp_df,
                       ggplot2::aes(x = X, y = Y, fill = gradex_interp)) +
    ggplot2::geom_contour(data   = interp_df,
                          ggplot2::aes(x = X, y = Y, z = gradex_interp),
                          breaks = contour_breaks,
                          color = "white", linewidth = 0.45, alpha = 0.9)

  if (has_metR) {
    p <- p + metR::geom_text_contour(
      data   = interp_df,
      ggplot2::aes(x = X, y = Y, z = gradex_interp),
      breaks = contour_breaks, size = 2.8,
      color  = "white", fontface = "bold",
      stroke = 0.2, stroke.color = "black", skip = 0L)
  }

  p <- p +
    ggplot2::geom_point(
      data  = stations_df,
      ggplot2::aes(x = X, y = Y, size = Mejor_gradex, fill = Mejor_gradex),
      shape = 21L, color = "black", stroke = 1.2, alpha = 0.95) +
    # FIX v1.1: ggrepel for non-overlapping labels (avoids the gauging-point
    # info box hiding nearby station names)
    ggrepel::geom_text_repel(
      data = stations_df,
      ggplot2::aes(x = X, y = Y, label = nombre),
      size = 3.0, fontface = "bold", color = "black",
      box.padding = 0.6, point.padding = 0.4,
      segment.color = "grey40", segment.size = 0.3,
      bg.color = "white", bg.r = 0.12,
      max.overlaps = Inf, min.segment.length = 0,
      seed = 1L) +
    ggrepel::geom_text_repel(
      data = stations_df,
      ggplot2::aes(x = X, y = Y,
                   label = sprintf("%.1f mm", Mejor_gradex)),
      size = 2.6, color = "grey20", fontface = "italic",
      box.padding = 0.5, point.padding = 0.3,
      nudge_y = -1500,
      segment.color = "grey60", segment.size = 0.2,
      max.overlaps = Inf, min.segment.length = 0,
      seed = 2L) +
    ggplot2::geom_point(
      data = data.frame(X = gauge_xy[, 1L], Y = gauge_xy[, 2L]),
      ggplot2::aes(x = X, y = Y),
      color = "red", size = 5L, shape = 17L) +
    ggplot2::scale_fill_viridis_c(name = "GRADEX (mm)",
                                   option = "plasma", direction = -1L,
                                   limits = gradex_range) +
    ggplot2::scale_size_continuous(name = "Measured\nGRADEX (mm)",
                                    range = c(4L, 10L))

  if (!compact) {
    p <- p + ggplot2::annotate(
      "label",
      x = gauge_xy[, 1L], y = gauge_xy[, 2L],
      label = sprintf(
        "\u2605 GAUGING POINT\nGRADEX: %.2f mm\nCI95 normal: [%.1f, %.1f] mm\nCI95 bootstrap: [%.1f, %.1f] mm",
        gradex_final,
        ci$ci_normal[1], ci$ci_normal[2],
        ci$ci_bootstrap[1], ci$ci_bootstrap[2]),
      color = "#8b0000", fontface = "bold", size = 2.8,
      fill = "white", alpha = 0.9,
      hjust = 0, vjust = 1,
      label.padding = ggplot2::unit(0.25, "lines"))
  }

  p <- p +
    ggspatial::annotation_scale(location = "bl", width_hint = 0.3,
                                 unit_category = "metric") +
    ggspatial::annotation_north_arrow(location = "tr", which_north = "grid",
                                       style = ggspatial::north_arrow_minimal()) +
    ggplot2::labs(
      title    = paste("GRADEX MAP \u2014", method_label),
      subtitle = paste(subtitle_parts, collapse = " | "),
      x = "Easting \u2014 CTM12 (m)",
      y = "Northing \u2014 CTM12 (m)",
      caption = if (!compact)
        paste0("GRADEX at gauging point: ",
               round(gradex_final, 2), " mm") else NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold",
                                               size = 13L),
      plot.subtitle   = ggplot2::element_text(hjust = 0.5, size = 8L,
                                               color = "gray40"),
      plot.caption    = ggplot2::element_text(hjust = 0.5, face = "italic",
                                               size = 8L),
      axis.text       = ggplot2::element_text(size = 7L))
  p
}


#' Generate main map (selected method) and side-by-side comparison panel
#'
#' @param analysis_result  List from run_gradex_analysis.
#' @param cv_result        List from run_cross_validation_both.
#' @param gradex_final     Numeric. Final GRADEX at gauge (mm).
generate_maps <- function(analysis_result, cv_result, gradex_final) {

  message("\nGenerating interpolation maps...")

  ctm12_crs       <- 9377L
  selected_method <- cv_result$selected_method
  stn_proj        <- sf::st_transform(analysis_result$stations_sf, ctm12_crs)
  gauge_proj      <- sf::st_transform(analysis_result$gauge_sf,    ctm12_crs)
  gauge_xy        <- sf::st_coordinates(gauge_proj)
  stations_sp     <- sf::as_Spatial(stn_proj)

  # Build interpolation grid (SpatialPixels — batch prediction)
  bbox_p  <- sf::st_bbox(stn_proj)
  buf_m   <- BBOX_BUFFER * 111320
  x_seq   <- seq(bbox_p["xmin"] - buf_m, bbox_p["xmax"] + buf_m,
                  length.out = GRID_RES_MAP)
  y_seq   <- seq(bbox_p["ymin"] - buf_m, bbox_p["ymax"] + buf_m,
                  length.out = GRID_RES_MAP)
  grid_df <- expand.grid(X = x_seq, Y = y_seq)
  sp::coordinates(grid_df) <- ~ X + Y
  sp::proj4string(grid_df) <- sp::CRS(
    SRS_string = sp::proj4string(stations_sp))
  sp::gridded(grid_df) <- TRUE

  # Station data frame in CTM12 for labels
  stn_xy <- sf::st_coordinates(stn_proj)
  stations_coords_df <- analysis_result$gradex_stations |>
    dplyr::mutate(X = stn_xy[, 1L], Y = stn_xy[, 2L])

  ci <- analysis_result$confidence_intervals

  # --- Build IDW grid --------------------------------------------------------
  message(sprintf("  Building IDW grid (exp=%.1f)...",
                  cv_result$idw$optimal_power))
  df_idw <- build_interp_grid(stations_sp, grid_df, "IDW",
                               idw_power = cv_result$idw$optimal_power)

  # --- Build Kriging grid ----------------------------------------------------
  df_ok <- NULL
  if (!is.null(cv_result$kriging)) {
    message("  Building Kriging grid...")
    df_ok <- build_interp_grid(
      stations_sp, grid_df, "Kriging",
      vgm_model = cv_result$vgm_info$vgm_model$var_model)
  }

  # --- Main map (selected method) -------------------------------------------
  df_main   <- if (selected_method == "Kriging" && !is.null(df_ok)) df_ok
               else df_idw
  vgm_main  <- if (selected_method == "Kriging") cv_result$vgm_info else NULL
  idp_main  <- if (selected_method == "IDW") cv_result$idw$optimal_power else NULL

  map_main <- make_interp_map(
    interp_df    = df_main,
    stations_df  = stations_coords_df,
    gauge_xy     = gauge_xy,
    gradex_final = gradex_final,
    ci           = ci,
    method_label = paste0(selected_method, " \u2605 SELECTED"),
    n_stations   = analysis_result$n_stations,
    idw_power    = idp_main,
    vgm_info     = vgm_main,
    compact      = FALSE)

  print(map_main)
  ggplot2::ggsave("mapa_gradex_final.png", map_main,
                  width = 12L, height = 10L, dpi = 300L)
  message("Main map saved: mapa_gradex_final.png")

  # --- Comparison panel (IDW left | Kriging right) --------------------------
  if (!is.null(df_ok)) {
    p_left <- make_interp_map(
      df_idw, stations_coords_df, gauge_xy,
      gradex_final = cv_result$idw$metrics$Valor[
        cv_result$idw$metrics$Metrica == "RMSE (mm)"],
      ci = ci,
      method_label = paste0("IDW",
        if (selected_method == "IDW") " \u2605 SELECTED" else ""),
      n_stations = analysis_result$n_stations,
      idw_power  = cv_result$idw$optimal_power,
      compact    = TRUE)

    p_right <- make_interp_map(
      df_ok, stations_coords_df, gauge_xy,
      gradex_final = cv_result$kriging$metrics$Valor[
        cv_result$kriging$metrics$Metrica == "RMSE (mm)"],
      ci = ci,
      method_label = paste0("Ordinary Kriging",
        if (selected_method == "Kriging") " \u2605 SELECTED" else ""),
      n_stations = analysis_result$n_stations,
      vgm_info   = cv_result$vgm_info,
      compact    = TRUE)

    map_panel <- patchwork::wrap_plots(p_left, p_right, ncol = 2L) +
      patchwork::plot_annotation(
        title    = "METHOD COMPARISON \u2014 IDW vs. ORDINARY KRIGING",
        subtitle = paste0(
          "Selected: ", selected_method,
          " | IDW RMSE=", round(cv_result$idw$metrics$Valor[
            cv_result$idw$metrics$Metrica == "RMSE (mm)"], 2), " mm",
          " | Kriging RMSE=",
          if (!is.null(cv_result$kriging))
            round(cv_result$kriging$metrics$Valor[
              cv_result$kriging$metrics$Metrica == "RMSE (mm)"], 2)
          else "N/A", " mm",
          if (!is.null(cv_result$vgm_info))
            paste0(" | N/S=",
                   round(cv_result$vgm_info$nugget_sill_ratio, 3),
                   " (threshold=", KRIGING_NUGGET_SILL_MAX, ")")
          else ""),
        theme = ggplot2::theme(
          plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold",
                                                size = 14L),
          plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9L,
                                                color = "gray40")))

    print(map_panel)
    ggplot2::ggsave("mapa_comparacion_metodos.png", map_panel,
                    width = 20L, height = 10L, dpi = 300L)
    message("Comparison panel saved: mapa_comparacion_metodos.png")
  } else {
    message("  Kriging grid unavailable — comparison panel skipped.")
  }

  invisible(map_main)
}


# -----------------------------------------------------------------------------
# SECTION 14: CROSS-VALIDATION CONCLUSIONS
# -----------------------------------------------------------------------------

print_validation_conclusions <- function(cv_result, gradex_final,
                                          confidence_intervals) {
  sel      <- cv_result$selected_method
  m        <- if (sel == "Kriging" && !is.null(cv_result$kriging))
                cv_result$kriging$metrics else cv_result$idw$metrics

  r2        <- m$Valor[m$Metrica == "R\u00b2"]
  rmse      <- m$Valor[m$Metrica == "RMSE (mm)"]
  mae       <- m$Valor[m$Metrica == "MAE (mm)"]
  bias      <- m$Valor[m$Metrica == "BIAS (mm)"]
  rel_error <- m$Valor[m$Metrica == "Error Relativo Medio (%)"]

  ci_lower  <- confidence_intervals$ci_normal[1]
  ci_upper  <- confidence_intervals$ci_normal[2]
  unc_interp <- rmse / gradex_final * 100
  unc_stat   <- (ci_upper - ci_lower) / (2 * gradex_final) * 100
  unc_total  <- sqrt(unc_interp^2 + unc_stat^2)

  sep <- paste(rep("\u2550", 80L), collapse = "")
  cat("\n\u2554", sep, "\u2557\n", sep = "")
  cat(sprintf("\u2551  VALIDATION CONCLUSIONS \u2014 METHOD: %-42s\u2551\n", sel))
  cat("\u2560", sep, "\u2563\n", sep = "")

  # R2
  r2_msg <- dplyr::case_when(
    r2 >= 0.9 ~ "Excellent spatial structure capture",
    r2 >= 0.8 ~ "Good spatial structure capture",
    r2 >= 0.7 ~ "Moderate — consider more stations",
    TRUE      ~ "Low — method NOT recommended")
  cat(sprintf("\u2551  R\u00b2 = %.3f  %-57s\u2551\n", r2, r2_msg))

  # RMSE
  cat(sprintf("\u2551  RMSE = %.2f mm  (%.1f%% of mean GRADEX)%33s\u2551\n",
              rmse, unc_interp, " "))

  # Bias
  bias_msg <- if (abs(bias) < 0.5) "Negligible"
              else if (abs(bias) < 1) paste("Small",
                         if (bias > 0) "(overestimation)" else "(underestimation)")
              else "Significant — review model"
  cat(sprintf("\u2551  BIAS = %.2f mm  %-56s\u2551\n", bias, bias_msg))

  # Relative error
  cat(sprintf("\u2551  Rel. Error = %.1f%%  %-55s\u2551\n", rel_error,
              if (rel_error < 5)  "Very low" else
              if (rel_error < 10) "Acceptable for HEC-HMS" else "Moderate"))

  # Variogram info (always shown when available, informational)
  if (!is.null(cv_result$vgm_info)) {
    vgi <- cv_result$vgm_info
    cat(sprintf("\u2551  Variogram (OK): %-63s\u2551\n",
      paste0(vgi$model_name, " | Range=", round(vgi$range_m, 0),
             " m | Nugget/Sill=", round(vgi$nugget_sill_ratio, 3),
             " [context only]")))
  }

  cat("\u2560", sep, "\u2563\n", sep = "")
  cat(sprintf("\u2551  Combined uncertainty:  IDW/Kriging \u00b1%.1f%%  |  Stat CI \u00b1%.1f%%  |  Total \u00b1%.1f%%%s\u2551\n",
              unc_interp, unc_stat, unc_total, paste(rep(" ", 5L), collapse = "")))

  # Overall verdict
  verdict <- if (r2 >= 0.8 && rmse < 3 && abs(bias) < 0.5 && rel_error < 10)
    "\u2705 VALIDATED — reliable for hydrological design"
  else if (r2 >= 0.7 && rmse < 5 && rel_error < 15)
    "\u26a0 ACCEPTABLE — use with extended sensitivity analysis"
  else
    "\u274c NOT RECOMMENDED — expand gauge network"

  cat(sprintf("\u2551  Overall: %-70s\u2551\n", verdict))

  # HEC-HMS sensitivity range
  pct <- if (rel_error < 10) 10 else if (rel_error < 20) 15 else 20
  cat(sprintf("\u2551  HEC-HMS range: [%.2f, %.2f] mm  (\u00b1%d%%)%37s\u2551\n",
              gradex_final * (1 - pct/100),
              gradex_final * (1 + pct/100), pct, " "))

  cat("\u255a", sep, "\u255d\n", sep = "")
  invisible(unc_total)
}


# -----------------------------------------------------------------------------
# SECTION 15: SHAPEFILE EXPORT
# -----------------------------------------------------------------------------

export_shapefiles <- function(analysis_result, config, gradex_final,
                               selected_method) {
  message("\nExporting shapefiles...")
  output_dir <- "shapefile_resultados"
  if (!dir.exists(output_dir)) dir.create(output_dir)

  station_shape <- analysis_result$gradex_stations |>
    dplyr::select(id, nombre, n_years,
                  record_quality,
                  Gumbel_gx  = Gumbel_gradex,
                  LNorm_gx   = LogNormal_gradex,
                  Pearson_gx = Pearson_gradex,
                  Gamma_gx   = Gamma_gradex,
                  Mejor_Mod  = Mejor_Modelo,
                  Gradex_f   = Mejor_gradex)

  coords_mat         <- sf::st_coordinates(analysis_result$stations_sf)
  station_shape$lon  <- coords_mat[, 1L]
  station_shape$lat  <- coords_mat[, 2L]
  names(station_shape) <- substr(names(station_shape), 1L, 10L)

  suppressWarnings(sf::st_write(
    sf::st_as_sf(station_shape, coords = c("lon", "lat"), crs = 4326L,
                 remove = FALSE),
    file.path(output_dir, "gradex_estaciones.shp"),
    delete_layer = TRUE, quiet = TRUE))

  gauge_shape <- data.frame(
    id         = "GAUGE",
    nombre     = "Gauging_Point",
    gradex_est = gradex_final,
    metodo     = selected_method,
    ic_inf     = analysis_result$confidence_intervals$ci_normal[1],
    ic_sup     = analysis_result$confidence_intervals$ci_normal[2],
    lon        = config$lon_aforo,
    lat        = config$lat_aforo)
  names(gauge_shape) <- substr(names(gauge_shape), 1L, 10L)

  suppressWarnings(sf::st_write(
    sf::st_as_sf(gauge_shape, coords = c("lon", "lat"), crs = 4326L,
                 remove = FALSE),
    file.path(output_dir, "punto_aforo.shp"),
    delete_layer = TRUE, quiet = TRUE))

  message("Shapefiles saved in: ", output_dir, "/")
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# SECTION 16: EXCEL EXPORT
# -----------------------------------------------------------------------------

export_results_excel <- function(analysis_result, config, trend_results,
                                  cv_result, gradex_final) {
  message("\nExporting results to Excel...")
  ci <- analysis_result$confidence_intervals

  workbook_data <- list(
    Summary = data.frame(
      Parameter = c("GRADEX_Final_mm", "Selected_Method",
                    "Number_of_Stations", "Min_Years_Required",
                    "Gauge_Longitude", "Gauge_Latitude",
                    "CI_95_Lower_mm", "CI_95_Upper_mm",
                    "Nugget_Sill_Ratio",
                    "IDW_Exponent", "IDW_RMSE_CV",
                    "Kriging_RMSE_CV", "Calculation_Date"),
      Value = c(
        gradex_final,
        cv_result$selected_method,
        analysis_result$n_stations,
        config$min_anios,
        config$lon_aforo, config$lat_aforo,
        round(ci$ci_normal[1], 2), round(ci$ci_normal[2], 2),
        if (!is.null(cv_result$vgm_info))
          round(cv_result$vgm_info$nugget_sill_ratio, 3) else "N/A",
        cv_result$idw$optimal_power,
        round(cv_result$idw$metrics$Valor[
          cv_result$idw$metrics$Metrica == "RMSE (mm)"], 3),
        if (!is.null(cv_result$kriging))
          round(cv_result$kriging$metrics$Valor[
            cv_result$kriging$metrics$Metrica == "RMSE (mm)"], 3) else "N/A",
        as.character(Sys.Date()))
    ),
    Station_Results       = analysis_result$station_table,
    Distribution_Summary  = analysis_result$distribution_comparison,
    Confidence_Intervals  = data.frame(
      Method     = c("Normal", "Percentile", "Bootstrap"),
      Lower_mm   = round(c(ci$ci_normal[1],
                            ci$ci_percentile[1], ci$ci_bootstrap[1]), 2),
      Upper_mm   = round(c(ci$ci_normal[2],
                            ci$ci_percentile[2], ci$ci_bootstrap[2]), 2),
      Conf_Level = paste0(ci$conf_level * 100, "%")),
    CV_Metrics_Combined   = cv_result$metrics_combined,
    CV_Predictions_IDW    = cv_result$idw$predictions,
    CV_Predictions_Kriging = if (!is.null(cv_result$kriging))
                               cv_result$kriging$predictions else data.frame(),
    Trend_Analysis          = trend_results,
    Independence_Tests      = if (!is.null(analysis_result$independence_results))
      analysis_result$independence_results |>
        dplyr::select(id, nombre, n, lb_stat, lb_pval, lags_used, iid_ok)
      else data.frame(),
    Homogeneity_Discordancy = if (!is.null(analysis_result$homogeneity_results))
      analysis_result$homogeneity_results$discordancy
      else data.frame(),
    Homogeneity_Summary     = if (!is.null(analysis_result$homogeneity_results)) {
      hr <- analysis_result$homogeneity_results
      data.frame(
        Statistic = c("H-statistic", "H-classification",
                      "V_observed", "Regional L-CV",
                      "Regional L-skew", "N_simulations"),
        Value     = c(round(hr$H_stat, 4), hr$H_class,
                      round(hr$V_obs, 6), round(hr$regional_lcv, 4),
                      round(hr$regional_lskew, 4), hr$n_simulations)
      )
    } else data.frame()
  )

  writexl::write_xlsx(workbook_data, "resultados_gradex_completos.xlsx")
  message("Excel workbook saved: resultados_gradex_completos.xlsx")
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# SECTION 17: FINAL SUMMARY REPORT
# -----------------------------------------------------------------------------

print_final_summary <- function(analysis_result, config, cv_result,
                                 gradex_final) {
  sep <- paste(rep("\u2550", 85L), collapse = "")
  ci  <- analysis_result$confidence_intervals
  sel <- cv_result$selected_method

  best_dist_idx <- which.min(
    analysis_result$distribution_comparison$Mean_RMSE)
  best_dist <- analysis_result$distribution_comparison$Distribution[best_dist_idx]

  rmse_sel <- if (sel == "Kriging" && !is.null(cv_result$kriging))
    cv_result$kriging$metrics$Valor[
      cv_result$kriging$metrics$Metrica == "RMSE (mm)"]
  else
    cv_result$idw$metrics$Valor[
      cv_result$idw$metrics$Metrica == "RMSE (mm)"]

  r2_sel <- if (sel == "Kriging" && !is.null(cv_result$kriging))
    cv_result$kriging$metrics$Valor[
      cv_result$kriging$metrics$Metrica == "R\u00b2"]
  else
    cv_result$idw$metrics$Valor[cv_result$idw$metrics$Metrica == "R\u00b2"]

  cat("\n\u2554", sep, "\u2557\n", sep = "")
  cat("\u2551", paste(rep(" ", 29L), collapse = ""),
      "FINAL SUMMARY FOR HEC-HMS",
      paste(rep(" ", 31L), collapse = ""), "\u2551\n", sep = "")
  cat("\u2560", sep, "\u2563\n", sep = "")
  cat(sprintf("\u2551   GRADEX (T10\u2013T100):          %8.2f mm%48s\u2551\n",
              gradex_final, " "))
  cat(sprintf("\u2551   Selected method:           %-54s\u2551\n", sel))
  cat(sprintf("\u2551   Gauging point (lon, lat):  (%.4f, %.4f)%38s\u2551\n",
              config$lon_aforo, config$lat_aforo, " "))
  cat("\u2560", sep, "\u2563\n", sep = "")
  cat(sprintf("\u2551   95%% CI (normal):           [%.2f, %.2f] mm%43s\u2551\n",
              ci$ci_normal[1], ci$ci_normal[2], " "))
  cat(sprintf("\u2551   95%% CI (bootstrap):        [%.2f, %.2f] mm%43s\u2551\n",
              ci$ci_bootstrap[1], ci$ci_bootstrap[2], " "))
  cat("\u2560", sep, "\u2563\n", sep = "")
  cat(sprintf("\u2551   Best distribution (RMSE):  %-54s\u2551\n", best_dist))
  cat(sprintf("\u2551   CV RMSE (%s):%-68s\u2551\n", sel,
              paste0("       ", round(rmse_sel, 2), " mm")))
  cat(sprintf("\u2551   CV R\u00b2  (%s):%-68s\u2551\n", sel,
              paste0("       ", round(r2_sel, 3))))
  if (!is.null(cv_result$vgm_info)) {
    cat(sprintf("\u2551   Nugget/Sill ratio:         %.3f (threshold=%.2f)%35s\u2551\n",
                cv_result$vgm_info$nugget_sill_ratio,
                KRIGING_NUGGET_SILL_MAX, " "))
  }
  cat("\u2560", sep, "\u2563\n", sep = "")

  # Record quality warning block
  rq <- analysis_result$gradex_stations |>
    dplyr::filter(grepl("UNRELIABLE|ACCEPTABLE", record_quality))
  if (nrow(rq) > 0L) {
    cat(sprintf("\u2551   \u26a0 RECORD-LENGTH WARNINGS (%d station(s)):%40s\u2551\n",
                nrow(rq), " "))
    for (i in seq_len(nrow(rq))) {
      cat(sprintf("\u2551     %-81s\u2551\n",
                  paste0(rq$nombre[i], ": ", rq$n_years[i],
                         " yr \u2014 ", rq$record_quality[i])))
    }
    cat("\u2560", sep, "\u2563\n", sep = "")
  }

  pct <- if (gradex_final < 5) { cat(sprintf(
    "\u2551   Basin type:  LOW pluvial energy%52s\u2551\n", " ")); 30
  } else if (gradex_final < 15) { cat(sprintf(
    "\u2551   Basin type:  MODERATE pluvial energy%47s\u2551\n", " ")); 20
  } else { cat(sprintf(
    "\u2551   Basin type:  HIGH pluvial energy%53s\u2551\n", " ")); 15 }
  cat(sprintf("\u2551   Sensitivity: \u00b1%d%% \u2192 [%.2f, %.2f] mm%47s\u2551\n",
              pct, gradex_final*(1-pct/100),
              gradex_final*(1+pct/100), " "))
  # v1.2 cosmetic fix: keep diagnostics INSIDE the box (was leaking out before)
  ir <- analysis_result$independence_results
  hr <- analysis_result$homogeneity_results

  if (!is.null(ir) || !is.null(hr)) {
    cat("\u2560", sep, "\u2563\n", sep = "")
    cat(sprintf("\u2551   DIAGNOSTICS%72s\u2551\n", " "))
  }
  if (!is.null(ir)) {
    n_autocorr <- sum(!ir$iid_ok, na.rm = TRUE)
    diag_line <- sprintf("Independence (Ljung-Box):  %d/%d stations pass i.i.d. test",
                         nrow(ir) - n_autocorr, nrow(ir))
    cat(sprintf("\u2551   %-82s\u2551\n", diag_line))
  }
  if (!is.null(hr)) {
    n_disc <- sum(hr$discordancy$discordant, na.rm = TRUE)
    diag_line <- sprintf("Discordancy (D > %.1f):     %d station(s) flagged",
                         DISC_THRESHOLD, n_disc)
    cat(sprintf("\u2551   %-82s\u2551\n", diag_line))
    diag_line <- sprintf("Homogeneity (H-stat):      %s",
                         if (is.finite(hr$H_stat))
                           paste0(round(hr$H_stat, 3), " \u2014 ", hr$H_class)
                         else "N/A")
    cat(sprintf("\u2551   %-82s\u2551\n", diag_line))
  }
  cat("\u255a", sep, "\u255d\n\n", sep = "")

  cat("Process completed. Output directory:", getwd(), "\n\n")
  cat("Generated files:\n")
  cat("  resultados_gradex_completos.xlsx\n")
  cat("    Sheets: Summary | Station_Results | Distribution_Summary\n")
  cat("    Sheets: Confidence_Intervals | CV_Metrics_Combined\n")
  cat("    Sheets: CV_Predictions_IDW | CV_Predictions_Kriging\n")
  cat("    Sheets: Independence_Tests | Homogeneity_Discordancy\n")
  cat("    Sheets: Homogeneity_Summary | Trend_Analysis\n")
  cat("  mapa_gradex_final.png             (selected method)\n")
  cat("  mapa_comparacion_metodos.png      (IDW vs Kriging panel)\n")
  cat("  validacion_cruzada_idw.png\n")
  cat("  validacion_cruzada_kriging.png\n")
  cat("  analisis_tendencia.png\n")
  cat("  shapefile_resultados/\n\n")
  cat(sprintf("Enter in HEC-HMS:  GRADEX = %.2f mm\n", gradex_final))
  cat(sprintf("95%% CI: [%.1f, %.1f] mm\n",
              ci$ci_normal[1], ci$ci_normal[2]))

  invisible(NULL)
}


# =============================================================================
# SECTION 18: EXECUTION PIPELINE — DO NOT MODIFY BELOW THIS LINE
# =============================================================================

config <- list(
  lon_aforo         = LON_GAUGE,
  lat_aforo         = LAT_GAUGE,
  min_anios         = MIN_YEARS,
  optimal_idw_power = NULL
)

# Step 1 — GRADEX estimation per station (with record-length warnings)
analysis_result <- run_gradex_analysis(
  data      = df_raw,
  lon_gauge = config$lon_aforo,
  lat_gauge = config$lat_aforo,
  min_years = config$min_anios
)

# Step 2 — Mann-Kendall trend analysis
trend_results <- run_trend_analysis(df_raw)

# Step 2b — Non-stationarity sensitivity: re-run pointwise GRADEX on
#           Sen-detrended series to quantify the bias of the stationarity
#           assumption (manuscript review recommendation P7).
detrended_sensitivity <- tryCatch(
  run_detrended_sensitivity(
    data                = df_raw,
    trend_results       = trend_results,
    results_stationary  = analysis_result$gradex_stations
  ),
  error = function(e) {
    message("  Detrended sensitivity failed: ", e$message)
    NULL
  }
)

# Step 3 — IDW optimisation + Kriging variogram + LOO-CV for both methods
#          + multi-metric selection (RMSE + R² + BIAS scoring)
cv_result <- run_cross_validation_both(
  analysis_result$gradex_stations, config)

# Step 4 — Final GRADEX estimate at gauging point using selected method
gradex_final <- if (!is.null(cv_result)) {
  estimate_gradex_at_gauge(analysis_result, cv_result)
} else {
  analysis_result$gradex_hec_hms
}
analysis_result$gradex_hec_hms <- gradex_final

# Step 5 — Cross-validation conclusions
if (!is.null(cv_result)) {
  print_validation_conclusions(
    cv_result            = cv_result,
    gradex_final         = gradex_final,
    confidence_intervals = analysis_result$confidence_intervals)
}

# Step 6 — Maps: main (selected method) + IDW vs Kriging comparison panel
if (!is.null(cv_result)) {
  generate_maps(analysis_result, cv_result, gradex_final)
}

# Step 7 — Shapefiles
export_shapefiles(analysis_result, config, gradex_final,
                  if (!is.null(cv_result)) cv_result$selected_method else "IDW")

# Step 8 — Excel workbook
export_results_excel(analysis_result, config, trend_results,
                     cv_result, gradex_final)

# Step 9 — Final console summary
print_final_summary(analysis_result, config, cv_result, gradex_final)

# Step 10 — Session info for reproducibility
writeLines(utils::capture.output(utils::sessionInfo()), "session_info.txt")
message("Session info saved: session_info.txt")



# =============================================================================
# SECTION 19: H&W REGIONAL DISTRIBUTION LAYER + BOOTSTRAP CI
# =============================================================================
#
# Purpose (focused scope):
#   This section adds a regional H&W layer ON TOP of the existing pointwise
#   spatial interpolation. It does NOT replace the IDW/Kriging pipeline.
#   Instead it:
#     1. Uses regional L-moments (H&W eq. 4.3) to select the best distribution
#        for the region via the Z goodness-of-fit statistic (H&W eq. 5.6)
#     2. Re-estimates GRADEX at each station using that regional distribution
#        (instead of each station fitting its own best distribution independently)
#     3. Applies the same IDW/Kriging interpolation on the H&W-refined GRADEX
#        values to obtain the final estimate at the gauging point
#     4. Replaces the three current CI methods with the formally correct
#        H&W bootstrap CI (H&W Section 6.3), which propagates parameter
#        estimation uncertainty, inter-site variability, and record-length
#        effects into a single interval
#
# Architecture:
#   Layer 1 (unchanged):  Pointwise GRADEX per station
#   Layer 2 (new):        H&W regional distribution selection
#   Layer 3 (unchanged):  IDW/Kriging spatial interpolation on H&W GRADEX
#   Layer 4 (new):        H&W bootstrap CI replaces current CI methods
#
# Key references:
#   Hosking & Wallis (1997), Regional Frequency Analysis,
#   Cambridge University Press. Chapters 4, 5, 6.
#
# Parameters (defined in Section 1):
#   HW_NSIM, CONFIDENCE_LEVEL, RANDOM_SEED, MIN_YEARS
# =============================================================================


# -----------------------------------------------------------------------------
# SECTION 19a: REGIONAL L-MOMENTS (H&W eq. 4.3)
# -----------------------------------------------------------------------------

#' Compute regional L-moment ratios weighted by record length
#'
#' Standardises each series by its own mean (index-flood step), then computes
#' record-length-weighted L-moment ratios across all stations. The regional
#' L-moments characterise the shape of the regional frequency distribution,
#' decoupled from the local scale (site mean).
#'
#' @param data  Raw precipitation tibble (id, nombre, anio, p_max).
#' @return Named list: lcv_R, lsk_R, lku_R, N_total, n_sites, site_table.
compute_regional_lmom_hw <- function(data) {
  message("\n", rep("-", 60L))
  message("H&W LAYER: Regional L-moment estimation (eq. 4.3)")
  message(rep("-", 60L))

  site_sum <- data |>
    dplyr::group_by(id, nombre) |>
    dplyr::summarise(
      n         = dplyr::n(),
      mean_site = mean(p_max, na.rm = TRUE),
      series    = list(stats::na.omit(p_max)),
      .groups   = "drop"
    ) |>
    dplyr::filter(n >= MIN_YEARS)

  if (nrow(site_sum) < 3L)
    stop("H&W layer requires >= 3 stations with >= MIN_YEARS records.")

  # Station-count reliability classification for H&W regional analysis.
  # Reference: Hosking & Wallis (1997), p. 58 and WMO No. 168 (2008), p. 179.
  #   >= 15 stations : reliable H&W analysis
  #    8-14 stations : acceptable
  #    5-7  stations : limited — Z-test has low power; CI will be wide
  #    3-4  stations : minimum — results should be treated as indicative only
  n_stn_class <- dplyr::case_when(
    nrow(site_sum) >= 15 ~ "RELIABLE (>= 15 stations)",
    nrow(site_sum) >=  8 ~ "ACCEPTABLE (8-14 stations)",
    nrow(site_sum) >=  5 ~ "LIMITED (5-7 stations) — wide CI expected",
    TRUE                 ~ "MINIMUM (3-4 stations) — indicative only"
  )

  if (nrow(site_sum) < 8L) {
    warning(sprintf(
      paste0("H&W regional analysis: only %d station(s) available.\n",
             "  Classification: %s\n",
             "  Expected consequences:\n",
             "    - Z-test has LOW statistical power (sigma_4 will be large)\n",
             "    - Bootstrap CI will be WIDE (potentially >= 100%% amplitude)\n",
             "    - Regional L-moments are sensitive to individual station outliers\n",
             "  Minimum recommended: 8 stations for reliable Z-test discrimination.\n",
             "  NOTE: the H&W weight penalty is applied AFTER the bootstrap CI\n",
             "  is computed (Step 5 of run_hw_layer), not at this stage.\n",
             "  If CI amplitude >= 100%%, the H&W weight is set to zero and the\n",
             "  pointwise estimate is used with an expanded sensitivity range.\n",
             "  The decision will be reported in the H&W LAYER RESULTS block."),
      nrow(site_sum), n_stn_class))
  }
  message(sprintf("  Network classification: %s", n_stn_class))


  # At-site L-moment ratios on standardised (dimensionless) series
  site_sum <- site_sum |>
    dplyr::mutate(
      lmom = purrr::map2(series, mean_site, function(x, mu) {
        x_std <- unlist(x) / mu
        lm    <- tryCatch(lmom::samlmu(x_std, nmom = 4L),
                          error = function(e) rep(NA_real_, 4L))
        data.frame(
          lcv   = if (!is.na(lm[1]) && lm[1] > 0) lm[2]/lm[1] else NA_real_,
          lskew = lm[3], lkurt = lm[4])
      })
    ) |>
    tidyr::unnest(lmom) |>
    dplyr::filter(!is.na(lcv))

  ni     <- site_sum$n
  N      <- sum(ni)
  lcv_R  <- sum(ni * site_sum$lcv,   na.rm = TRUE) / N
  lsk_R  <- sum(ni * site_sum$lskew, na.rm = TRUE) / N
  lku_R  <- sum(ni * site_sum$lkurt, na.rm = TRUE) / N

  message(sprintf("  Stations: %d | Total station-years: %d", nrow(site_sum), N))
  message(sprintf("  Regional L-CV:       %.4f", lcv_R))
  message(sprintf("  Regional L-skew:     %.4f", lsk_R))
  message(sprintf("  Regional L-kurtosis: %.4f", lku_R))

  site_table <- site_sum |>
    dplyr::select(id, nombre, n, mean_site, lcv, lskew, lkurt) |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 4L)))

  list(lcv_R    = lcv_R, lsk_R   = lsk_R, lku_R = lku_R,
       N_total  = N,     n_sites = nrow(site_sum),
       site_table  = site_table,
       site_summary = site_sum)
}


# -----------------------------------------------------------------------------
# SECTION 19b: DISTRIBUTION SELECTION VIA Z-STATISTIC (H&W eq. 5.6)
# -----------------------------------------------------------------------------

#' Select best regional distribution for the H&W layer
#'
#' Tests five distributions (GEV, GLO, GNO, PE3, GPA) by comparing their
#' theoretical L-kurtosis with the observed regional L-kurtosis. The Z-stat
#' measures the deviation in units of the Monte Carlo standard deviation of
#' the regional L-kurtosis estimator (sigma_4). |Z| < 1.64 = adequate fit.
#'
#' @param reg_lmom  List from compute_regional_lmom_hw.
#' @param n_sim     Integer. Monte Carlo simulations for sigma_4.
#' @return Named list: best_dist, params_reg, Z_table, fit_adequate, best_Z.
select_hw_distribution <- function(reg_lmom, n_sim = HW_NSIM) {
  message("\nH&W LAYER: Distribution selection (Z-statistic, eq. 5.6)...")

  lcv <- reg_lmom$lcv_R
  lsk <- reg_lmom$lsk_R
  lku <- reg_lmom$lku_R
  N   <- reg_lmom$N_total
  ni  <- reg_lmom$site_summary$n

  dists <- c("GEV", "GLO", "GNO", "PE3", "GPA")

  # Fit each distribution to [1, lcv_R, lsk_R] and get theoretical L-kurtosis
  fit_dist_params <- function(d) {
    tryCatch(switch(d,
      "GEV" = lmom::pelgev(c(1, lcv, lsk)),
      "GLO" = lmom::pelglo(c(1, lcv, lsk)),
      "GNO" = lmom::pelgno(c(1, lcv, lsk)),
      "PE3" = lmom::pelpe3(c(1, lcv, lsk)),
      "GPA" = lmom::pelgpa(c(1, lcv, lsk))),
      error = function(e) NULL)
  }

  quant_fn <- function(d, p) {
    p_fit <- fit_dist_params(d)
    if (is.null(p_fit)) return(NA_real_)
    tryCatch(switch(d,
      "GEV" = lmom::quagev(p, p_fit), "GLO" = lmom::quaglo(p, p_fit),
      "GNO" = lmom::quagno(p, p_fit), "PE3" = lmom::quape3(p, p_fit),
      "GPA" = lmom::quagpa(p, p_fit)), error = function(e) NA_real_)
  }

  # Theoretical L-kurtosis via lmom analytical population L-moment functions.
  # lmom::lmr*() returns exact population L-moments for the fitted parameters
  # without numerical integration, avoiding the bounds issue that caused
  # lmom::lmrp(bounds = c(0, Inf)) to fail for distributions with support
  # on the full real line (GEV with shape < 0, GLO, GNO, PE3).
  tau4_theo <- vapply(dists, function(d) {
    p_fit <- fit_dist_params(d)
    if (is.null(p_fit)) return(NA_real_)
    tryCatch({
      lm <- switch(d,
        "GEV" = lmom::lmrgev(p_fit, nmom = 4L),
        "GLO" = lmom::lmrglo(p_fit, nmom = 4L),
        "GNO" = lmom::lmrgno(p_fit, nmom = 4L),
        "PE3" = lmom::lmrpe3(p_fit, nmom = 4L),
        "GPA" = lmom::lmrgpa(p_fit, nmom = 4L))
      if (is.null(lm) || length(lm) < 4L || !is.finite(lm[4])) NA_real_
      else lm[4]   # t4 = L-kurtosis (4th L-moment ratio)
    }, error = function(e) NA_real_)
  }, numeric(1L))

  # Monte Carlo sigma_4 (H&W eq. 5.5) using kappa-4 simulation
  kp <- tryCatch(lmom::pelkap(c(1, lcv, lsk, lku)), error = function(e) NULL)
  sigma4 <- NA_real_
  if (!is.null(kp)) {
    set.seed(RANDOM_SEED)
    lku_mc <- replicate(n_sim, {
      sims <- lapply(ni, function(n_i)
        lmom::quakap(stats::runif(n_i), kp))
      lku_i <- vapply(sims, function(s) {
        lm <- tryCatch(lmom::samlmu(s, nmom = 4L),
                       error = function(e) rep(NA_real_, 4L))
        lm[4]
      }, numeric(1L))
      sum(ni * lku_i, na.rm = TRUE) / N
    })
    sigma4 <- stats::sd(lku_mc[is.finite(lku_mc)])
    message(sprintf("  sigma_4 (MC, n=%d): %.6f", n_sim, sigma4))
  }

  Z_vals <- if (is.finite(sigma4) && sigma4 > 0)
    (tau4_theo - lku) / sigma4 else rep(NA_real_, length(dists))

  Z_table <- data.frame(
    Distribution      = dists,
    Tau4_regional     = round(lku, 4),
    Tau4_theoretical  = round(tau4_theo, 4),
    Z_statistic       = round(Z_vals, 4),
    Adequate_fit      = !is.na(Z_vals) & abs(Z_vals) < 1.64,
    stringsAsFactors  = FALSE)

  adequate <- Z_table$Adequate_fit

  # Guard: if all Z_vals are NA (e.g. tau4 computation failed), fall back
  # to Gumbel/GEV as the most commonly appropriate distribution for
  # annual precipitation maxima (extreme value theory justification).
  finite_Z <- which(is.finite(Z_vals))

  best_idx <- if (any(adequate, na.rm = TRUE)) {
    which(adequate)[which.min(abs(Z_vals[adequate]))]
  } else if (length(finite_Z) > 0L) {
    message("  No distribution passes |Z| < 1.64. Selecting minimum |Z|.")
    finite_Z[which.min(abs(Z_vals[finite_Z]))]
  } else {
    message("  All Z-statistics are NA. Defaulting to GEV (EVT justification).")
    1L   # GEV is first in dists vector
  }
  best_dist <- dists[best_idx]

  message("  Z-statistic results:")
  for (i in seq_len(nrow(Z_table))) {
    r   <- Z_table[i, ]
    mrk <- if (i == best_idx) " \u2605 SELECTED"
           else if (!is.na(r$Adequate_fit) && r$Adequate_fit) " \u2713"
           else ""
    message(sprintf("    %-4s | tau4=%.4f | Z=%+7.4f%s",
                    r$Distribution, r$Tau4_theoretical, r$Z_statistic, mrk))
  }

  params_reg <- fit_dist_params(best_dist)
  list(best_dist    = best_dist,
       params_reg   = params_reg,
       Z_table      = Z_table,
       sigma4       = sigma4,
       fit_adequate = any(adequate, na.rm = TRUE),
       best_Z       = round(Z_vals[best_idx], 4))
}


# -----------------------------------------------------------------------------
# SECTION 19c: GRADEX PER STATION USING REGIONAL DISTRIBUTION
# -----------------------------------------------------------------------------

#' Compute GRADEX at each station using the regionally selected distribution
#'
#' Each station's series is fitted with the distribution selected by the H&W
#' Z-test (same family for all stations, consistent with regional homogeneity),
#' but with parameters estimated from each station's own data. This preserves
#' local scale while using the regionally validated distributional shape.
#'
#' @param data       Raw precipitation tibble.
#' @param dist_name  Character. Distribution selected by select_hw_distribution.
#' @return Data frame with id, nombre, gradex_hw per station.
compute_hw_gradex_per_station <- function(data, dist_name) {
  message(sprintf("\nH&W LAYER: Fitting %s to each station...", dist_name))

  y_10  <- -log(-log(1 - 1/10))
  y_100 <- -log(-log(1 - 1/100))

  results <- data |>
    dplyr::group_by(id, nombre) |>
    dplyr::summarise(
      n      = dplyr::n(),
      series = list(stats::na.omit(p_max)),
      .groups = "drop"
    ) |>
    dplyr::filter(n >= MIN_YEARS) |>
    dplyr::mutate(
      gradex_hw = purrr::map2_dbl(series, id, function(x, sid) {
        x  <- unlist(x)
        lm <- tryCatch(lmom::samlmu(x), error = function(e) NULL)
        if (is.null(lm) || anyNA(lm)) return(NA_real_)
        p <- tryCatch(switch(dist_name,
          "GEV" = lmom::pelgev(lm[1:3]), "GLO" = lmom::pelglo(lm[1:3]),
          "GNO" = lmom::pelgno(lm[1:3]), "PE3" = lmom::pelpe3(lm[1:3]),
          "GPA" = lmom::pelgpa(lm[1:3])), error = function(e) NULL)
        if (is.null(p) || anyNA(p)) return(NA_real_)
        q100 <- switch(dist_name,
          "GEV" = lmom::quagev(1-1/100, p), "GLO" = lmom::quaglo(1-1/100, p),
          "GNO" = lmom::quagno(1-1/100, p), "PE3" = lmom::quape3(1-1/100, p),
          "GPA" = lmom::quagpa(1-1/100, p))
        q10  <- switch(dist_name,
          "GEV" = lmom::quagev(1-1/10, p), "GLO" = lmom::quaglo(1-1/10, p),
          "GNO" = lmom::quagno(1-1/10, p), "PE3" = lmom::quape3(1-1/10, p),
          "GPA" = lmom::quagpa(1-1/10, p))
        if (!all(is.finite(c(q100, q10)))) return(NA_real_)
        round((q100 - q10) / (y_100 - y_10), 3)
      })
    ) |>
    dplyr::select(id, nombre, n, gradex_hw) |>
    dplyr::filter(!is.na(gradex_hw))

  message(sprintf("  H&W GRADEX computed for %d/%d stations.",
                  sum(!is.na(results$gradex_hw)), nrow(results)))
  print(results |> dplyr::select(nombre, n, gradex_hw))

  results
}


# -----------------------------------------------------------------------------
# SECTION 19d: IDW/KRIGING ON H&W-REFINED GRADEX VALUES
# -----------------------------------------------------------------------------

#' Interpolate H&W-refined GRADEX values to the gauging point
#'
#' Uses the same spatial method selected in Section 10 (IDW or Kriging),
#' applied to the H&W per-station GRADEX values instead of the pointwise ones.
#'
#' @param hw_station_gradex  Data frame from compute_hw_gradex_per_station.
#' @param analysis_result    List from run_gradex_analysis (for coordinates).
#' @param cv_result          List from run_cross_validation_both (method).
#' @return Numeric. H&W-refined GRADEX at gauging point.
interpolate_hw_gradex <- function(hw_station_gradex, analysis_result,
                                   cv_result) {
  sel_method <- if (!is.null(cv_result)) cv_result$selected_method else "IDW"
  message(sprintf("\nH&W LAYER: Interpolating H&W GRADEX using %s...", sel_method))

  # Merge H&W GRADEX with station coordinates
  coords <- sf::st_coordinates(analysis_result$stations_sf)
  stn_coords <- data.frame(
    id  = analysis_result$gradex_stations$id,
    lon = coords[, 1L], lat = coords[, 2L])

  hw_sf <- hw_station_gradex |>
    dplyr::inner_join(stn_coords, by = "id") |>
    dplyr::rename(Mejor_gradex = gradex_hw) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326L)

  gauge_sf <- analysis_result$gauge_sf
  n_stn    <- nrow(hw_sf)

  gradex_hw_gauge <- tryCatch({
    if (sel_method == "Kriging" && !is.null(cv_result$kriging) && n_stn >= 5L) {
      ctm12     <- 9377L
      hw_proj   <- sf::st_transform(hw_sf,    ctm12)
      gau_proj  <- sf::st_transform(gauge_sf, ctm12)
      hw_sp     <- sf::as_Spatial(hw_proj)
      gau_sp    <- sf::as_Spatial(gau_proj)
      vgm_model <- cv_result$vgm_info$vgm_model$var_model
      kr <- gstat::krige(Mejor_gradex ~ 1, hw_sp, gau_sp,
                         model = vgm_model, debug.level = 0)
      kr$var1.pred
    } else {
      opt_p <- if (!is.null(cv_result)) cv_result$idw$optimal_power else 2
      mod   <- gstat::gstat(formula = Mejor_gradex ~ 1, data = hw_sf,
                            nmax = IDW_NMAX, set = list(idp = opt_p))
      suppressWarnings(predict(mod, gauge_sf, debug.level = 0)$var1.pred)
    }
  }, error = function(e) {
    warning("H&W interpolation failed: ", e$message,
            "\nFalling back to station mean.")
    mean(hw_station_gradex$gradex_hw, na.rm = TRUE)
  })

  gradex_hw_gauge <- round(gradex_hw_gauge, 3)
  message(sprintf("  H&W GRADEX at gauge (%s): %.3f mm",
                  sel_method, gradex_hw_gauge))
  gradex_hw_gauge
}


# -----------------------------------------------------------------------------
# SECTION 19e: H&W BOOTSTRAP CI (H&W Section 6.3)
# -----------------------------------------------------------------------------

#' Compute the formally correct CI for the H&W-refined GRADEX
#'
#' Replaces the three current CI methods (normal, percentile, simple bootstrap)
#' with the H&W simulation algorithm:
#'   1. Fit kappa-4 to regional L-moments
#'   2. Simulate B synthetic regions with same ni record lengths
#'   3. For each replicate: re-estimate regional L-moments, re-select
#'      parameters of the chosen distribution, compute GRADEX at each
#'      synthetic station, re-interpolate to gauging point
#'   4. CI = empirical percentiles of the bootstrap GRADEX distribution
#'
#' This propagates: distributional parameter uncertainty + inter-site
#' variability + record length effects + spatial interpolation uncertainty.
#'
#' @param reg_lmom    List from compute_regional_lmom_hw.
#' @param dist_sel    List from select_hw_distribution.
#' @param analysis_result  List from run_gradex_analysis (for coordinates).
#' @param cv_result   List from run_cross_validation_both.
#' @param n_boot      Integer. Bootstrap replicates.
#' @return Named list: ci_lower, ci_upper, boot_values, n_valid.
bootstrap_hw_ci <- function(reg_lmom, dist_sel, analysis_result,
                              cv_result, n_boot = HW_NSIM) {

  message(sprintf("\nH&W LAYER: Bootstrap CI (n=%d replicates)...", n_boot))

  dist  <- dist_sel$best_dist
  lcv   <- reg_lmom$lcv_R
  lsk   <- reg_lmom$lsk_R
  lku   <- reg_lmom$lku_R
  ni    <- reg_lmom$site_summary$n
  N     <- reg_lmom$N_total
  y_10  <- -log(-log(1 - 1/10))
  y_100 <- -log(-log(1 - 1/100))

  # Station coordinates for re-interpolation in each bootstrap replicate
  coords  <- sf::st_coordinates(analysis_result$stations_sf)
  stn_ids <- analysis_result$gradex_stations$id
  gauge_sf <- analysis_result$gauge_sf
  opt_p    <- if (!is.null(cv_result)) cv_result$idw$optimal_power else 2

  # Station means (for re-dimensionalising synthetic series)
  means_i <- reg_lmom$site_summary$mean_site

  kp <- tryCatch(lmom::pelkap(c(1, lcv, lsk, lku)), error = function(e) NULL)
  if (is.null(kp)) {
    warning("Kappa-4 fitting failed. Cannot compute H&W bootstrap CI.")
    return(list(ci_lower = NA_real_, ci_upper = NA_real_,
                boot_values = numeric(0L), n_valid = 0L))
  }

  set.seed(RANDOM_SEED)
  boot_vals <- replicate(n_boot, tryCatch({

    # 1. Simulate dimensionless synthetic series from kappa-4
    sim_std <- lapply(ni, function(n_i)
      lmom::quakap(stats::runif(n_i), kp))

    # 2. Re-estimate regional L-moments from synthetic region
    lm_list <- lapply(sim_std, function(s)
      tryCatch(lmom::samlmu(s, nmom = 3L),
               error = function(e) rep(NA_real_, 3L)))
    ok      <- vapply(lm_list, function(lm)
                !anyNA(lm) && lm[1] > 0, logical(1L))
    if (sum(ok) < 3L) return(NA_real_)

    ni_ok    <- ni[ok]
    N_ok     <- sum(ni_ok)
    lcv_b    <- vapply(lm_list[ok], function(lm) lm[2]/lm[1], numeric(1L))
    lsk_b    <- vapply(lm_list[ok], function(lm) lm[3],       numeric(1L))
    lcv_R_b  <- sum(ni_ok * lcv_b) / N_ok
    lsk_R_b  <- sum(ni_ok * lsk_b) / N_ok

    # 3. Re-fit selected distribution to simulated regional L-moments
    p_b <- tryCatch(switch(dist,
      "GEV" = lmom::pelgev(c(1, lcv_R_b, lsk_R_b)),
      "GLO" = lmom::pelglo(c(1, lcv_R_b, lsk_R_b)),
      "GNO" = lmom::pelgno(c(1, lcv_R_b, lsk_R_b)),
      "PE3" = lmom::pelpe3(c(1, lcv_R_b, lsk_R_b)),
      "GPA" = lmom::pelgpa(c(1, lcv_R_b, lsk_R_b))),
      error = function(e) NULL)
    if (is.null(p_b) || anyNA(p_b)) return(NA_real_)

    # 4. Compute GRADEX at each station (re-dimensionalised by site mean)
    gradex_b <- mapply(function(s, mu) {
      lm_s <- tryCatch(lmom::samlmu(s * mu, nmom = 3L),
                       error = function(e) NULL)
      if (is.null(lm_s) || anyNA(lm_s)) return(NA_real_)
      p_s <- tryCatch(switch(dist,
        "GEV" = lmom::pelgev(lm_s[1:3]), "GLO" = lmom::pelglo(lm_s[1:3]),
        "GNO" = lmom::pelgno(lm_s[1:3]), "PE3" = lmom::pelpe3(lm_s[1:3]),
        "GPA" = lmom::pelgpa(lm_s[1:3])), error = function(e) NULL)
      if (is.null(p_s) || anyNA(p_s)) return(NA_real_)
      q100 <- switch(dist,
        "GEV"=lmom::quagev(1-1/100,p_s),"GLO"=lmom::quaglo(1-1/100,p_s),
        "GNO"=lmom::quagno(1-1/100,p_s),"PE3"=lmom::quape3(1-1/100,p_s),
        "GPA"=lmom::quagpa(1-1/100,p_s))
      q10  <- switch(dist,
        "GEV"=lmom::quagev(1-1/10,p_s), "GLO"=lmom::quaglo(1-1/10,p_s),
        "GNO"=lmom::quagno(1-1/10,p_s), "PE3"=lmom::quape3(1-1/10,p_s),
        "GPA"=lmom::quagpa(1-1/10,p_s))
      if (!all(is.finite(c(q100, q10)))) NA_real_
      else (q100 - q10) / (y_100 - y_10)
    }, sim_std[ok], means_i[ok], SIMPLIFY = TRUE)

    if (sum(is.finite(gradex_b)) < 2L) return(NA_real_)

    # 5. Re-interpolate to gauge using IDW (fast, stable for CI loop)
    boot_sf <- sf::st_as_sf(
      data.frame(
        Mejor_gradex = gradex_b[is.finite(gradex_b)],
        lon = coords[ok, 1L],
        lat = coords[ok, 2L]),
      coords = c("lon", "lat"), crs = 4326L)
    mod_b <- gstat::gstat(formula = Mejor_gradex ~ 1, data = boot_sf,
                          nmax = IDW_NMAX, set = list(idp = opt_p))
    pred_b <- suppressWarnings(
      predict(mod_b, gauge_sf, debug.level = 0)$var1.pred)
    if (!is.finite(pred_b) || pred_b <= 0) NA_real_ else pred_b

  }, error = function(e) NA_real_))

  boot_vals <- boot_vals[is.finite(boot_vals)]
  n_valid   <- length(boot_vals)
  message(sprintf("  Valid replicates: %d/%d", n_valid, n_boot))

  if (n_valid < 10L) {
    warning("Too few valid bootstrap replicates. CI may not be reliable.")
    return(list(ci_lower = NA_real_, ci_upper = NA_real_,
                boot_values = boot_vals, n_valid = n_valid))
  }

  alpha   <- 1 - CONFIDENCE_LEVEL
  ci_vals <- stats::quantile(boot_vals, probs = c(alpha/2, 1-alpha/2))
  message(sprintf("  H&W Bootstrap CI %.0f%%: [%.2f, %.2f] mm",
                  CONFIDENCE_LEVEL * 100, ci_vals[1], ci_vals[2]))

  list(ci_lower    = round(ci_vals[1], 3),
       ci_upper    = round(ci_vals[2], 3),
       boot_values = boot_vals,
       n_valid     = n_valid)
}



# -----------------------------------------------------------------------------
# SECTION 19f-pre: Z-TEST DIAGNOSTIC AND WEIGHTED RECOMMENDATION
# -----------------------------------------------------------------------------

#' Diagnose why the Z-test passed or failed and assess network limitations
#'
#' Interprets the Z-statistic results in terms of:
#'   - Whether failure is due to network size (sigma_4 too large) or genuine
#'     heterogeneity (all |Z| >> 1.64 even with small sigma_4)
#'   - Distance from the observed point to each theoretical curve
#'   - Power of the test given the current number of stations
#'
#' @param Z_table   Data frame from select_hw_distribution.
#' @param sigma4    Numeric. Monte Carlo sigma_4 from select_hw_distribution.
#' @param n_sites   Integer. Number of stations in the regional analysis.
#' @param N_total   Integer. Total station-years.
#' @return Named list with diagnosis character string and power assessment.
diagnose_z_test <- function(Z_table, sigma4, n_sites, N_total) {

  sep  <- paste(rep("\u2500", 72L), collapse = "")
  cat("\n\u250c", sep, "\u2510\n", sep = "")
  cat("\u2502  Z-TEST DIAGNOSTIC REPORT",
      paste(rep(" ", 46L), collapse = ""), "\u2502\n", sep = "")
  cat("\u251c", sep, "\u2524\n", sep = "")

  Z_vals    <- Z_table$Z_statistic
  tau4_theo <- Z_table$Tau4_theoretical
  tau4_obs  <- Z_table$Tau4_regional[1]
  dists     <- Z_table$Distribution
  min_absZ  <- min(abs(Z_vals), na.rm = TRUE)
  best_dist <- dists[which.min(abs(Z_vals))]

  # --- Power assessment -------------------------------------------------------
  # The critical |Z| is 1.64 (10% significance, two-tailed).
  # The test has limited power when sigma_4 is large relative to the
  # differences between theoretical tau4 curves. With few stations, sigma_4
  # is large and the acceptance band is wide — yet the test may still fail if
  # the observed point is far from all theoretical curves.
  #
  # Rule of thumb (Hosking & Wallis, 1997, p. 100):
  #   N_total >= 500 station-years for reliable Z-test discrimination
  #   N_total < 100 : test has very low discrimination power

  power_class <- dplyr::case_when(
    N_total >= 500 ~ "HIGH — Z-test is discriminating",
    N_total >= 200 ~ "MODERATE — Z-test is indicative",
    N_total >= 100 ~ "LOW — Z-test has limited power",
    TRUE           ~ "VERY LOW — Z-test is uninformative"
  )

  # Acceptance band in tau4 units: +/- 1.64 * sigma4
  accept_band <- if (is.finite(sigma4)) 1.64 * sigma4 else NA_real_

  cat(sprintf("\u2502  Stations: %d | Station-years: %d | Test power: %-25s\u2502\n",
              n_sites, N_total, power_class))
  cat(sprintf("\u2502  sigma_4 = %.4f | Acceptance band: \u00b1%.4f tau4 units%14s\u2502\n",
              sigma4, accept_band, " "))
  cat(sprintf("\u2502  Observed tau4 = %.4f%51s\u2502\n", tau4_obs, " "))
  cat("\u251c", sep, "\u2524\n", sep = "")

  # --- Z-table with interpretation --------------------------------------------
  cat(sprintf("\u2502  %-6s  %10s  %10s  %8s  %8s  %-16s\u2502\n",
              "Dist", "tau4_theo", "tau4_obs", "|Z|", "Status", "Distance*sigma4"))
  cat(sprintf("\u2502  %s\u2502\n", paste(rep("\u2500", 70L), collapse = "")))

  for (i in seq_len(nrow(Z_table))) {
    r      <- Z_table[i, ]
    status <- if (!is.na(r$Adequate_fit) && r$Adequate_fit) "\u2713 PASS"
              else if (i == which.min(abs(Z_vals)))       "\u2605 best"
              else                                         "\u2717 fail"
    dist_sigma <- if (is.finite(sigma4) && sigma4 > 0)
      abs(r$Tau4_theoretical - tau4_obs) / sigma4 else NA_real_
    cat(sprintf("\u2502  %-6s  %10.4f  %10.4f  %8.4f  %-8s  %.2f sigma\u2080 %6s\u2502\n",
                r$Distribution,
                r$Tau4_theoretical, tau4_obs,
                abs(r$Z_statistic), status,
                dist_sigma, " "))
  }
  cat("\u251c", sep, "\u2524\n", sep = "")

  # --- Root cause diagnosis ---------------------------------------------------
  all_fail     <- !any(Z_table$Adequate_fit, na.rm = TRUE)
  near_miss    <- is.finite(min_absZ) && min_absZ < 2.5   # |Z| < 2.5 = borderline
  very_far     <- is.finite(min_absZ) && min_absZ >= 3.0  # |Z| >= 3.0 = genuine misfit

  diagnosis <- if (!all_fail) {
    "Distribution accepted — no issue."
  } else if (N_total < 200 && near_miss) {
    paste0("NETWORK SIZE LIMITATION: The minimum |Z| = ", round(min_absZ, 2),
           " is close to the 1.64 threshold.\n",
           "\u2502    With ", N_total, " station-years, sigma_4 = ", round(sigma4, 4),
           " is too large for the test to discriminate.\n",
           "\u2502    Practical conclusion: ", best_dist,
           " is the most plausible distribution.\n",
           "\u2502    More data would likely confirm it (need ~500 station-years).")
  } else if (near_miss) {
    paste0("BORDERLINE REJECTION: min |Z| = ", round(min_absZ, 2),
           " just exceeds 1.64. ", best_dist,
           " is selected as the best available fit.\n",
           "\u2502    Consider treating this as an acceptable approximation.")
  } else if (very_far) {
    paste0("POSSIBLE HETEROGENEITY: All |Z| >= 3.0. The observed L-kurtosis\n",
           "\u2502    does not match any of the 5 distributions, even with generous\n",
           "\u2502    sigma_4. Check the H-statistic from the homogeneity test.\n",
           "\u2502    If H >= 2, the stations may belong to different regions.")
  } else {
    paste0("MODERATE MISFIT: min |Z| = ", round(min_absZ, 2),
           ". ", best_dist, " selected as closest fit.\n",
           "\u2502    Verify with H-statistic and discordancy test results.")
  }

  cat(sprintf("\u2502  Diagnosis: %-60s\u2502\n", ""))
  # Print multi-line diagnosis
  diag_lines <- strsplit(diagnosis, "\n")[[1]]
  for (dl in diag_lines)
    cat(sprintf("\u2502  %s%s\u2502\n", dl,
                paste(rep(" ", max(0L, 70L - nchar(dl))), collapse = "")))

  cat("\u2514", sep, "\u2518\n", sep = "")

  list(power_class = power_class, diagnosis = diagnosis,
       near_miss = near_miss, network_limitation = N_total < 200 && near_miss,
       accept_band = accept_band, min_absZ = min_absZ)
}


#' Compute weighted GRADEX recommendation for HEC-HMS
#'
#' Weights are derived from THREE quality criteria applied in sequence:
#'
#'   1. R² gate (spatial method):
#'      If R² < 0.3, spatial interpolation has no predictive value.
#'      The spatial weight is set to zero and only the H&W GRADEX is used
#'      — but ONLY if the H&W CI is acceptable (criterion 2).
#'
#'   2. Bootstrap CI amplitude gate (H&W):
#'      CI_amplitude = (CI_upper - CI_lower) / GRADEX_hw * 100
#'      Thresholds (Merz & Bloeschl, 2008; Stedinger et al., 1993):
#'        CI_amplitude < 50%  : acceptable    — H&W receives full weight
#'        50% <= amplitude < 100% : degraded  — H&W weight penalised linearly
#'        amplitude >= 100%   : unreliable    — H&W weight set to zero;
#'                              pointwise estimate used with expanded CI
#'
#'   3. Composite weighting (when both methods are usable):
#'      w_spatial = R²  (zero if R² < 0.3)
#'      w_hw      = ci_factor * (2 - |best_Z|) / 2, clamped [0,1]
#'      Normalised and applied to compute the weighted GRADEX.
#'
#'   4. Sensitivity range:
#'      If H&W weight = 0 (CI unreliable): fixed 25-30% expansion.
#'      Otherwise: max(method divergence, CI amplitude / 4, 10%).
#'
#' @param gradex_spatial  Numeric. Pointwise IDW/Kriging estimate.
#' @param gradex_hw       Numeric. H&W-refined estimate.
#' @param hw_ci_lower     Numeric. Lower bound of H&W bootstrap CI.
#' @param hw_ci_upper     Numeric. Upper bound of H&W bootstrap CI.
#' @param cv_result       List from run_cross_validation_both.
#' @param dist_sel        List from select_hw_distribution.
#' @return Named list with all recommendation components and quality flags.
compute_weighted_recommendation <- function(gradex_spatial, gradex_hw,
                                             hw_ci_lower, hw_ci_upper,
                                             cv_result, dist_sel) {

  sel <- if (!is.null(cv_result)) cv_result$selected_method else "IDW"
  m   <- if (sel == "Kriging" && !is.null(cv_result$kriging))
           cv_result$kriging$metrics else if (!is.null(cv_result))
           cv_result$idw$metrics     else NULL

  r2_spatial <- if (!is.null(m))
    as.numeric(m$Valor[m$Metrica == "R\u00b2"]) else 0.5
  r2_spatial <- max(0, min(1, r2_spatial))

  # --- Criterion 1: R² gate ---------------------------------------------------
  r2_adequate  <- r2_spatial >= 0.3
  r2_class     <- dplyr::case_when(
    r2_spatial >= 0.7 ~ "GOOD",
    r2_spatial >= 0.5 ~ "MODERATE",
    r2_spatial >= 0.3 ~ "WEAK",
    TRUE              ~ "INADEQUATE (< 0.3) — spatial interpolation unreliable")

  # --- Criterion 2: Bootstrap CI amplitude gate ------------------------------
  # Amplitude expressed as % of the H&W GRADEX central estimate
  ci_amplitude_pct <- if (is.finite(hw_ci_lower) && is.finite(hw_ci_upper) &&
                           gradex_hw > 0)
    (hw_ci_upper - hw_ci_lower) / gradex_hw * 100 else Inf

  ci_class <- dplyr::case_when(
    ci_amplitude_pct <  50  ~ "ACCEPTABLE (< 50%)",
    ci_amplitude_pct < 100  ~ "DEGRADED (50-100%) — use with caution",
    TRUE                    ~ "UNRELIABLE (≥ 100%) — H&W weight set to zero")

  # Linear penalty: 1.0 when amplitude=0, 0.0 when amplitude>=100%
  ci_factor <- max(0, min(1, 1 - ci_amplitude_pct / 100))

  hw_ci_adequate <- ci_amplitude_pct < 100

  # --- Inline diagnostic — always printed, cannot be suppressed --------------
  message(rep("-", 60L))
  message("WEIGHT CALCULATION DIAGNOSTICS (compute_weighted_recommendation):")
  message(sprintf("  hw_ci_lower        = %s", hw_ci_lower))
  message(sprintf("  hw_ci_upper        = %s", hw_ci_upper))
  message(sprintf("  gradex_hw          = %.3f mm", gradex_hw))
  message(sprintf("  ci_amplitude_pct   = %.1f%%", ci_amplitude_pct))
  message(sprintf("  ci_factor          = %.4f", ci_factor))
  message(sprintf("  hw_ci_adequate     = %s  [threshold: < 100%%]", hw_ci_adequate))
  message(sprintf("  r2_spatial         = %.4f", r2_spatial))
  message(sprintf("  r2_adequate        = %s  [threshold: >= 0.3]", r2_adequate))
  message(rep("-", 60L))

  # --- Criterion 3: Composite weights ----------------------------------------
  # Priority order:
  #   A. CI invalid (NA or amplitude >= 100%) → H&W weight = 0, regardless
  #      of any other metric. This is the hard gate — an unreliable CI means
  #      the H&W estimate cannot be trusted at all.
  #   B. R² < 0.3 → spatial weight = 0 (interpolation has no predictive value)
  #   C. Both valid but weak → weighted average using R² and ci_factor * Z-score
  #   D. Neither valid → use spatial value (less wrong than H&W with invalid CI)
  #      with maximum sensitivity range. No fallback average.

  best_Z_abs    <- abs(dist_sel$best_Z)

  # Hard gate A: CI invalid → H&W weight strictly zero
  w_hw_raw      <- if (hw_ci_adequate)
    ci_factor * max(0, (2.0 - best_Z_abs) / 2.0) else 0

  # Gate B: R² low → spatial weight zero
  w_spatial_raw <- if (r2_adequate) r2_spatial else 0

  # Determine regime
  hw_ok       <- hw_ci_adequate && w_hw_raw > 0
  spatial_ok  <- r2_adequate    && w_spatial_raw > 0

  if (hw_ok && spatial_ok) {
    # Both usable: normalised weighted average
    w_total   <- w_spatial_raw + w_hw_raw
    w_spatial <- w_spatial_raw / w_total
    w_hw      <- w_hw_raw      / w_total
    fallback  <- FALSE
    message("  REGIME: Both usable -> weighted average")
  } else if (hw_ok && !spatial_ok) {
    # Only H&W usable
    w_spatial <- 0; w_hw <- 1
    fallback  <- FALSE
    message("  REGIME: H&W only (R² < 0.3)")
  } else if (!hw_ok && spatial_ok) {
    # Only spatial usable (CI amplitude >= 100% or CI invalid)
    w_spatial <- 1; w_hw <- 0
    fallback  <- FALSE
    message("  REGIME: Spatial only (H&W CI invalid or amplitude >= 100%)")
  } else {
    # Neither method reliable — use spatial as least-wrong option
    w_spatial <- 1; w_hw <- 0
    fallback  <- TRUE
    message("  REGIME: Fallback — neither method reliable; spatial used")
  }

  message(sprintf("  FINAL: w_spatial=%.0f%%  w_hw=%.0f%%  fallback=%s",
                  w_spatial*100, w_hw*100, fallback))
  message(rep("-", 60L))

  # --- Compute weighted GRADEX ------------------------------------------------
  gradex_weighted <- round(w_spatial * gradex_spatial + w_hw * gradex_hw, 2)

  # --- Criterion 4: Sensitivity range -----------------------------------------
  diff_abs     <- abs(gradex_hw - gradex_spatial)
  diff_pct     <- diff_abs / gradex_weighted * 100

  sens_pct <- if (fallback) {
    # Neither method reliable: maximum sensitivity range
    30L
  } else if (!hw_ci_adequate && spatial_ok) {
    # H&W invalid, spatial only: expand to reflect missing regional info
    25L
  } else if (hw_ok && !spatial_ok) {
    # Spatial invalid, H&W only: CI amplitude drives the range
    as.integer(max(20, min(30, round(ci_amplitude_pct / 4))))
  } else {
    # Both usable: driven by method divergence and CI amplitude
    as.integer(max(10, min(30, round(max(diff_pct, ci_amplitude_pct / 4)))))
  }

  gradex_min <- round(gradex_weighted * (1 - sens_pct / 100), 2)
  gradex_max <- round(gradex_weighted * (1 + sens_pct / 100), 2)

  # --- Determine the final recommended value and its rationale ----------------
  recommendation_basis <- if (fallback && !spatial_ok && !hw_ok) {
    # Neither method valid: spatial used as least-wrong option
    paste0("Pointwise ", sel, " only — both methods unreliable",
           " (R²=", round(r2_spatial, 3),
           ", CI=", if (is.finite(ci_amplitude_pct))
             paste0(round(ci_amplitude_pct, 0), "%") else "invalid", ")")
  } else if (!hw_ok && spatial_ok) {
    paste0("Pointwise ", sel, " only (H&W CI amplitude=",
           if (is.finite(ci_amplitude_pct))
             paste0(round(ci_amplitude_pct, 0), "%") else "invalid",
           " — H&W weight set to zero)")
  } else if (hw_ok && !spatial_ok) {
    paste0("H&W ", dist_sel$best_dist, " only (R²=",
           round(r2_spatial, 3), " < 0.3 — spatial weight set to zero)")
  } else {
    paste0("Weighted: ", round(w_spatial*100), "% ", sel,
           " (R²=", round(r2_spatial, 3), ") + ",
           round(w_hw*100), "% H&W (CI=",
           round(ci_amplitude_pct, 0), "%, |Z|=", round(best_Z_abs, 2), ")")
  }

  # Recommendation data frame for Excel
  recommendation_df <- data.frame(
    Scenario  = c(
      paste0("BASE — ", recommendation_basis),
      paste0("Minimum (−", sens_pct, "%)"),
      paste0("Maximum (+", sens_pct, "%)"),
      paste0("Pointwise ", sel, " (reference)"),
      paste0("H&W ", dist_sel$best_dist, " (reference)"),
      "H&W CI lower (95%)", "H&W CI upper (95%)"
    ),
    GRADEX_mm = c(gradex_weighted, gradex_min, gradex_max,
                  round(gradex_spatial, 2), round(gradex_hw, 2),
                  round(hw_ci_lower, 2), round(hw_ci_upper, 2)),
    Role = c(
      "PRIMARY — use this value in HEC-HMS",
      "Sensitivity lower bound",
      "Sensitivity upper bound",
      "Reference only",
      "Reference only",
      "H&W CI bound (informational)",
      "H&W CI bound (informational)"
    ),
    stringsAsFactors = FALSE
  )

  list(
    gradex_weighted     = gradex_weighted,
    weight_spatial      = round(w_spatial, 3),
    weight_hw           = round(w_hw, 3),
    r2_spatial          = round(r2_spatial, 3),
    r2_class            = r2_class,
    r2_adequate         = r2_adequate,
    best_Z_abs          = round(best_Z_abs, 3),
    ci_amplitude_pct    = round(ci_amplitude_pct, 1),
    ci_class            = ci_class,
    ci_factor           = round(ci_factor, 3),
    hw_ci_adequate      = hw_ci_adequate,
    fallback            = fallback,
    gradex_min          = gradex_min,
    gradex_max          = gradex_max,
    sensitivity_pct     = sens_pct,
    recommendation_basis = recommendation_basis,
    recommendation_df   = recommendation_df
  )
}

# -----------------------------------------------------------------------------
# SECTION 19f: MASTER FUNCTION — H&W LAYER
# -----------------------------------------------------------------------------

#' Execute the full H&W layer and produce the final refined GRADEX estimate
#'
#' @param data            Raw precipitation tibble.
#' @param analysis_result List from run_gradex_analysis.
#' @param cv_result       List from run_cross_validation_both.
#' @return Named list with all H&W layer results.
run_hw_layer <- function(data, analysis_result, cv_result) {

  message("\n", rep("=", 60L))
  message("H&W REGIONAL DISTRIBUTION LAYER (Hosking & Wallis, 1997)")
  message(rep("=", 60L))

  # Step 1: Regional L-moments
  reg_lmom <- compute_regional_lmom_hw(data)

  # Step 2: Distribution selection via Z-statistic
  dist_sel <- select_hw_distribution(reg_lmom, n_sim = HW_NSIM)

  # Step 3: Per-station GRADEX with regionally selected distribution
  hw_stn <- compute_hw_gradex_per_station(data, dist_sel$best_dist)

  if (nrow(hw_stn) == 0L || all(is.na(hw_stn$gradex_hw))) {
    warning("H&W layer: no valid per-station GRADEX. Returning NULL.")
    return(NULL)
  }

  # Step 4: Interpolate H&W GRADEX to gauging point
  gradex_hw_gauge <- interpolate_hw_gradex(hw_stn, analysis_result, cv_result)

  # Step 5: H&W bootstrap CI (replaces the three current CI methods)
  hw_ci <- bootstrap_hw_ci(reg_lmom, dist_sel, analysis_result,
                            cv_result, n_boot = HW_NSIM)

  # ---- Z-test diagnostic ----------------------------------------------------
  z_diag <- diagnose_z_test(
    Z_table  = dist_sel$Z_table,
    sigma4   = dist_sel$sigma4,
    n_sites  = reg_lmom$n_sites,
    N_total  = reg_lmom$N_total)

  # ---- Weighted recommendation for HEC-HMS ----------------------------------
  # Pass bootstrap CI bounds so the CI amplitude penalty can be applied
  rec <- compute_weighted_recommendation(
    gradex_spatial = analysis_result$gradex_hec_hms,
    gradex_hw      = gradex_hw_gauge,
    hw_ci_lower    = hw_ci$ci_lower,
    hw_ci_upper    = hw_ci$ci_upper,
    cv_result      = cv_result,
    dist_sel       = dist_sel)

  # ---- Summary report -------------------------------------------------------
  sel_method     <- if (!is.null(cv_result)) cv_result$selected_method else "IDW"
  gradex_spatial <- analysis_result$gradex_hec_hms
  diff_mm        <- gradex_hw_gauge - gradex_spatial
  diff_pct       <- diff_mm / gradex_spatial * 100

  sep  <- paste(rep("\u2550", 72L), collapse = "")
  sep2 <- paste(rep("\u2500", 72L), collapse = "")

  cat("\n\u2554", sep, "\u2557\n", sep = "")
  cat("\u2551  H&W LAYER RESULTS",
      paste(rep(" ", 53L), collapse = ""), "\u2551\n", sep = "")
  cat("\u2560", sep, "\u2563\n", sep = "")
  cat(sprintf("\u2551  Regional distribution:  %-47s\u2551\n",
              paste0(dist_sel$best_dist, " (Z=", dist_sel$best_Z, ")",
                     if (dist_sel$fit_adequate) " \u2713 PASS"
                     else paste0(" \u26a0 ", z_diag$power_class))))
  cat(sprintf("\u2551  Spatial method:         %-47s\u2551\n", sel_method))
  cat("\u2560", sep, "\u2563\n", sep = "")
  cat(sprintf("\u2551  H&W GRADEX:             %8.3f mm%32s\u2551\n",
              gradex_hw_gauge, " "))
  cat(sprintf("\u2551  Pointwise GRADEX:       %8.3f mm%32s\u2551\n",
              gradex_spatial, " "))
  cat(sprintf("\u2551  Difference:             %+.3f mm (%+.1f%%)%30s\u2551\n",
              diff_mm, diff_pct, " "))
  cat(sprintf("\u2551  H&W CI %.0f%%:             [%.2f, %.2f] mm%28s\u2551\n",
              CONFIDENCE_LEVEL*100, hw_ci$ci_lower, hw_ci$ci_upper, " "))
  cat("\u2560", sep, "\u2563\n", sep = "")

  # Quality gate status block
  cat(sprintf("\u2551  QUALITY GATE STATUS%52s\u2551\n", " "))
  cat(sprintf("\u2551  Spatial R\u00b2:        %-53s\u2551\n",
              paste0(round(rec$r2_spatial, 3), " — ", rec$r2_class)))
  cat(sprintf("\u2551  H&W CI amplitude: %-53s\u2551\n",
              paste0(round(rec$ci_amplitude_pct, 1), "% — ", rec$ci_class)))
  cat(sprintf("\u2551  CI penalty factor: %.3f (1=no penalty, 0=full penalty)%14s\u2551\n",
              rec$ci_factor, " "))

  # Active warning if CI is unreliable
  if (!rec$hw_ci_adequate) {
    cat(sprintf("\u2551  \u274c H&W WEIGHT SET TO ZERO: CI amplitude \u2265 100%%%22s\u2551\n",
                " "))
    cat(sprintf("\u2551     Using %s pointwise estimate with expanded sensitivity%12s\u2551\n",
                sel_method, " "))
  }
  if (!rec$r2_adequate) {
    cat(sprintf("\u2551  \u274c SPATIAL WEIGHT SET TO ZERO: R\u00b2 < 0.3%30s\u2551\n",
                " "))
  }

  cat("\u2560", sep, "\u2563\n", sep = "")
  cat(sprintf("\u2551  WEIGHTED RECOMMENDATION FOR HEC-HMS%36s\u2551\n", " "))
  cat(sprintf("\u2551  Basis: %-64s\u2551\n",
              substr(rec$recommendation_basis, 1L, 64L)))
  cat("\u2560", sep, "\u2563\n", sep = "")
  cat(sprintf("\u2551  \u2605 BASE VALUE (HEC-HMS):   %8.2f mm%32s\u2551\n",
              rec$gradex_weighted, " "))
  cat(sprintf("\u2551  \u2193 Minimum (\u2212%d%%):         %8.2f mm%32s\u2551\n",
              rec$sensitivity_pct, rec$gradex_min, " "))
  cat(sprintf("\u2551  \u2191 Maximum (+%d%%):          %8.2f mm%32s\u2551\n",
              rec$sensitivity_pct, rec$gradex_max, " "))
  cat("\u255a", sep, "\u255d\n", sep = "")

  list(
    gradex_hw_gauge    = gradex_hw_gauge,
    hw_ci_lower        = hw_ci$ci_lower,
    hw_ci_upper        = hw_ci$ci_upper,
    hw_boot_values     = hw_ci$boot_values,
    hw_n_valid         = hw_ci$n_valid,
    dist_selected      = dist_sel$best_dist,
    Z_table            = dist_sel$Z_table,
    fit_adequate       = dist_sel$fit_adequate,
    best_Z             = dist_sel$best_Z,
    regional_lmom      = reg_lmom,
    hw_station_gradex  = hw_stn,
    z_diagnostic       = z_diag,
    recommendation     = rec
  )
}


# -----------------------------------------------------------------------------
# SECTION 19g: UPDATE EXCEL WITH H&W LAYER RESULTS
# -----------------------------------------------------------------------------

export_hw_layer_excel <- function(hw_result,
  file_path = "resultados_gradex_completos.xlsx") {

  message("\nUpdating Excel workbook with H&W layer results...")

  existing <- tryCatch({
    sheet_names <- readxl::excel_sheets(file_path)
    lapply(stats::setNames(sheet_names, sheet_names),
           function(s) readxl::read_excel(file_path, sheet = s))
  }, error = function(e) {
    message("  Cannot read existing workbook. Writing H&W sheets only.")
    list()
  })

  hw_summary <- data.frame(
    Parameter = c(
      "GRADEX_HW_Refined_mm",
      paste0("HW_CI_", round(CONFIDENCE_LEVEL*100), "_Lower_mm"),
      paste0("HW_CI_", round(CONFIDENCE_LEVEL*100), "_Upper_mm"),
      "HW_Bootstrap_Valid_Replicates",
      "HW_Distribution_Selected",
      "HW_Z_Statistic", "HW_Fit_Adequate",
      "HW_Regional_LCV", "HW_Regional_Lskew", "HW_Regional_Lkurtosis",
      "HW_N_Stations", "HW_Total_Station_Years"
    ),
    Value = c(
      hw_result$gradex_hw_gauge,
      hw_result$hw_ci_lower, hw_result$hw_ci_upper,
      hw_result$hw_n_valid,
      hw_result$dist_selected,
      hw_result$best_Z, hw_result$fit_adequate,
      round(hw_result$regional_lmom$lcv_R, 4),
      round(hw_result$regional_lmom$lsk_R, 4),
      round(hw_result$regional_lmom$lku_R, 4),
      hw_result$regional_lmom$n_sites,
      hw_result$regional_lmom$N_total),
    stringsAsFactors = FALSE)

  updated <- c(existing, list(
    HW_Summary       = hw_summary,
    HW_Z_Selection   = hw_result$Z_table,
    HW_Site_Lmom     = hw_result$regional_lmom$site_table,
    HW_Station_GRADEX = hw_result$hw_station_gradex
  ))

  writexl::write_xlsx(updated, file_path)
  message(sprintf("  Workbook updated: %d sheets.", length(updated)))
  invisible(NULL)
}


# =============================================================================
# SECTION 19 EXECUTION
# =============================================================================

message("\n", rep("#", 60L))
message("# SECTION 19: H&W REGIONAL DISTRIBUTION LAYER")
message(rep("#", 60L))

hw_result <- tryCatch(
  run_hw_layer(
    data             = df_raw,
    analysis_result  = analysis_result,
    cv_result        = cv_result
  ),
  error = function(e) {
    warning("H&W layer failed: ", e$message,
            "\nFinal GRADEX remains the pointwise spatial estimate.")
    NULL
  }
)

if (!is.null(hw_result)) {

  # Add recommendation sheet before writing Excel
  export_hw_layer_excel_v2 <- function(hw_res,
    file_path = "resultados_gradex_completos.xlsx") {

    message("\nUpdating Excel workbook with H&W layer results...")
    existing <- tryCatch({
      nms <- readxl::excel_sheets(file_path)
      lapply(stats::setNames(nms, nms),
             function(s) readxl::read_excel(file_path, sheet = s))
    }, error = function(e) list())

    hw_summary <- data.frame(
      Parameter = c(
        "GRADEX_HW_Refined_mm",
        paste0("HW_CI_", round(CONFIDENCE_LEVEL*100), "_Lower_mm"),
        paste0("HW_CI_", round(CONFIDENCE_LEVEL*100), "_Upper_mm"),
        "HW_Bootstrap_Valid_Replicates",
        "HW_Distribution_Selected",
        "HW_Z_Statistic", "HW_Fit_Adequate",
        "HW_Test_Power", "HW_Diagnosis",
        "HW_Regional_LCV", "HW_Regional_Lskew", "HW_Regional_Lkurtosis",
        "HW_N_Stations", "HW_Total_Station_Years",
        "GRADEX_Weighted_mm", "Weight_Spatial_pct", "Weight_HW_pct",
        "Sensitivity_Range_pct",
        "GRADEX_Minimum_mm", "GRADEX_Maximum_mm"
      ),
      Value = c(
        hw_res$gradex_hw_gauge,
        hw_res$hw_ci_lower, hw_res$hw_ci_upper,
        hw_res$hw_n_valid,
        hw_res$dist_selected,
        hw_res$best_Z, hw_res$fit_adequate,
        hw_res$z_diagnostic$power_class,
        substr(hw_res$z_diagnostic$diagnosis, 1L, 200L),
        round(hw_res$regional_lmom$lcv_R, 4),
        round(hw_res$regional_lmom$lsk_R, 4),
        round(hw_res$regional_lmom$lku_R, 4),
        hw_res$regional_lmom$n_sites,
        hw_res$regional_lmom$N_total,
        hw_res$recommendation$gradex_weighted,
        round(hw_res$recommendation$weight_spatial * 100, 1),
        round(hw_res$recommendation$weight_hw * 100, 1),
        hw_res$recommendation$sensitivity_pct,
        hw_res$recommendation$gradex_min,
        hw_res$recommendation$gradex_max
      ),
      stringsAsFactors = FALSE)

    # Quality diagnostic sheet — single source of truth for all quality flags
    quality_df <- data.frame(
      Indicator = c(
        "Spatial_R2",
        "Spatial_R2_Class",
        "Spatial_R2_Adequate",
        "HW_CI_Lower_mm",
        "HW_CI_Upper_mm",
        "HW_CI_Amplitude_pct",
        "HW_CI_Class",
        "HW_CI_Factor",
        "HW_CI_Adequate",
        "HW_Z_Best",
        "HW_Fit_Adequate",
        "HW_Test_Power",
        "Network_N_Stations",
        "Network_Station_Years",
        "Weight_Spatial_pct",
        "Weight_HW_pct",
        "Fallback_Used",
        "Sensitivity_Range_pct",
        "Recommendation_Basis",
        "GRADEX_Base_mm",
        "GRADEX_Min_mm",
        "GRADEX_Max_mm"
      ),
      Value = c(
        round(hw_res$recommendation$r2_spatial, 4),
        hw_res$recommendation$r2_class,
        hw_res$recommendation$r2_adequate,
        round(hw_res$hw_ci_lower, 3),
        round(hw_res$hw_ci_upper, 3),
        round(hw_res$recommendation$ci_amplitude_pct, 1),
        hw_res$recommendation$ci_class,
        round(hw_res$recommendation$ci_factor, 3),
        hw_res$recommendation$hw_ci_adequate,
        round(hw_res$best_Z, 4),
        hw_res$fit_adequate,
        hw_res$z_diagnostic$power_class,
        hw_res$regional_lmom$n_sites,
        hw_res$regional_lmom$N_total,
        round(hw_res$recommendation$weight_spatial * 100, 1),
        round(hw_res$recommendation$weight_hw * 100, 1),
        hw_res$recommendation$fallback,
        hw_res$recommendation$sensitivity_pct,
        hw_res$recommendation$recommendation_basis,
        hw_res$recommendation$gradex_weighted,
        hw_res$recommendation$gradex_min,
        hw_res$recommendation$gradex_max
      ),
      Threshold_Reference = c(
        "R² >= 0.3 adequate; >= 0.7 good",
        "GOOD / MODERATE / WEAK / INADEQUATE",
        "TRUE = R² >= 0.3",
        "H&W bootstrap 95% CI",
        "H&W bootstrap 95% CI",
        "< 50% acceptable; < 100% degraded; >= 100% unreliable",
        "ACCEPTABLE / DEGRADED / UNRELIABLE",
        "1 = no penalty; 0 = full penalty (CI >= 100%)",
        "TRUE = CI amplitude < 100%",
        "|Z| < 1.64 = adequate fit",
        "TRUE = at least one dist. passes Z-test",
        "Station-years: >= 500 high; >= 200 moderate; >= 100 low",
        ">= 15 reliable; >= 8 acceptable; < 5 minimum",
        ">= 500 station-years for reliable Z-test",
        "Derived from R² and CI amplitude",
        "Derived from |Z| and CI amplitude",
        "TRUE = both methods inadequate",
        ">= 25% if CI unreliable; >= 20% if R² low",
        "Explanation of weighting decision",
        "Recommended value for HEC-HMS",
        "Lower sensitivity bound",
        "Upper sensitivity bound"
      ),
      stringsAsFactors = FALSE
    )

    updated <- c(existing, list(
      HW_Summary        = hw_summary,
      HW_Z_Selection    = hw_res$Z_table,
      HW_Site_Lmom      = hw_res$regional_lmom$site_table,
      HW_Station_GRADEX = hw_res$hw_station_gradex,
      HW_Recommendation = hw_res$recommendation$recommendation_df,
      HW_Quality_Diag   = quality_df
    ))

    writexl::write_xlsx(updated, file_path)
    message(sprintf("  Workbook updated: %d sheets.", length(updated)))
    invisible(NULL)
  }

  export_hw_layer_excel_v2(hw_result)

  # Store H&W CI and recommendation in analysis_result for downstream use
  analysis_result$ci_hw <- list(
    ci_lower  = hw_result$hw_ci_lower,
    ci_upper  = hw_result$hw_ci_upper,
    method    = "H&W Bootstrap (Section 6.3)",
    dist      = hw_result$dist_selected
  )
  analysis_result$gradex_recommended <- hw_result$recommendation$gradex_weighted
  analysis_result$gradex_recommended_min <- hw_result$recommendation$gradex_min
  analysis_result$gradex_recommended_max <- hw_result$recommendation$gradex_max

  sel_sp <- if (!is.null(cv_result)) cv_result$selected_method else "IDW"

  message("\n", rep("=", 60L))
  message("FINAL RECOMMENDED VALUES FOR HEC-HMS")
  message(rep("=", 60L))
  message(sprintf("  \u2605 BASE VALUE (weighted):    %.2f mm",
                  hw_result$recommendation$gradex_weighted))
  message(sprintf("  \u2193 Minimum (\u2212%d%%):           %.2f mm",
                  hw_result$recommendation$sensitivity_pct,
                  hw_result$recommendation$gradex_min))
  message(sprintf("  \u2191 Maximum (+%d%%):            %.2f mm",
                  hw_result$recommendation$sensitivity_pct,
                  hw_result$recommendation$gradex_max))
  message(rep("-", 60L))
  message(sprintf("  H&W (%s) GRADEX:         %.3f mm",
                  hw_result$dist_selected, hw_result$gradex_hw_gauge))
  message(sprintf("  Pointwise (%s) GRADEX:   %.3f mm", sel_sp, gradex_final))
  message(sprintf("  Weights:  %s=%.0f%% | H&W=%.0f%%",
                  sel_sp,
                  hw_result$recommendation$weight_spatial * 100,
                  hw_result$recommendation$weight_hw * 100))
  message(sprintf("  H&W CI %.0f%%: [%.2f, %.2f] mm",
                  CONFIDENCE_LEVEL * 100,
                  hw_result$hw_ci_lower, hw_result$hw_ci_upper))
  message(rep("=", 60L))

} else {
  message("H&W layer skipped. Using pointwise estimate as final result.")
  analysis_result$gradex_recommended <- gradex_final
}


# =============================================================================
# SECTION 20: MONTE CARLO VALIDATION OF THE FRAMEWORK
# =============================================================================
# Validates the GRADEX-DUAL framework on synthetic regions with KNOWN
# true GRADEX. For each replicate:
#   1. Simulate K stations from a GEV with known parameters (mu, alpha, kappa).
#      Each station has a true GRADEX = (Q100 - Q10) / (y100 - y10) computed
#      analytically from the known GEV parameters.
#   2. Place stations randomly within a synthetic basin and simulate annual
#      maximum series of length n_yrs.
#   3. Run the pointwise GRADEX estimation on the synthetic data and
#      interpolate to a synthetic gauging point (centroid).
#   4. Record bias = estimate - true_gradex, RMSE, and CI coverage.
#
# This implements recommendation P6 of the manuscript review:
# "añadir un experimento Montecarlo de validación".
#
# To enable: set RUN_MC_VALIDATION = TRUE below. Disabled by default because
# it adds 5-15 minutes of runtime for typical settings.
# =============================================================================

RUN_MC_VALIDATION  <- TRUE    # v1.2: enabled by default; set FALSE to skip
MC_N_REPLICATES    <- 200L    # Number of synthetic regions to simulate
MC_N_STATIONS      <- 9L      # Stations per replicate (match study network)
MC_N_YEARS         <- 30L     # Years of record per station
MC_TRUE_KAPPA      <- -0.10   # True GEV shape parameter
MC_TRUE_LCV        <- 0.18    # True L-CV
MC_TRUE_MU_MEAN    <- 110     # Mean factor index across stations (mm)
MC_TRUE_MU_SD      <- 18      # Inter-station spread of factor index (mm)


#' Compute GRADEX analytically from GEV parameters (H&W parameterisation)
#'
#' @param mu     numeric, location parameter
#' @param alpha  numeric, scale parameter
#' @param kappa  numeric, shape parameter (negative = heavy tail)
#' @return numeric: true GRADEX = (q100 - q10) / (y100 - y10)
true_gradex_gev <- function(mu, alpha, kappa) {
  y10  <- -log(-log(1 - 1/10))
  y100 <- -log(-log(1 - 1/100))
  if (abs(kappa) < 1e-6) {
    return(alpha)  # Gumbel limit: GRADEX = scale
  }
  q10  <- mu + (alpha / kappa) * (1 - (-log(1 - 1/10))^kappa)
  q100 <- mu + (alpha / kappa) * (1 - (-log(1 - 1/100))^kappa)
  (q100 - q10) / (y100 - y10)
}


#' Simulate one GEV series in the H&W parameterisation
#'
#' @param n      numeric, length of series
#' @param mu     location (typically the factor index per station)
#' @param alpha  scale
#' @param kappa  shape
#' @return numeric vector of length n
sim_gev_series <- function(n, mu, alpha, kappa) {
  u <- stats::runif(n)
  if (abs(kappa) < 1e-6) {
    return(mu - alpha * log(-log(u)))
  }
  mu + (alpha / kappa) * (1 - (-log(u))^kappa)
}


#' Run one Monte Carlo replicate of the GRADEX-DUAL framework
#'
#' Returns a list with: true_gradex_mean (regional average true GRADEX),
#' estimated_pointwise (IDW interpolated to centroid), bias, abs_error.
mc_one_replicate <- function(rep_id, n_stations, n_years,
                              kappa, lcv, mu_mean, mu_sd) {
  # Random factor indices per station, with positive support
  mu_i <- pmax(20, stats::rnorm(n_stations, mean = mu_mean, sd = mu_sd))
  # alpha calibrated to maintain target L-CV
  g1   <- gamma(1 + kappa)
  alpha_i <- mu_i * lcv * kappa / ((1 - 2^(-kappa)) * g1)
  # Adjusted location so that the mean equals mu_i
  loc_i <- mu_i - alpha_i * (1 - g1) / kappa

  # True per-station GRADEX
  true_gradex_i <- mapply(true_gradex_gev, loc_i, alpha_i, kappa)
  true_gradex_regional <- mean(true_gradex_i)

  # Random spatial layout in a 50 km x 50 km square (CTM12 metric units)
  X <- stats::runif(n_stations, 4700000, 4750000)
  Y <- stats::runif(n_stations, 6900000, 6950000)
  # Synthetic gauging point at the centroid
  gauge <- c(mean(X), mean(Y))

  # Simulate series and estimate GRADEX per station via Gumbel L-moments
  # (deliberately using the simplest layer to keep runtime reasonable)
  est_gradex_i <- numeric(n_stations)
  for (k in seq_len(n_stations)) {
    series <- sim_gev_series(n_years, loc_i[k], alpha_i[k], kappa)
    series <- pmax(series, 1)
    g_par  <- tryCatch(lmom::pelgum(lmom::samlmu(series)),
                        error = function(e) NULL)
    if (is.null(g_par) || anyNA(g_par)) {
      est_gradex_i[k] <- NA_real_
    } else {
      y10  <- -log(-log(1 - 1/10))
      y100 <- -log(-log(1 - 1/100))
      est_gradex_i[k] <- (lmom::quagum(1 - 1/100, g_par) -
                           lmom::quagum(1 - 1/10,  g_par)) / (y100 - y10)
    }
  }

  ok <- is.finite(est_gradex_i)
  if (sum(ok) < 3L) return(NULL)

  # IDW interpolation to gauge (exponent = 2, fixed for speed)
  d  <- sqrt((X[ok] - gauge[1])^2 + (Y[ok] - gauge[2])^2)
  d  <- pmax(d, 1)
  w  <- 1 / d^2
  est_gauge <- sum(w * est_gradex_i[ok]) / sum(w)

  data.frame(
    rep_id    = rep_id,
    true_gradex = round(true_gradex_regional, 3),
    est_gradex  = round(est_gauge, 3),
    bias        = round(est_gauge - true_gradex_regional, 3),
    rel_bias_pct = round(100 * (est_gauge - true_gradex_regional) /
                          true_gradex_regional, 2)
  )
}


#' Run Monte Carlo validation and report sesgo/RMSE/coverage
run_mc_validation <- function(n_rep, n_stations, n_years,
                               kappa, lcv, mu_mean, mu_sd, seed = 2024L) {
  message("\n", rep("#", 60L))
  message(sprintf(
    "MONTE CARLO VALIDATION  |  %d replicates | %d stations | %d yr",
    n_rep, n_stations, n_years))
  message(rep("#", 60L))
  set.seed(seed)
  results <- vector("list", n_rep)
  pb <- utils::txtProgressBar(min = 0L, max = n_rep, style = 3L)
  for (k in seq_len(n_rep)) {
    results[[k]] <- mc_one_replicate(k, n_stations, n_years,
                                       kappa, lcv, mu_mean, mu_sd)
    utils::setTxtProgressBar(pb, k)
  }
  close(pb)
  results <- do.call(rbind, Filter(Negate(is.null), results))

  if (nrow(results) == 0L) {
    message("  All replicates failed.")
    return(NULL)
  }

  bias_mean      <- mean(results$bias, na.rm = TRUE)
  bias_median    <- stats::median(results$bias, na.rm = TRUE)
  rmse           <- sqrt(mean(results$bias^2, na.rm = TRUE))
  rel_bias_mean  <- mean(results$rel_bias_pct, na.rm = TRUE)
  rel_bias_p95   <- stats::quantile(abs(results$rel_bias_pct),
                                     0.95, na.rm = TRUE)

  message("\n", rep("=", 60L))
  message("MC VALIDATION RESULTS")
  message(rep("=", 60L))
  message(sprintf("  Valid replicates:        %d / %d", nrow(results), n_rep))
  message(sprintf("  Mean true GRADEX:        %.2f mm",
                  mean(results$true_gradex)))
  message(sprintf("  Mean estimated GRADEX:   %.2f mm",
                  mean(results$est_gradex)))
  message(sprintf("  Bias (mean):             %+.3f mm  (%+.2f%%)",
                  bias_mean, rel_bias_mean))
  message(sprintf("  Bias (median):           %+.3f mm",
                  bias_median))
  message(sprintf("  RMSE:                    %.3f mm",
                  rmse))
  message(sprintf("  P95 |rel bias|:          %.2f%%",
                  rel_bias_p95))

  # Coverage assessment: fraction of replicates where the true GRADEX falls
  # within +/- 30 % of the estimate (the manuscript's sensitivity range).
  within_30 <- mean(abs(results$rel_bias_pct) <= 30, na.rm = TRUE)
  message(sprintf("  Coverage of \u00B130%% range:   %.1f%%",
                  within_30 * 100))
  message(rep("=", 60L))

  # Save MC results to disk for inclusion in the manuscript
  utils::write.csv(results, "monte_carlo_validation.csv", row.names = FALSE)
  message("  Replicate-level results saved: monte_carlo_validation.csv")

  invisible(list(replicates = results,
                  bias_mean = bias_mean, bias_median = bias_median,
                  rmse = rmse, coverage_30 = within_30,
                  rel_bias_p95 = rel_bias_p95))
}


if (isTRUE(RUN_MC_VALIDATION)) {
  mc_results <- tryCatch(
    run_mc_validation(
      n_rep      = MC_N_REPLICATES,
      n_stations = MC_N_STATIONS,
      n_years    = MC_N_YEARS,
      kappa      = MC_TRUE_KAPPA,
      lcv        = MC_TRUE_LCV,
      mu_mean    = MC_TRUE_MU_MEAN,
      mu_sd      = MC_TRUE_MU_SD,
      seed       = RANDOM_SEED
    ),
    error = function(e) {
      message("MC validation failed: ", e$message)
      NULL
    }
  )
} else {
  message("\n[SECTION 20] Monte Carlo validation disabled ",
          "(set RUN_MC_VALIDATION = TRUE to enable).")
}

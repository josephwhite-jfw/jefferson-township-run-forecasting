# ─────────────────────────────────────────────────────────────────────────────
# Jefferson Township – Monthly run models (R)
#   - Loads panel via here::here()
#   - Poisson & NegBin with 6‑month lags (sqft per 1,000; apprbld per $1M)
#   - Month fixed effects + nursing home beds control
#   - YoY "change model" (OLS on deltas)
#   - Marginal effects at the mean + simple scenario calculator
# ─────────────────────────────────────────────────────────────────────────────

# 0) Packages -----------------------------------------------------------------
req <- c("here","tidyverse","lubridate","janitor","broom","MASS","modelsummary")
inst <- rownames(installed.packages())
if (any(!req %in% inst)) install.packages(setdiff(req, inst))
invisible(lapply(req, library, character.only = TRUE))

# 1) Paths & read -------------------------------------------------------------
PANEL_PATH <- here::here("data","clean","panel_monthly_with_parcels.csv")
panel <- readr::read_csv(PANEL_PATH, show_col_types = FALSE) |> janitor::clean_names()

# 2) Prep variables -----------------------------------------------------------
# Expected engineered columns from your panel builder:
#  - res_area_sqft_k_lag6, com_area_sqft_k_lag6, ind_area_sqft_k_lag6
#  - res_apprbld_m_lag6,   com_apprbld_m_lag6,   ind_apprbld_m_lag6
#  - nh_total_certified_beds
#  - calls_yoy, *_k_yoy, *_m_yoy (for change model)
stopifnot("total_calls" %in% names(panel))

panel <- panel |>
  mutate(
    month = as.Date(month),
    month_num = factor(lubridate::month(month)),          # seasonality
    nh_beds_100 = nh_total_certified_beds / 100
  )

# Keep a modeling df for the lagged-level models
lag_vars <- c("res_area_sqft_k_lag6","com_area_sqft_k_lag6","ind_area_sqft_k_lag6",
              "res_apprbld_m_lag6","com_apprbld_m_lag6","ind_apprbld_m_lag6")
need_cols <- c("total_calls","nh_beds_100","month_num", lag_vars)
missing <- setdiff(need_cols, names(panel))
if (length(missing)) {
  stop(glue::glue("Missing expected columns in panel: {paste(missing, collapse=', ')}"))
}

mod_df <- panel |> dplyr::select(all_of(need_cols)) |> drop_na()

# 3) Poisson & Negative Binomial (lagged levels) ------------------------------
form <- as.formula(
  "total_calls ~ nh_beds_100 + res_area_sqft_k_lag6 + com_area_sqft_k_lag6 + ind_area_sqft_k_lag6 +
                 res_apprbld_m_lag6 + com_apprbld_m_lag6 + ind_apprbld_m_lag6 +
                 month_num"
)

pois_fit <- glm(form, data = mod_df, family = poisson(link = "log"))
nb_fit   <- MASS::glm.nb(form, data = mod_df)  # log link by default

# Marginal effects at the mean: dE[y]/dx = beta * mean(y) for log-link
mu <- mean(mod_df$total_calls, na.rm = TRUE)
me_table <- function(fit, label_map) {
  coef <- broom::tidy(fit) |> filter(term != "(Intercept)", !str_starts(term,"month_num"))
  coef |>
    mutate(
      me_month = estimate * mu,
      me_year  = me_month * 12
    ) |>
    transmute(
      variable = dplyr::recode(term, !!!label_map),
      `ΔRuns/Month (at mean)` = round(me_month, 3),
      `ΔRuns/Year (×12)`      = round(me_year, 3),
      `p.value`               = round(p.value, 4),
      `std.err`               = round(std.error, 6)
    )
}

labels <- c(
  nh_beds_100 = "Nursing home beds (+100)",
  res_area_sqft_k_lag6 = "Residential sqft (+1,000), lag 6m",
  com_area_sqft_k_lag6 = "Commercial sqft (+1,000), lag 6m",
  ind_area_sqft_k_lag6 = "Industrial sqft (+1,000), lag 6m",
  res_apprbld_m_lag6   = "Residential appraised (+$1M), lag 6m",
  com_apprbld_m_lag6   = "Commercial appraised (+$1M), lag 6m",
  ind_apprbld_m_lag6   = "Industrial appraised (+$1M), lag 6m"
)

pois_me <- me_table(pois_fit, labels)
nb_me   <- me_table(nb_fit,   labels)

cat("\n── Poisson marginal effects at mean (monthly & annualized) ──\n")
print(pois_me)
cat("\n── NegBin marginal effects at mean (monthly & annualized) ──\n")
print(nb_me)

# 4) CHANGE model (YoY deltas) -----------------------------------------------
# Focus on growth: runs change vs. YoY changes in sqft & appraised & beds
change_vars <- c("calls_yoy",
                 "res_area_sqft_k_yoy","com_area_sqft_k_yoy","ind_area_sqft_k_yoy",
                 "res_apprbld_m_yoy","com_apprbld_m_yoy","ind_apprbld_m_yoy",
                 "nh_beds_yoy")

missing_change <- setdiff(change_vars, names(panel))
if (length(missing_change)) {
  warning(paste("Some YoY columns missing; change model may be skipped:", 
                paste(missing_change, collapse=", ")))
}

chg_df <- panel |> dplyr::select(any_of(change_vars)) |> drop_na()
if (nrow(chg_df) >= 20) {
  chg_fit <- lm(calls_yoy ~ res_area_sqft_k_yoy + com_area_sqft_k_yoy + ind_area_sqft_k_yoy +
                  res_apprbld_m_yoy + com_apprbld_m_yoy + ind_apprbld_m_yoy +
                  nh_beds_yoy, data = chg_df)
  
  chg_tab <- broom::tidy(chg_fit) |>
    filter(term != "(Intercept)") |>
    transmute(
      predictor = dplyr::recode(term,
                                res_area_sqft_k_yoy = "Residential Δsqft (+1,000) YoY",
                                com_area_sqft_k_yoy = "Commercial Δsqft (+1,000) YoY",
                                ind_area_sqft_k_yoy = "Industrial Δsqft (+1,000) YoY",
                                res_apprbld_m_yoy   = "Residential Δappraised (+$1M) YoY",
                                com_apprbld_m_yoy   = "Commercial Δappraised (+$1M) YoY",
                                ind_apprbld_m_yoy   = "Industrial Δappraised (+$1M) YoY",
                                nh_beds_yoy         = "Nursing beds Δ(+1) YoY"
      ),
      `ΔRuns/Month (from change)` = round(estimate, 3),
      `ΔRuns/Year (×12)`         = round(estimate * 12, 3),
      `p.value`                  = round(p.value, 4),
      `std.err`                  = round(std.error, 6)
    )
  
  cat("\n── Change model (OLS on YoY deltas): monthly & annualized effects ──\n")
  print(chg_tab)
} else {
  chg_fit <- NULL
  chg_tab <- tibble(note = "Not enough non-missing YoY rows to run change model.")
  print(chg_tab)
}

# 5) Correlations among development predictors (diagnostics) ------------------
dev_cols <- c("res_area_sqft_k_lag6","com_area_sqft_k_lag6","ind_area_sqft_k_lag6",
              "res_apprbld_m_lag6","com_apprbld_m_lag6","ind_apprbld_m_lag6")
corr_tbl <- mod_df |> dplyr::select(all_of(dev_cols)) |> cor(use = "pairwise.complete.obs")
cat("\n── Correlation among lagged development predictors ──\n")
print(round(corr_tbl, 3))

# 6) Scenario calculator (simple helpers) -------------------------------------
# A) Using CHANGE model (recommended for planning if available)
scenario_change <- function(delta_com_sqft_k = 0, delta_res_sqft_k = 0, delta_ind_sqft_k = 0,
                            delta_com_app_m = 0, delta_res_app_m = 0, delta_ind_app_m = 0,
                            delta_nh_beds = 0) {
  if (is.null(chg_fit)) stop("Change model not available (insufficient YoY rows).")
  newx <- tibble(
    com_area_sqft_k_yoy = delta_com_sqft_k,
    res_area_sqft_k_yoy = delta_res_sqft_k,
    ind_area_sqft_k_yoy = delta_ind_sqft_k,
    com_apprbld_m_yoy   = delta_com_app_m,
    res_apprbld_m_yoy   = delta_res_app_m,
    ind_apprbld_m_yoy   = delta_ind_app_m,
    nh_beds_yoy         = delta_nh_beds
  )
  dm <- model.matrix(formula(chg_fit), newx)
  # remove intercept col if present
  if ("(Intercept)" %in% colnames(dm)) dm <- dm[, setdiff(colnames(dm), "(Intercept)"), drop = FALSE]
  monthly <- as.numeric(dm %*% coef(chg_fit)[colnames(dm)])
  list(
    monthly_change = monthly,
    yearly_change  = monthly * 12
  )
}

# B) Using Poisson at-mean marginal effects (levels/lag model)
#    Provide deltas in same units as predictors: sqft in 1,000; appr in $1M; beds in 100s.
scenario_poisson <- function(delta_res_sqft_k_lag6 = 0, delta_com_sqft_k_lag6 = 0, delta_ind_sqft_k_lag6 = 0,
                             delta_res_app_m_lag6 = 0,  delta_com_app_m_lag6 = 0,  delta_ind_app_m_lag6 = 0,
                             delta_nh_beds_100 = 0) {
  b <- coef(pois_fit)
  # map variables to betas (ignore month FE)
  pick <- function(name) ifelse(name %in% names(b), b[[name]], 0)
  me_month <- ( pick("res_area_sqft_k_lag6") * delta_res_sqft_k_lag6 +
                  pick("com_area_sqft_k_lag6") * delta_com_sqft_k_lag6 +
                  pick("ind_area_sqft_k_lag6") * delta_ind_sqft_k_lag6 +
                  pick("res_apprbld_m_lag6")   * delta_res_app_m_lag6 +
                  pick("com_apprbld_m_lag6")   * delta_com_app_m_lag6 +
                  pick("ind_apprbld_m_lag6")   * delta_ind_app_m_lag6 +
                  pick("nh_beds_100")          * delta_nh_beds_100 ) * mu
  list(monthly_change = me_month, yearly_change = me_month * 12)
}

# 7) Example scenarios (comment/uncomment) ------------------------------------
# cat("\nExample scenario (CHANGE model): +80k commercial sqft YoY, +10 NH beds YoY\n")
# print(scenario_change(delta_com_sqft_k = 80, delta_nh_beds = 10))
#
# cat("\nExample scenario (Poisson at-mean): +50k com sqft (lagged), +$5M res appraised (lagged)\n")
# print(scenario_poisson(delta_com_sqft_k_lag6 = 50, delta_res_app_m_lag6 = 5))

# 8) Nicely formatted model tables (optional) ---------------------------------
modelsummary::modelsummary(
  list("Poisson (lagged levels)" = pois_fit, "NegBin (lagged levels)" = nb_fit, "Change model (OLS YoY)" = chg_fit),
  gof_omit = "IC|Log|Adj|Pseudo",
  output = "markdown"
)

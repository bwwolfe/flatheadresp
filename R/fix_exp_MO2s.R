#' Fix \eqn{MO_{2}}  values for experimental chambers based on corrected fish masses and densities
#'
#' @description
#' Adjusts mass-specific \eqn{MO_{2}} values in experiment data by recalculating based on
#' corrected fish masses and/or densities. Can operate on data loaded from an
#' experiment directory (`path`) or on user-supplied data frames.
#'
#' @details
#' If `path` is provided, metadata and \eqn{MO_{2}}  data are automatically loaded using
#' `get_exp_metadata()` and either `calc_exp_mo2s()` or `get_exp_mo2s()` depending
#' on `uncorrected`. If `path` is `NULL`, both `mo2_data` and `metadata` must be
#' supplied. Chambers can be selectively updated by providing `NA` for unchanged
#' masses or densities.
#'
#' @param path Character string. Path to experiment directory. If provided, overrides
#' `mo2_data` and `metadata`.
#' @param mo2_data Data frame of \eqn{MO_{2}}  values (from `calc_exp_mo2s()` or `get_exp_mo2s()`).
#' Required if `path` is `NULL`.
#' @param metadata Data frame of chamber metadata. Required if `path` is `NULL`.
#' @param new_masses Numeric vector or list of new fish masses (kg). Can be named
#' by chamber ID (e.g., `list("1" = 0.345)`) or unnamed (order corresponds to chambers).
#' Use `NA` for chambers that should not be changed.
#' @param new_densities Numeric scalar or vector of fish densities (kg/L). Same naming rules
#' as `new_masses`. Defaults to `1.0`.
#' @param mo2_col Column name for \eqn{MO_{2}}  values in `mo2_data`. Default `"MO2"`.
#' @param chamber_col Column name for chamber IDs. Default `"chamber"`.
#' @param mass_col Column name in metadata for original fish mass (kg).
#' @param vresp_col Column name in metadata for respirometer volume (L).
#' @param keep_original Logical. If `TRUE`, preserves original \eqn{MO_{2}} values in a new
#' column (`MO2_orig`). Default `TRUE`.
#' @param uncorrected Logical. If `TRUE`, loads \eqn{MO_{2}}  data using `get_exp_mo2s()`
#' instead of `calc_exp_mo2s()`. Default `FALSE`.
#' @param calc_report Logical. If `TRUE`, will print the `calc_cycle_mo2s()` summary to
#' console. Default `FALSE`.
#'
#' @return A data frame identical to `calc_exp_mo2s()` output but with corrected \eqn{MO_{2}}  values.
#'
#' @examples
#' # Example using system file from package
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment", package = "flatheadresp")
#' corrected <- fix_exp_mo2s(path = exp_dir_path,
#'                           new_masses = c("1" = 0.345, "3" = 0.512),
#'                           new_densities = 1.02)
#'
#' @export

fix_exp_mo2s <- function(path = NULL,
                         mo2_data = NULL,
                         metadata = NULL,
                         new_masses = NULL,
                         new_densities = 1.0,
                         mo2_col = "MO2",
                         chamber_col = "chamber",
                         mass_col = "Mass.of.fish..kg",
                         vresp_col = "Volume.respirometer..L",
                         keep_original = TRUE,
                         uncorrected = FALSE,
                         calc_report = FALSE) {

  ## --- Input handling ---
  if (!is.null(path)) {
    if (!is.null(mo2_data) || !is.null(metadata)) {
      cli::cli_warn("'path' is provided, so provided 'mo2_data' and/or `metadata` will be overwritten. Set 'path = NULL' to provide these arguments manually.")
    }
    metadata <- get_exp_metadata(path)
    mo2_data <- if (isTRUE(uncorrected)) {
      get_exp_mo2s(path)
    } else {
     calc_exp_mo2s(path, report = calc_report)
    }


  } else {
    if (is.null(mo2_data) || !is.data.frame(mo2_data)) {
      stop("When 'path' is NULL, 'mo2_data' must be provided as a dataframe.")
    }
    if (is.null(metadata) || !is.data.frame(metadata)) {
      stop("When 'path' is NULL, 'metadata' must be provided as a dataframe.")
    }
  }

  ## --- Chamber consistency check ---
  chambers_meta <- sort(unique(metadata[[chamber_col]]))
  chambers_data <- sort(unique(mo2_data[[chamber_col]]))
  if (length(chambers_meta) != length(chambers_data)) {
    cli::cli_warn("Number of chambers differ: {length(chambers_meta)} in metadata vs {length(chambers_data)} in mo2_data.")
  }

  if(new_densities == 1.0 && is.null(new_masses)) stop ("No new_densities or new_masses provided. Provide at least one of these or there is nothing to fix.")

  ## --- Normalize new_masses ---
  n_chambers <- length(chambers_meta)
  if (is.null(new_masses)) new_masses <- rep(NA, n_chambers)
  if (is.null(names(new_masses))) {
    if (length(new_masses) != n_chambers) {
      stop("Unnamed new_masses must have length equal to number of chambers in metadata.")
    }
    names(new_masses) <- as.character(chambers_meta)
  } else {
    all_names <- as.character(chambers_meta)
    nm_filled <- rep(NA_real_, length(all_names))
    names(nm_filled) <- all_names
    nm_filled[names(new_masses)] <- unlist(new_masses)
    new_masses <- nm_filled
  }


  ## --- Normalize new_densities ---
  if (length(new_densities) == 1) {
    # Single value: replicate for all chambers
    dens_map <- stats::setNames(rep(new_densities, n_chambers), as.character(chambers_meta))
  } else if (is.null(names(new_densities))) {
    # Unnamed vector: must match chamber count
    if (length(new_densities) != n_chambers) {
      stop("Unnamed new_densities must match number of chambers.")
    }
    dens_map <- stats::setNames(new_densities, as.character(chambers_meta))
  } else {
    # Named vector or list: fill missing chambers with NA first
    all_names <- as.character(chambers_meta)
    dens_filled <- rep(NA_real_, length(all_names))
    names(dens_filled) <- all_names
    dens_filled[names(new_densities)] <- unlist(new_densities)
    dens_map <- dens_filled
  }

  dens_map[is.na(dens_map)] <- 1.0

  ## --- Check required columns ---
  req_cols <- c(mo2_col, chamber_col)
  if (!all(req_cols %in% names(mo2_data))) {
    stop("Missing required columns in mo2_data: ",
         paste(setdiff(req_cols, names(mo2_data)), collapse = ", "))
  }

  ## --- Preserve original MO2if requested ---
  mo2_orig_col <- paste0(mo2_col, "_orig")
  if (isTRUE(keep_original) && !mo2_orig_col %in% names(mo2_data)) {
  mo2_data[[mo2_orig_col]] <- mo2_data[[mo2_col]]
  }

  ## --- Apply corrections ---
  cli::cli_h2("Applying MO2 corrections per chamber")

  for (ch_name in as.character(chambers_meta)) {
    ch_id <- as.numeric(ch_name)
    meta <- metadata[metadata[[chamber_col]] == ch_id, ]
    if (nrow(meta) != 1) {
      cli::cli_warn("No metadata found for chamber {ch_id}. Skipping.")
      next
    }

    mfish <- meta[[mass_col]]
    vresp <- meta[[vresp_col]]
    new_mass <- new_masses[[ch_name]]
    dens <- dens_map[[ch_name]]


    # --- Decide the effective new_mass and its status (robust NA-safe) ---
    if (!is.na(new_mass) && new_mass == mfish) {
      # Mass explicitly equal to the original
      mass_status <- "not changed"
      # new_mass stays as-is (equal to mfish)
    } else if (is.na(new_mass)) {
      # Mass not provided -> use original
      new_mass    <- mfish
      mass_status <- "not changed"
    } else {
      # Mass changed
      mass_status <- sprintf("%.4f kg", new_mass)
    }

    # --- Decide the effective density and its status (mirrors mass logic) ---
    # Normalize NA densities to 1.0 (unchanged)
    if (is.na(dens)) {
      dens        <- 1.0
      dens_status <- "not changed"
    } else if (dens == 1.0) {
      dens_status <- "not changed"
    } else {
      dens_status <- sprintf("%.3f kg/L", dens)
    }

    # Compute original and new rRespFish
    vfish_orig <- (mfish) / 1  #AquaResp assumes 1 kg/L (neutral fw buoyancy)
    vreal_orig <- vresp - vfish_orig
    rRespFish_orig <- vreal_orig / mfish

    vfish_new <- (new_mass) / dens
    vreal_new <- vresp - vfish_new
    rRespFish_new <- vreal_new / new_mass

    idx <- mo2_data[[chamber_col]] == ch_id
    if (!any(idx)) {
      cli::cli_warn("No rows found for chamber {ch_id} in mo2_data. Skipping.")
      next
    }

    old_mo2 <- mo2_data[[mo2_col]][idx]
    mo2_data[[mo2_col]][idx] <- old_mo2 / rRespFish_orig * rRespFish_new
    new_mo2 <- mo2_data[[mo2_col]][idx]

    valid <- is.finite(old_mo2) & old_mo2 != 0 & is.finite(new_mo2)
    pct_diff <-
      if (any(valid)) {
        mean((new_mo2[valid] - old_mo2[valid]) / old_mo2[valid] * 100, na.rm = TRUE)
        } else NA_real_

    # Determine colors for new values
    mass_color_new <-
      if (new_mass > mfish) "green" else
        if (new_mass < mfish) "red" else "cyan"
    rresp_color_new <-
      if (rRespFish_new > rRespFish_orig) "green" else
        if (rRespFish_new < rRespFish_orig) "red" else "black"
    dens_color_new <- if (!is.na(dens) && dens != 1.0) {
      if (dens > 1.0) "green" else "red"
    } else "cyan"

    diff_color <- if (!is.na(pct_diff)) {
      if (pct_diff > 0) "green" else if (pct_diff < 0) "red" else "blue"
    } else "blue"

    cli::cli_div(theme = list(
      span.mass_new = list(color = mass_color_new),
      span.rresp_new = list(color = rresp_color_new),
      span.dens_new = list(color = dens_color_new),
      span.diff = list(color = diff_color),
      h3 = list(color = "cyan", "font-weight" = "bold")
    ))

    # Header: original mass neutral, new mass colored
    cli::cli_h3("Chamber {ch_name}: {mfish} kg -> {.mass_new {mass_status}}")

    # rRespFish info: original neutral, new colored
    cli::cli_text("\tRespirometer ratio (rRespFish): {round(rRespFish_orig,2)} -> {.rresp_new {round(rRespFish_new,2)}} L/kg")

    # Density info if changed: original neutral, new colored
    if (!is.na(dens) && dens != 1.0) {
      cli::cli_text("\tDensity: original = 1.000 kg/L, new = {.dens_new {dens_status}}")
    }

    # Mean MO2 % difference colored
    cli::cli_text("\tMean MO2 % difference: {.diff {round(pct_diff,2)}}%")
    cat("\n")

    cli::cli_end()

  }

  mo2_data <- dplyr::relocate(.data = mo2_data,
                              .data[[mo2_orig_col]],
                              .after = dplyr::any_of(mo2_col))

  ## --- Append original/new mass and density columns if changes occurred ---
  mass_changed <- any(!is.na(new_masses) & new_masses != metadata[[mass_col]])
  dens_changed <- any(dens_map != 1.0)

  if (mass_changed) {
    mo2_data$Mass_orig <- metadata[[mass_col]][match(mo2_data[[chamber_col]], metadata[[chamber_col]])]
    mo2_data$Mass_new <- new_masses[as.character(mo2_data[[chamber_col]])]
  }

  if (dens_changed) {
    mo2_data$Density_orig <- 1.0
    mo2_data$Density_new <- dens_map[as.character(mo2_data[[chamber_col]])]
  }

  return(mo2_data)
}

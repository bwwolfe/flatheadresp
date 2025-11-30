#' Recalculate MO2 cycle outputs from AquaResp experiment data
#'
#' Recalculates MO2s and related metrics from raw cycle oxygen data with the
#' internal function `cycle_summary()`. It also returns the original AquaResp
#' values from `get_exp_mo2s()`
#'
#' @param path Character. Experiment directory.
#' @param chambers Numeric. The chamber numbers to include, negative subsets will exclude those chambers.
#' @param tolerance Numeric (non-negative). Tolerance for differences in
#'   `delta.po2` and `R.2` to flag corrections. Defaults to `0.001`, as there may
#'   be small rounding differences in the statistics calculated by this package
#'   vs the original AquaResp logic (coded in Python).
#' @param index_cols Character vector of index column names to join on.
#'   Defaults to `c("cycle", "chamber")`.
#' @param mass_col Character. Name of the chamber animal mass (kg) column in metadata. Defaults to
#' `"Mass.of.fish..kg"`. Note that the returned mass column is the cleaner `Mass_kg`
#' @param report Logical. If `TRUE` (default), a summary is printed to the console.
#' @details
#'
#' A logical `corrected` flag is inserted after the main columns, indicating
#'   if `delta.po2` or `R.2` differ by more than `tolerance`. The original values
#'   from AquaResp are appended at the end, suffixed with `"_uncorrected"`.
#'
#' \strong{Reporting:} When `calc_exp_mo2s()` runs, it reports several quality-control metrics to the console:
#'
#' \itemize{
#'   \item \strong{Number of cycles corrected:} Cycles adjusted due to bugs in AquaResp code.
#'
#'   \item \strong{Cycles with minimum \eqn{PO_{2} \le 0\%} saturation:}
#'   These values are treated as missing and removed before calculating \eqn{MO_{2}} or \eqn{R^{2}}.
#'
#'   \item \strong{Cycles with minimum \eqn{PO_{2} < 0\%}:}
#'   Often corresponds to erroneous values (e.g., \eqn{-300\%}) when an oxygen probe was not functioning.
#'   If present, these cycles may indicate disturbances worth investigating.
#'
#'   \item \strong{Cycles with corrected \eqn{R^{2}} values:}
#'   Corrections occur when differences exceed a tolerance (default = \eqn{0.001}).
#'   Even without data removal, minor discrepancies can arise between AquaResp's Python calculations
#'   and this package's implementation.
#'
#'   \item \strong{Additional metrics:}
#'     \itemize{
#'       \item Greatest rate of missing \eqn{PO_{2}} measurements in a cycle
#'             (typically \eqn{\sim 0.1\%}; higher values may indicate the computer is missing
#'             a lot of oxygen sensor data and may need review).
#'       \item Number of cycles where \emph{all} chambers have \eqn{R^{2} < 0.95},
#'             which may indicate disturbances or insufficient measurement window durations,
#'             if they are common. This metric is also provided on a per-chamber basis.
#'     }
#' }
#' @return A data.frame with the column order described above.
#' @export
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment", package = "flatheadresp")
#' calc_exp_mo2s(path = exp_dir_path)

calc_exp_mo2s <- function(path,
                          chambers = NULL,
                          tolerance = 0.001,
                          report = TRUE,
                          index_cols = c(cycle = "cycle", chamber = "chamber"),
                          mass_col = "Mass.of.fish..kg") {
  # --- basic checks ---
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is.numeric(tolerance), length(tolerance) == 1, tolerance >= 0)
  stopifnot(is.character(index_cols), length(index_cols) >= 1)

  # --- read inputs ---
  orig <- get_exp_mo2s(path)
  corr <- cycle_summary(path, long = TRUE)

  # --- ensure required index columns exist ---
  if (!all(index_cols %in% names(orig))) {
    stop("Missing index column(s) in get_exp_mo2s(): ",
         paste(setdiff(index_cols, names(orig)), collapse = ", "))
  }
  if (!all(index_cols %in% names(corr))) {
    stop("Missing index column(s) in cycle_summary(): ",
         paste(setdiff(index_cols, names(corr)), collapse = ", "))
  }

  chamber_col <- index_cols["chamber"]
  cycle_col   <- index_cols["cycle"]

  # --- resolve chambers (allow negatives for exclusion) ---
  all_chambers <- sort(unique(orig[[chamber_col]]))
  if (is.null(chambers)) {
    sel_chambers <- all_chambers
  } else {
    chambers <- as.integer(chambers)
    if (any(chambers < 0)) {
      exclude <- abs(chambers[chambers < 0])
      sel_chambers <- setdiff(all_chambers, exclude)
    } else {
      sel_chambers <- chambers
    }
    missing <- setdiff(abs(chambers), all_chambers)
    if (length(missing) > 0) {
      cli::cli_abort("Invalid chamber(s): {missing}. Available: {all_chambers}")
    }
  }

  # --- filter data to selected chambers ---
  orig <- orig[orig[[chamber_col]] %in% sel_chambers, ]
  corr <- corr[corr[[chamber_col]] %in% sel_chambers, ]

  # --- determine name sets ---
  orig_non_index <- setdiff(names(orig), index_cols)
  first_three_uncorrected <- utils::head(orig_non_index, 3)
  overlap <- intersect(setdiff(names(orig), index_cols),
                       setdiff(names(corr), index_cols))

  # --- align corrected rows to original via index_cols ---
  corr_aligned <- merge(orig[, index_cols, drop = FALSE],
                        corr,
                        by = index_cols,
                        all.x = TRUE,
                        sort = FALSE)

  # --- build main output ---
  out <- orig[, index_cols, drop = FALSE]
  if (length(first_three_uncorrected) > 0) {
    out[first_three_uncorrected] <- orig[first_three_uncorrected]
  }



  remaining_orig_cols <- setdiff(orig_non_index, first_three_uncorrected)
  for (col in remaining_orig_cols) {
    out[[col]] <- if (col %in% names(corr_aligned)) corr_aligned[[col]] else orig[[col]]
  }
  corr_only_cols <- setdiff(setdiff(names(corr_aligned), index_cols), names(orig))
  if (length(corr_only_cols) > 0) {
    out[corr_only_cols] <- corr_aligned[corr_only_cols]
  }

  # --- corrected flag ---
  dpo2_corr <- corr_aligned[["delta.po2"]]
  dpo2_orig <- orig[["delta.po2"]]
  r2_corr   <- corr_aligned[["R.2"]]
  r2_orig   <- orig[["R.2"]]

  diff_dpo2 <- if (!is.null(dpo2_corr) && !is.null(dpo2_orig)) {
    na2false(abs(dpo2_corr - dpo2_orig) > tolerance)
  } else rep(FALSE, nrow(out))

  diff_r2 <- if (!is.null(r2_corr) && !is.null(r2_orig)) {
    na2false(abs(r2_corr - r2_orig) > tolerance)
  } else rep(FALSE, nrow(out))

  corrected_flag <- diff_dpo2 | diff_r2

  # --- original overlap columns appended ---
  orig_overlap_df <- if (length(overlap) > 0) {
    ss <- orig[overlap]
    names(ss) <- paste0(overlap, "_uncorrected")
    ss
  } else NULL

  final <- cbind(out, corrected = corrected_flag, orig_overlap_df)

  if(isTRUE(report)) {
    # --- summary metrics ---
    min_col <- if ("minimum.po2" %in% names(orig)) "minimum.po2" else
      if ("min.po2" %in% names(orig)) "min.po2" else NA_character_
    min_unc <- if (!is.na(min_col)) orig[[min_col]] else rep(NA_real_, nrow(orig))

    chambers <- sel_chambers
    n_cycles_total <- length(unique(orig[[cycle_col]]))
    n_corrected <- length(unique(orig[[cycle_col]][corrected_flag]))
    n_min_le0 <- length(unique(orig[[cycle_col]][!is.na(min_unc) & min_unc <= 0]))
    n_min_lt0 <- length(unique(orig[[cycle_col]][!is.na(min_unc) & min_unc < 0]))
    r2_changed_rows <- corrected_flag & diff_r2
    n_r2_diff <- length(unique(orig[[cycle_col]][r2_changed_rows]))

    r2_low_rows <- if (!is.null(r2_corr)) na2false(r2_corr < 0.95) else rep(FALSE, nrow(out))

    pct0_vec <- if ("pct0" %in% names(corr_aligned)) corr_aligned[["pct0"]] else NA_real_
    max_pct0 <- suppressWarnings(max(pct0_vec, na.rm = TRUE))
    if (!is.finite(max_pct0)) max_pct0 <- NA_real_
    pct0_str <- if (is.na(max_pct0)) "NA" else sprintf("%.2f%%", max_pct0 * 100)

    # --- CLI reporting ---
    cli::cli_div(theme = list(
      "span.normie"  = list(color = "black"),
      "span.total"   = list(color = "black", "font-weight" = "bold"),
      "span.ok"      = list(color = "green", "font-weight" = "bold"),
      "span.bad"     = list(color = "red",   "font-weight" = "bold"),
      "span.neutral" = list(color = "blue",  "font-weight" = "bold"),
      h3             = list(color = "cyan",  "font-weight" = "bold")
    ))

    header_text <- ifelse(is.null(chambers),
                          "MO2 Summary for all chambers {.strong {paste(chambers, collapse = ', ')}}",
                          "MO2 Summary for selected chambers: {.strong {paste(chambers, collapse = ', ')}}")

    cli::cli_h2(header_text)

    cli::cli_h3("With a po2 and r2 difference tolerance of {tolerance}:")

    pad_bracket <- function(x, width, cls = NULL) {
      markup <- if (is.null(cls)) sprintf("{%s}", x) else sprintf("{.%s %s}", cls, x)
      pad <- max(1, width - nchar(as.character(x)))
      paste0(markup, strrep("\u00A0", pad))
    }

    counts <- c(n_cycles_total, n_corrected, n_min_le0, n_min_lt0, n_r2_diff)
    count_width <- 2 + max(2, nchar(as.character(counts)))

    cli::cli_bullets(c(
      "*" = sprintf("%scycle(s) total.", pad_bracket(n_cycles_total, count_width, "neutral")),
      "*" = sprintf("%scycle(s) had corrected values.", pad_bracket(n_corrected, count_width, if (n_corrected > 0) "ok" else "neutral")),
      "*" = sprintf("%scycle(s) had an uncorrected minimum.po2 <= 0.", pad_bracket(n_min_le0, count_width, if (n_min_le0 > 0) "bad" else "ok")),
      "*" = sprintf("%scycle(s) had an uncorrected minimum.po2 < 0.", pad_bracket(n_min_lt0, count_width, if (n_min_lt0 > 0) "bad" else "ok")),
      "*" = sprintf("%scycle(s) R^2 were corrected (changed > %s).", pad_bracket(n_r2_diff, count_width, if (n_r2_diff > 0) "bad" else "neutral"), tolerance)
    ))

    if (pct0_str == "NA") {
      cli::cli_alert_info("Max percent of pO2 measurements missing (<= 0) in a cycle: NA")
    } else {
      cli::cli_h3(sprintf("Max percent of pO2 measurements missing (<= 0) in a cycle: {.%s %s}", if (max_pct0 > 0) "bad" else "ok", pct0_str))
    }


    # Build header text
    cycles_header_text <- ifelse(
      is.null(chambers),
      "in which all chambers' corrected R.2 < 0.95.",
      "in which all selected chambers' corrected R.2 < 0.95."
    )

    crit_idx <- which(r2_low_rows)

    if (length(crit_idx) == 0) {
      # Combine header and bullet into one string
      cli::cli_h3(sprintf("%s {.neutral {0}} cycles", cycles_header_text))
    } else {
      cycles_meeting <- orig[[cycle_col]][crit_idx]
      chambers_meeting <- orig[[chamber_col]][crit_idx]
      ch_by_cycle <- split(chambers_meeting, cycles_meeting)
      cycles_all <- names(Filter(function(x) {
        ux <- sort(unique(x))
        identical(ux, sort(unique(chambers)))
      }, ch_by_cycle))

      n_all <- length(cycles_all)
      cls_all <- if (n_all > 0) "bad" else "ok"

      # Build the inline text
      if (n_all <= 10) {
        inline_text <- sprintf("{.%s {%s} cycle(s)}", cls_all, n_all)
      } else {
        inline_text <- sprintf("{.%s {%s} cycle(s)}", cls_all, n_all)
      }

      # Combine header and inline text
      cli::cli_h3(sprintf("%s %s", inline_text, cycles_header_text))
    }
    # detail the cycles that had low r2, if this any
    if (n_all > 0)  {
      if (n_all <= 10) {
        cli::cli_bullets(sprintf("Cycle(s) {.%s {%s}}", cls_all, n_all, paste(cycles_all, collapse = ", ")))
      } else {
        cli::cli_bullets(sprintf("Cycle(s) {.%s {%s}}", cls_all, n_all, paste0(paste(cycles_all[1:10], collapse = ", "), " ...")))
      }
    }
    #  CHAMBER-LEVEL: cycles with corrected R.2 < 0.95 =====
    cli::cli_h3("Cycles with corrected R.2 < 0.95 by chamber:")

    for (ch in chambers) {
      idx <-
        which(r2_low_rows & orig[[chamber_col]] == ch)  # no min_unc gating
      cycles <- orig[[cycle_col]][idx]
      n_ch <- length(idx)
      cls_ch <- if (n_ch > 0)
        "bad"
      else
        "neutral"

      if (n_ch == 0) {
        cli::cli_bullets(sprintf("Chamber %s: {.%s {%s}} cycles", ch, cls_ch, n_ch))
      } else if (n_ch <= 10) {
        cli::cli_bullets(sprintf(
          "Chamber %s: {.%s {%s}} cycle(s) (%s)",
          ch,
          cls_ch,
          n_ch,
          paste(cycles, collapse = ", ")
        ))
      } else {
        cli::cli_bullets(sprintf(
          "Chamber %s: {.%s {%s}} cycle(s) (%s)",
          ch,
          cls_ch,
          n_ch,
          paste0(paste(cycles[1:10], collapse = ", "), " ...")
        ))
      }
    }

    cat("\n")
  }

  metadata <- get_exp_metadata(path)

  final$Mass_kg <-
    metadata[[mass_col]][match(final[[chamber_col]], metadata[[chamber_col]])]

  as.data.frame(final)
}

# helper: convert NA logicals to FALSE
na2false <- function(x) {
  x[is.na(x)] <- FALSE
  x
}

#' Q10 temperature correction of metabolic rate (MO2)
#'
#' @description
#' Applies a Q10 correction to a column of MO2 values by estimating the original
#' temperature for each row from a time–temperature reference dataset. By default,
#' it uses the package dataset `temp_correct_data`, but you can supply a custom
#' dataset via `temp_data`. The estimated temperature is appended as a new
#' column (default: `est_temp`), and the corrected MO2 is inserted immediately
#' after the original MO2 column. Note: if an experiment longer than 25 hours is
#' converted, beyond this the temperature estimated will be flat (`temp_correct_data`
#' is 25 hrs long).
#'
#' @details
#' The correction is:
#' \deqn{\mathrm{MO2}_{\mathrm{new}} = \mathrm{MO2}\times Q_{10}^{(T_{\mathrm{new}} - T_{\mathrm{old}})/10}}
#'
#' where \code{T_old} is the estimated temperature taken from `temp_data`
#' (or the package dataset `temp_correct_data` if `temp_data` is `NULL`)
#' at the **nearest** \code{TIME.HOURS} to each row's time. No interpolation is used.
#' If times fall outside the range of the reference, the nearest end-point is used.
#'
#' @param mo2_data A \code{data.frame} (or tibble) containing at least a time column
#'   and an MO2 column.
#' @param new_temp Numeric. The target temperature (°C) to correct to. Can be a single
#'   value or a numeric vector the same length as \code{nrow(mo2_data)}.
#' @param q10 Numeric scalar giving the Q10 value to use (default: \code{2}).
#' @param time_col Character name of the time column in \code{mo2_data} that is in
#'   decimal hours since experiment start (default: \code{"TIME.HOURS"}).
#' @param mo2_col Character name of the MO2 column in \code{mo2_data} to correct
#'   (default: \code{"MO2"}).
#' @param est_temp_col Name for the appended estimated temperature column in the output
#'   (default: \code{"est_temp"}).
#' @param out_col Optional name for the corrected MO2 column. If \code{NULL} (default),
#'   it becomes \code{"<mo2_col>_q10_<new_temp>"} for constant \code{new_temp}, or
#'   \code{"<mo2_col>_q10"} if \code{new_temp} varies by row.
#' @param temp_data Optional. A data.frame/tibble with columns \code{TIME.HOURS} and
#'   \code{est_temp} to use as the reference time–temperature series. If \code{NULL},
#'   the function loads the package dataset \code{temp_correct_data}.
#'
#' @return
#' The input \code{mo2_data} with:
#' \itemize{
#'   \item a new column \code{est_temp_col} giving the estimated original temperature (°C),
#'   \item a new numeric column with the Q10-corrected MO2, inserted immediately after \code{mo2_col}.
#' }
#' @importFrom stats approx
#' @examples
#' \dontrun{
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment", package = "flatheadresp")
#' mo2_data <- calc_exp_mo2s(path = exp_dir_path)
#'
#' # Default: uses the package dataset temp_correct_data
#' res1 <- correct_for_temperature(mo2_data, new_temp = 16, q10 = 2)
#'
#' library(ggplot2)
#' ggplot(res1|>subset(chamber != 4)) +
#'   geom_path(aes(TIME.HOURS, MO2, colour = factor(chamber))) +
#'   geom_path(aes(TIME.HOURS, MO2_q10_16, colour = factor(chamber), group = chamber),
#'             linetype = "dotted") +
#'   theme_classic()

#' # Custom reference (must have TIME.HOURS and est_temp)
#' custom_ref <- data.frame(TIME.HOURS = seq(0, 25, by = 0.1),
#'                          est_temp   = 16 + 0.02 * seq(0, 25, by = 0.1))
#' res2 <- correct_for_temperature(mo2_data, new_temp = 16, q10 = 2,
#'                                 temp_data = custom_ref)
#'}
#' @seealso \code{\link{allometric_correct}}, \code{\link{temp_correct_data}}
#' @export
correct_for_temperature <- function(mo2_data,
                                    new_temp,
                                    q10,
                                    time_col = "TIME.HOURS",
                                    mo2_col = "MO2",
                                    est_temp_col = "est_temp",
                                    out_col = NULL,
                                    temp_data = NULL) {
  # --- Basic checks on mo2_data ---
  if (!is.data.frame(mo2_data)) stop("`mo2_data` must be a data.frame or tibble.")
  if (!time_col %in% names(mo2_data)) stop(sprintf("Time column `%s` not found.", time_col))
  if (!mo2_col %in% names(mo2_data))  stop(sprintf("MO2 column `%s` not found.", mo2_col))
  if (!is.numeric(mo2_data[[time_col]])) stop(sprintf("`%s` must be numeric (decimal hours).", time_col))
  if (!is.numeric(mo2_data[[mo2_col]]))  stop(sprintf("`%s` must be numeric.", mo2_col))

  if (!is.numeric(q10) || length(q10) != 1L || is.na(q10) || q10 <= 0) {
    stop("`q10` must be a single positive numeric value.")
  }
  if (!(length(new_temp) == 1L || length(new_temp) == nrow(mo2_data))) {
    stop("`new_temp` must be a single number or a vector the same length as `mo2_data`.")
  }

  # --- Obtain reference temperature series ---
  ref <- NULL
  if (is.null(temp_data)) {
    # Try to load package dataset temp_correct_data into a private env
    # Works regardless of LazyData setting
    tmp_env <- new.env(parent = emptyenv())
    pkgname <- utils::packageName()
    # If called outside a package namespace, fall back to current search path
    if (!is.null(pkgname)) {
      utils::data("temp_correct_data", package = pkgname, envir = tmp_env)
    } else {
      utils::data("temp_correct_data", envir = tmp_env)
    }
    if (!exists("temp_correct_data", envir = tmp_env, inherits = FALSE)) {
      stop("`temp_correct_data` dataset not found. ",
           "Either provide `temp_data=` or ensure the package dataset is available.")
    }
    ref <- get("temp_correct_data", envir = tmp_env)
  } else {
    ref <- temp_data
  }

  if (!is.data.frame(ref)) stop("`temp_data` must be a data.frame/tibble when provided.")
  required_cols <- c("TIME.HOURS", "est_temp")
  if (!all(required_cols %in% names(ref))) {
    stop("Reference data must contain columns: ", paste(required_cols, collapse = ", "))
  }
  if (!is.numeric(ref$TIME.HOURS) || !is.numeric(ref$est_temp)) {
    stop("`TIME.HOURS` and `est_temp` in reference data must be numeric.")
  }

  # --- Nearest-neighbour temperature assignment ---
  ord   <- order(ref$TIME.HOURS)
  t_x   <- ref$TIME.HOURS[ord]
  t_y   <- ref$est_temp[ord]
  times <- mo2_data[[time_col]]

  # find estimated temperature by linear interpolation (useful if the reference
  # temp data is low temporal resolution)

  est_temp_vals <-
    approx(x = t_x, y = t_y, xout = times, ties = "ordered", rule = 2)$y

  # Append estimated temperature (overwrite if name collides)
  mo2_data[[est_temp_col]] <- est_temp_vals

  # --- Q10 correction ---
  if (length(new_temp) == 1L) new_temp <- rep(new_temp, length(est_temp_vals))
  factor    <- q10 ^ ((new_temp - est_temp_vals) / 10)
  corrected <- mo2_data[[mo2_col]] * factor

  # --- Name & insert corrected column after mo2_col ---
  if (is.null(out_col)) {
    if (length(unique(new_temp)) == 1L) {
      nt_label <- sub("\\.?0+$", "", format(unique(new_temp), trim = TRUE))
      out_col  <- paste0(mo2_col, "_q10_", nt_label)
    } else {
      out_col  <- paste0(mo2_col, "_q10")
    }
  }

  mo2_idx <- match(mo2_col, names(mo2_data))
  left    <- mo2_data[, 1:mo2_idx, drop = FALSE]
  right   <- if (mo2_idx < ncol(mo2_data)) mo2_data[, (mo2_idx + 1):ncol(mo2_data), drop = FALSE] else NULL
  new_df  <- stats::setNames(data.frame(corrected, check.names = FALSE), out_col)
  out     <- if (is.null(right)) cbind(left, new_df) else cbind(left, new_df, right)

  # Preserve tibble class if present
  if (inherits(mo2_data, "tbl_df")) {
    out <- tibble::as_tibble(out)
  }
  out
}

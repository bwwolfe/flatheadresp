#' Allometric correction of mass-specific metabolic rate (MO2)
#'
#' @description
#' Computes an allometrically corrected metabolic rate from a data frame that
#' contains body mass (default column: \code{Mass_kg}) and mass-specific metabolic
#' rate (default column: \code{MO2}, in \eqn{\mathrm{mg}\ \mathrm{O_2}\ \mathrm{h}^{-1}\ \mathrm{kg}^{-1}}).
#' The correction uses an exponent \code{b} (the allometric scaling exponent) to yield
#' MO2 values in units of \eqn{\mathrm{mg}\ \mathrm{O_2}\ \mathrm{h}^{-1}\ \mathrm{kg}^{-b}}.
#'
#' @details
#' The corrected metabolic rate is computed as:
#'
#' \deqn{\mathrm{MO2}_b = \mathrm{MO2} \times \mathrm{Mass}^{(1 - b)}}
#'
#' where \code{MO2} is mass-specific metabolic rate
#' (\eqn{\mathrm{mg}\ \mathrm{O_2}\ \mathrm{h}^{-1}\ \mathrm{kg}^{-1}}),
#' \code{Mass} is body mass in kilograms, and \code{b} is the allometric scaling exponent.
#' This produces output in \eqn{\mathrm{mg}\ \mathrm{O_2}\ \mathrm{h}^{-1}\ \mathrm{kg}^{-b}},
#' aligning the mass exponent with \code{b}.
#'
#' A new column is created and inserted \emph{immediately after} the original \code{MO2}
#' column. The new column name appends the exponent to the MO2 base name, e.g.
#' \code{MO2.79} when \code{b = 0.79}. Trailing zeros are trimmed (e.g., \code{1.00} becomes \code{1}).
#'
#' @param mo2_data A \code{data.frame} containing at least the mass and MO2 columns.
#' (Output from `calc_exp_mo2s()` or `fix_exp_mo2s()` would work.)
#' @param b A single numeric value giving the allometric scaling exponent.
#'   Defaults to \code{1} i.e., isometric scaling, equivalent of normal mass-specific units.
#' @param mass_col A character string giving the column name for body mass
#'   (default: \code{"Mass_kg"}). Must be numeric and in kilograms.
#' @param mo2_col A character string giving the column name for mass-specific
#'   metabolic rate (default: \code{"MO2"}). Must be numeric and in units
#'   \eqn{\mathrm{mg}\ \mathrm{O_2}\ \mathrm{h}^{-1}\ \mathrm{kg}^{-1}}.
#' @param beta Optional alias for the exponent \code{b}. If provided (non-\code{NULL}),
#'   it overrides \code{b}. Use this if you prefer the more explicit name in scripts.
#' @param mass_specific Logical. If `TRUE` (default), `MO2` is expected to be in units of
#'   \eqn{\mathrm{mg}\ \mathrm{O_2}\ \mathrm{h}^{-1}\ \mathrm{kg}} i.e. the typical format used
#'   by AquaResp and the rest of this package. If `FALSE`, units are \eqn{\mathrm{mg}\ \mathrm{O_2}\ \mathrm{h}^{-1}}
#' @return
#' Returns the input \code{mo2_data} with one additional numeric column named
#' \code{paste0(mo2_col, ".", format(b))} containing the allometrically corrected MO2
#' (\eqn{\mathrm{mg}\ \mathrm{O_2}\ \mathrm{h}^{-1}\ \mathrm{kg}^{-b}}), inserted directly after \code{mo2_col}.
#'
#' @section Units:
#' \itemize{
#'   \item Input \code{MO2}: \eqn{\mathrm{mg}\ \mathrm{O_2}\ \mathrm{h}^{-1}\ \mathrm{kg}^{-1}}
#'   \item Input \code{Mass_kg}: \eqn{\mathrm{kg}}
#'   \item Output \code{MO2.b}: \eqn{\mathrm{mg}\ \mathrm{O_2}\ \mathrm{h}^{-1}\ \mathrm{kg}^{-b}}
#' }
#'
#' @note
#' This function does \strong{not} round values; apply rounding at presentation/export
#' time if desired (e.g., \code{round()} or \code{signif()}).
#'
#' @examples
#' \dontrun{
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment", package = "flatheadresp")
#' mo2_data <- calc_exp_mo2s(path = exp_dir_path)
#'
#' ## Default (b = 1): returns MO2_1 which equals MO2 in units mg O2 h^-1 kg^-1
#' head(allometric_correct(mo2_data))[,c(1,2,)]
#'
#' ## Using a typical fish scaling exponent, see Clarke and Johnston (1999)
#' allometric_correct(mo2_data, b = 0.79)
#'}
#' @references
#' Clarke, A. & Johnston, N. M. (1999). Scaling of metabolic rate with body mass
#' and temperature in teleost fish. \emph{Journal of Animal Ecology}, \strong{68}(5), 893–905.
#' doi:10.1046/j.1365-2656.1999.00337.x.
#' @seealso
#' \code{\link[base]{round}}, \code{\link[base]{signif}}
#'
#' @aliases allometric_correction
#' @export

allometric_correct <- function(mo2_data,
                               b = 1,
                               mass_col = "Mass_kg",
                               mo2_col = "MO2",
                               beta = NULL,
                               mass_specific = TRUE) {
  # Alias: if beta specified, it overrides b
  if (!is.null(beta)) b <- beta

  # --- Basic checks ---
  if (!is.data.frame(mo2_data)) stop("`mo2_data` must be a data.frame or tibble.")
  if (!mass_col %in% names(mo2_data)) stop(sprintf("Mass column `%s` not found.", mass_col))
  if (!mo2_col %in% names(mo2_data)) stop(sprintf("MO2 column `%s` not found.", mo2_col))

  if (!is.numeric(mo2_data[[mass_col]])) stop(sprintf("`%s` must be numeric.", mass_col))
  if (!is.numeric(mo2_data[[mo2_col]]))  stop(sprintf("`%s` must be numeric.", mo2_col))
  if (!is.numeric(b) || length(b) != 1L || is.na(b)) stop("`b` must be a single numeric value.")

  # --- Compute corrected MO2: mg O2 per hour per kg^b ---
  if(mass_specific){
  # MO2_b = MO2 * Mass^(1 - b)
  corrected <- mo2_data[[mo2_col]] * (mo2_data[[mass_col]]^(1 - b))
  } else {
    corrected <- mo2_data[[mo2_col]] / (mo2_data[[mass_col]]^b)
  }
  # --- Build column name like 'MO2_0.79' ---
  # Keep user's b as string, trimming trailing zeros
  b_label <- sub("\\.?0+$", "", format(b, trim = TRUE))  # e.g., 0.790 -> "0.79"; 1.0 -> "1"
  new_name <- paste0(mo2_col, "_", b_label)

  # --- Insert right after MO2 ---
  mo2_idx <- match(mo2_col, names(mo2_data))

  # Handle case where MO2 is the last column
  if (mo2_idx == ncol(mo2_data)) {
    mo2_data[[new_name]] <- corrected
    return(mo2_data)
  }

  # General case: splice in after MO2
  left  <- mo2_data[ , 1:mo2_idx, drop = FALSE]
  right <- mo2_data[ , (mo2_idx + 1):ncol(mo2_data), drop = FALSE]
  new_col_df <- stats::setNames(data.frame(corrected, check.names = FALSE), new_name)

  out <- cbind(left, new_col_df, right)
  # Preserve tibble class if present
  if (inherits(mo2_data, "tbl_df")) {
    out <- tibble::as_tibble(out)
  }
  out
}

#' Pretty-print a matrix with colour-mapped values, flags, and a legend
#'
#' Prints a numeric matrix to the console with row/column labels, padded spacing,
#' and per-cell colour based on value magnitude (scaled min \eqn{\rightarrow} max).
#' Two optional flag matrices can add CLI styles to each value:
#'   - \code{r2_mat}: numeric \eqn{R^2} values converted internally to a binary mask using \code{r2_threshold}.
#'   - \code{missing_mat}: logical/numeric mask for missing \eqn{O_2} values.
#'
#' This version:
#'  - Ensures all values are mapped to a visible colour even on low-colour consoles.
#'  - Handles \eqn{R^2} thresholding internally and auto-generates the legend label.
#'  - In high-colour mode (\eqn{\ge 256} colours), applies a dark base background with value-dependent colour.
#'  - Supports an optional extra column printed at the end of each row.
#'
#' @param x Numeric matrix. Values to print (NAs allowed).
#' @param r2_mat Optional numeric matrix (or vector) of \eqn{R^2} values with the same dimensions as \code{x}.
#'   Converted to a logical mask where \eqn{R^2 <} \code{r2_threshold}.
#' @param missing_mat Optional logical/numeric matrix (or vector) with the same dimensions as \code{x}.
#'   Non-zero (or \code{TRUE}) entries mark cells for styling.
#' @param r2_threshold Numeric. Threshold for \eqn{R^2} flagging (default \code{0.95}).
#' @param digits Integer. Decimal places to print (default \code{2}).
#' @param pad Integer. Left padding spaces per cell (default \code{1}).
#' @param palette Character vector of hex colours for the gradient (low \eqn{\rightarrow} high),
#'   used when 256/truecolor is available. Default:
#'   \code{c("#F49898", "#EF6F6C", "#D92523", "#C90202")}.
#' @param na_string Character. Representation for \code{NA} cells (default \code{"NA "}).
#' @param show_border Logical. Draw a vertical separator after the row label (default \code{TRUE}).
#' @param legend Logical. Print a legend below the matrix (default \code{TRUE}).
#' @param legend_cols Integer. Number of blocks in the legend colour bar (default \code{48}).
#' @param missing_label Character. Label for the missing flag in the legend (default \code{"missing pO2s"}).
#' @param show_flag_counts Logical. Include counts of flagged cells in the legend (default \code{TRUE}).
#' @param extra_col Vector with length equal to \code{nrow(x)}. When supplied, an extra column is printed
#'   to the right of the matrix, one entry per row.
#' @param extra_col_label Character. Optional header label for \code{extra_col}. If provided, it appears
#'   in the header row aligned with the extra column.
#'
#' @return Invisibly returns \code{NULL}.
#' @importFrom cli make_ansi_style combine_ansi_styles col_white style_underline style_bold
#' @export
#'
#' @details
#' **Colour and styling:**
#'
#' The function uses the \code{cli} package for ANSI styling. Colour strategy adapts to the terminal:
#'
#' - **High-colour mode** (\eqn{\ge 256} colours):
#'   - Values are mapped via a Lab-space gradient built from \code{palette}.
#'   - A dark base background (\code{"#121212"}) is applied behind value cells for consistent contrast.
#'   - Missing-flag style: white background with black text.
#'   - \eqn{R^2}-flag style: bold and underline and italic with green text.
#'   - Both flags: combined styles including a dark green accent.
#'   - \code{NA} cells: bright yellow background with blurred blue text.
#'
#' - **Low-colour mode** (< 256 colours):
#'   - Value scaling is disabled; all numeric cells use black background with white text.
#'   - \eqn{R^2}-flagged cells: blue background with white text.
#'   - Missing-flagged cells: yellow background with black text.
#'   - Both flags: magenta background with white text.
#'   - \code{NA} cells: red background with blue text.
#'
#' To ensure full colour support when auto-detection is uncertain, you can set:
#' \preformatted{
#' options(cli.default_num_colors = 256L)
#' }
#'
#' **Flagging:**
#'
#' - \code{r2_mat} (matrix or vector) is validated to match \code{x}. It is converted into a binary mask
#'   where \eqn{R^2 <} \code{r2_threshold} and \code{!is.na(r2_mat)}.
#' - \code{missing_mat} is coerced into a logical mask of the same dimensions as \code{x}; non-zero/\code{TRUE}
#'   entries mark cells for the missing style.
#'
#' **Extra column:**
#'
#' - If \code{extra_col} is supplied, a right-hand column is printed per row.
#' - If \code{extra_col_label} is supplied, a header label is printed above that column.
#'
#' **Input validation:**
#'
#' - \code{x} must be a numeric matrix; otherwise, the function stops with an error.
#' - \code{r2_mat} and \code{missing_mat} (if provided as vectors) are reshaped to match \code{nrow(x)} \eqn{\times} \code{ncol(x)}.
#' - Dimension names of \code{x} are propagated to internal masks when present.
#'
#' @examples
#' set.seed(12)
#' x <- matrix(runif(20, 0, 1), nrow = 5,
#'             dimnames = list(paste0("row", 1:5), paste0("col", 1:4)))
#' r2 <- matrix(runif(20, 0.90, 1.0), nrow = 5)   # mostly high R^2
#' miss <- matrix(sample(c(FALSE, TRUE), 20, replace = TRUE, prob = c(0.85, 0.15)), nrow = 5)
#'
#' print_matrix_exp(
#'   x,
#'   r2_mat = r2,
#'   missing_mat = miss,
#'   digits = 2,
#'   pad = 1
#' )
print_matrix_exp <- function(
    x,
    r2_mat = NULL,
    missing_mat = NULL,
    r2_threshold = 0.95,
    digits = 2,
    pad = 1,
    palette = c("#F49898", "#EF6F6C", "#D92523", "#C90202"),
    na_string = "NA ",
    show_border = TRUE,
    legend = TRUE,
    legend_cols = 48,
    missing_label = "missing pO2s",
    show_flag_counts = TRUE,
    extra_col = NULL,
    extra_col_label = NULL
) {

  # ---- Validate input matrix ----
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("x must be a numeric matrix.")
  }
  nr <- nrow(x); nc <- ncol(x)

  # ---- Validate extra_col ----
  if (!is.null(extra_col)) {
    if (length(extra_col) != nr) {
      stop("extra_col must have the same length as the number of rows in x.")
    }
    extra_col <- as.character(extra_col)
    extra_width <- max(nchar(extra_col), nchar(extra_col_label %||% ""), 1) + 6
    extra_col_label <- style_bold(extra_col_label)
  } else {
    extra_width <- 0
  }

  # ---- Handle missing flag matrix ----
  coerce_flag <- function(f) {
    if (is.null(f)) {
      return(matrix(FALSE, nrow = nr, ncol = nc, dimnames = dimnames(x)))
    }
    if (!is.matrix(f)) {
      if (length(f) != nr * nc) stop("flag must match dimensions of x.")
      f <- matrix(f, nrow = nr, ncol = nc, byrow = FALSE)
    }
    if (!all(dim(f) == dim(x))) stop("flag must have same dimensions as x.")
    f <- f != 0
    if (!is.null(dimnames(x))) dimnames(f) <- dimnames(x)
    f
  }

  # ---- Convert r2 matrix to binary mask using threshold ----
  if (is.null(r2_mat)) {
    f1 <- matrix(FALSE, nrow = nr, ncol = nc, dimnames = dimnames(x))
  } else {
    if (!is.matrix(r2_mat)) {
      if (length(r2_mat) != nr * nc) stop("r2_mat must match dimensions of x.")
      r2_mat <- matrix(r2_mat, nrow = nr, ncol = nc, byrow = FALSE)
    }
    if (!all(dim(r2_mat) == dim(x))) stop("r2_mat must have same dimensions as x.")
    if (!is.numeric(r2_mat)) stop("r2_mat must be numeric.")
    f1 <- (r2_mat < r2_threshold) & !is.na(r2_mat)
    if (!is.null(dimnames(x))) dimnames(f1) <- dimnames(x)
  }

  f2 <- coerce_flag(missing_mat)

  r2_label <- sprintf("R2 < %s", formatC(r2_threshold, format = "f", digits = 2))

  # ---- Prepare labels and formatting ----
  rn <- rownames(x); if (is.null(rn)) rn <- as.character(seq_len(nr))
  cn <- colnames(x); if (is.null(cn)) cn <- as.character(seq_len(nc))
  fmt_num <- function(v) ifelse(is.na(v), na_string, formatC(v, format = "f", digits = digits))
  formatted <- matrix(fmt_num(x), nrow = nr, ncol = nc)
  cell_width <- max(nchar(formatted), nchar(na_string), 1) + 2 * pad
  row_lab_width <- max(nchar(rn)) + 2L

  # ---- Compute value range for scaling ----
  finite_vals <- as.numeric(x[is.finite(x)])
  vmin <- if (length(finite_vals)) min(finite_vals) else 0
  vmax <- if (length(finite_vals)) max(finite_vals) else 1
  span <- if (vmax > vmin) (vmax - vmin) else 1
  norm01 <- function(v) pmax(0, pmin(1, (v - vmin) / span))

  # ---- Determine colour strategy ----
  ncols_supported <- cli::num_ansi_colors()
  has_256_plus <- ncols_supported >= 256

  if (has_256_plus) {
    ramp <- grDevices::colorRamp(palette, space = "Lab")
    to_hex <- function(t) {
      rgb <- ramp(t) / 255
      grDevices::rgb(rgb[,1], rgb[,2], rgb[,3])
    }
    style_cache <- new.env(parent = emptyenv())
    style_for_value <- function(t) {
      hex <- to_hex(t)
      s <- style_cache[[hex]]
      if (is.null(s)) {
        s <- make_ansi_style(hex, colors = ncols_supported)
        style_cache[[hex]] <- s
      }
      s
    }
    legend_colors <- function(n) to_hex(seq(0, 1, length.out = n))
    base_bg <- make_ansi_style("#121212", bg = TRUE)
    flag_style_missing <- combine_ansi_styles(cli::bg_br_white, cli::col_black)
    flag_style_r2 <- combine_ansi_styles(style_bold, style_underline, cli::style_italic, cli::col_green)
    flag_style_both <- combine_ansi_styles(flag_style_missing, flag_style_r2, make_ansi_style("#004200"))
    style_na <- combine_ansi_styles(cli::bg_br_yellow, cli::style_blurred, make_ansi_style("#0000FF"))
  } else {
    style_for_value <- function(t) function(s) cli::bg_black(cli::col_white(s))
    legend_colors <- function(n) rep("black", n)
    flag_style_r2 <- combine_ansi_styles(cli::bg_blue, cli::col_white)
    flag_style_missing <- combine_ansi_styles(cli::bg_yellow, cli::col_black)
    flag_style_both <- combine_ansi_styles(cli::bg_magenta, cli::col_white)
    style_na <- combine_ansi_styles(cli::bg_red, cli::col_blue)
    base_bg <- cli::bg_black
  }

  pad_cell_parts <- function(s) {
    ws <- cell_width - nchar(s)
    list(left = strrep(" ", max(ws, 0)), text = s)
  }

  # ---- Print header ----
  header <- combine_ansi_styles("bold", "underline")
  cat(header(sprintf("%-*s", row_lab_width, "")), sep = "")
  header_cells <- vapply(cn, function(cc) {
    parts <- pad_cell_parts(cc)
    paste0(parts$left, parts$text)
  }, character(1))
  cat(header(paste(header_cells, collapse = " ")), sep = "")
  if (!is.null(extra_col) && !is.null(extra_col_label)) {
    cat(header(paste0("   |   ", sprintf("%-*s", extra_width, extra_col_label))))
  }
  cat("\n")

  # ---- Print body ----
  for (i in seq_len(nr)) {
    row_lab <- sprintf(paste0("%-", row_lab_width, "s"), rn[i])
    cat(style_bold(row_lab), sep = "")
    if (show_border) cat("| ")
    line_parts <- character(nc)
    for (j in seq_len(nc)) {
      val <- x[i, j]
      txt <- formatted[i, j]
      parts <- pad_cell_parts(txt)
      if (is.na(val)) {
        styled_text <- style_na(parts$text)
      } else {
        t <- norm01(val)
        styled_text <- if (has_256_plus) base_bg(style_for_value(t)(parts$text)) else style_for_value(t)(parts$text)
        if (f1[i, j] && f2[i, j]) {
          styled_text <- flag_style_both(parts$text)
        } else if (f1[i, j]) {
          styled_text <- flag_style_r2(parts$text)
        } else if (f2[i, j]) {
          styled_text <- flag_style_missing(parts$text)
        }
      }
      line_parts[j] <- paste0(parts$left, styled_text)
    }
    cat(paste(line_parts, collapse = " "), sep = "")
    if (!is.null(extra_col)) {
      cat(" | ", sprintf("%-*s", extra_width, extra_col[i]), sep = "")
    }
    cat("\n")
  }

  # ---- Print legend ----
  if (legend) {
    legend_char <- "\u25A0" # blacksquare symbol

    blocks <- if (has_256_plus) {
      vapply(legend_colors(legend_cols), function(col) {
        base_bg(make_ansi_style(col, colors = ncols_supported)(legend_char))
      }, character(1))
    } else {
      vapply(rep("black", legend_cols), function(col) {
        cli::bg_black(cli::col_white(legend_char))
      }, character(1))
    }

    # Format min and max values
    min_label <- formatC(vmin, format = "f", digits = digits)
    max_label <- formatC(vmax, format = "f", digits = digits)

    # Compute spacing for labels
    block_line <- paste(blocks, collapse = "")
    left_space <- nchar(min_label) + 1
    right_space <- nchar(max_label)
    total_width <- left_space + legend_cols + right_space
    title_space <- strrep(" ", floor(total_width/2)-3)

    cat("\n", title_space, header("Legend"), title_space, "\n")

    # Position "min" and "max" labels above their respective values
    min_pos <- 1
    max_pos <- total_width - nchar("max")
    label_line <- paste0(
      strrep(" ", min_pos ),
      cli::style_italic(style_underline("min")),
      strrep(" ", max_pos - min_pos - nchar("min")),
      cli::style_italic(style_underline("max"))
    )

    # Print value scale with labels and blocks

    cat("  ", label_line, "\n", sep = "")
    cat("  ", min_label, " ", block_line, " ", max_label, "\n", sep = "")

    if (!has_256_plus) cat("  (no colour scaling in low-colour mode)\n")

    cat("  ", header("Flags:"), "  ",
        sprintf("%s \u2192 %s", r2_label, flag_style_r2("value")), "; ",
        sprintf("%s \u2192 %s", missing_label, flag_style_missing("value")), "; ",
        sprintf("both \u2192 %s", flag_style_both("value")), "\n", sep = "")

    if (show_flag_counts) {
      cat(style_bold("  Counts: "),
          sprintf("%s = %d,     %s = %d,     both = %d",
                  r2_label, sum(f1, na.rm = TRUE),
                  missing_label, sum(f2, na.rm = TRUE),
                  sum(f1 & f2, na.rm = TRUE)),
          "\n", sep = "")
    }
  }
}

#' Pretty-print \eqn{MO_2} matrix for an AquaResp experiment
#'
#' Reads \eqn{MO_2} data from an AquaResp experiment directory and prints a formatted matrix
#' with optional flags and start times for each cycle.
#'
#' @param path Character. Directory location of the AquaResp experiment on disk.
#' @param chamber_col Character. Column name for chamber IDs (default \code{"chamber"}).
#' @param cycle_col Character. Column name for cycle IDs (default \code{"cycle"}).
#' @param mo2_col Character. Column name for \eqn{MO_2} values (default \code{"MO2"}).
#' @param r2_col Character. Column name for \eqn{R^2} values (default \code{"R.2"}).
#' @param min_po2_col Character. Column name for minimum \eqn{PO_2} values (default \code{"minimum.po2"}).
#' @param max_po2_col Character. Column name for maximum \eqn{PO_2} values (default \code{"max.po2"}).
#' @param r2_threshold Numeric. Threshold for \eqn{R^2} flagging (default \code{0.95}).
#' @param ... Additional arguments passed to \code{\link{print_matrix_exp}}.
#'
#' @return Invisibly returns \code{NULL}.
#' @export
#'
#' @details
#' This function:
#' \itemize{
#'   \item Loads AquaResp experiment data using \code{\link{get_exp_MO2s}}.
#'   \item Validates required columns: chamber, cycle, \eqn{MO_2}, and \eqn{R^2}.
#'   \item Checks for duplicate (cycle, chamber) combinations and stops if found.
#'   \item Constructs:
#'     \enumerate{
#'       \item A matrix of \eqn{MO_2} values indexed by cycle and chamber.
#'       \item A matrix of \eqn{R^2} values for flagging.
#'       \item A missingness mask based on \eqn{PO_2} values.
#'     }
#'   \item Computes start times for each cycle and passes them as an extra column.
#'   \item Calls \code{\link{print_matrix_exp}} to render the matrix with:
#'     \itemize{
#'       \item Colour-mapped values.
#'       \item Flags for low \eqn{R^2} and missing \eqn{PO_2}.
#'       \item A legend showing colour scale and flag meanings.
#'       \item An extra column for cycle start times.
#'     }
#' }
#'
#' @section Flags:
#' \itemize{
#'   \item \eqn{R^2} flag: cells where \eqn{R^2 <} \code{r2_threshold}.
#'   \item Missing flag: cells where \eqn{PO_2} values are non-finite or \eqn{\le 0}. For example, when oxygen registers as 0 due to the computer missing a measurement, or when -300 is registered when the probe was out of position, will be flagged as missing values as these will affect \eqn{R^2} and \eqn{MO_2} estimates.
#' }
#'
#' @examples
#' # Example usage:
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' print_exp_mo2s(path = exp_dir_path)
print_exp_mo2s <- function(
    path,
    chamber_col   = "chamber",
    cycle_col     = "cycle",
    mo2_col       = "MO2",
    r2_col        = "R.2",
    min_po2_col   = "minimum.po2",
    max_po2_col   = "max.po2",
    r2_threshold  = 0.95,
    ...
) {
  # ---- Load experiment MO2 data ----
  exp_mo2s <- get_exp_MO2s(path)
  if (!is.data.frame(exp_mo2s)) stop("get_exp_MO2s() did not return a data.frame.")

  # ---- Validate required columns ----
  required_cols <- c(chamber_col, cycle_col, mo2_col, r2_col)
  missing_cols <- setdiff(required_cols, names(exp_mo2s))
  if (length(missing_cols)) {
    stop("Missing required columns in exp_mo2s: ", paste(missing_cols, collapse = ", "))
  }

  # ---- Strict uniqueness check ----
  key <- paste(exp_mo2s[[cycle_col]], exp_mo2s[[chamber_col]], sep = "\r")
  dup_idx <- duplicated(key) | duplicated(key, fromLast = TRUE)
  if (any(dup_idx)) {
    dup_df <- stats::aggregate(
      rep(1L, length(key)),
      by = list(cycle = exp_mo2s[[cycle_col]], chamber = exp_mo2s[[chamber_col]]),
      FUN = length
    )
    dup_df <- dup_df[dup_df$x > 1L, , drop = FALSE]
    sample_rows <- utils::head(dup_df, 10L)
    msg <- paste0(
      "Found duplicate rows for (cycle, chamber); strict mode requires unique keys.\n",
      "Examples (up to 10):\n",
      paste(sprintf("  cycle=%s, chamber=%s, n=%d",
                    sample_rows$cycle, sample_rows$chamber, sample_rows$x), collapse = "\n"),
      if (nrow(dup_df) > 10L) sprintf("\n ...and %d more.", nrow(dup_df) - 10L) else ""
    )
    stop(msg)
  }

  # ---- Factor levels ----
  cycle_levels <- sort(unique(exp_mo2s[[cycle_col]]))
  chamber_levels <- sort(unique(exp_mo2s[[chamber_col]]))
  f_cycle <- factor(exp_mo2s[[cycle_col]], levels = cycle_levels)
  f_chamber <- factor(exp_mo2s[[chamber_col]], levels = chamber_levels)

  # ---- MO2 matrix ----
  mo2_mat <- tapply(exp_mo2s[[mo2_col]], list(f_cycle, f_chamber),
                    FUN = function(x) if (length(x)) x[1] else NA_real_)
  dimnames(mo2_mat) <-
    list(cycle = paste("Cycle", as.character(cycle_levels)),
         chamber = paste0(" Ch", as.character(chamber_levels))
    )

  # ---- r2 matrix ----
  r2_mat <- tapply(exp_mo2s[[r2_col]], list(f_cycle, f_chamber),
                   FUN = function(x) if (length(x)) x[1] else NA_real_)
  dimnames(r2_mat) <- dimnames(mo2_mat)

  # ---- Missingness from PO2 ----
  has_po2 <- all(c(min_po2_col, max_po2_col) %in% names(exp_mo2s))
  if (has_po2) {
    min_mat <- tapply(exp_mo2s[[min_po2_col]], list(f_cycle, f_chamber),
                      FUN = function(x) if (length(x)) x[1] else NA_real_)
    max_mat <- tapply(exp_mo2s[[max_po2_col]], list(f_cycle, f_chamber),
                      FUN = function(x) if (length(x)) x[1] else NA_real_)
    dimnames(min_mat) <- dimnames(max_mat) <- dimnames(mo2_mat)
    missing_mat <- !(is.finite(min_mat) & is.finite(max_mat) & (min_mat > 0) & (max_mat > 0))
  } else {
    missing_mat <- !(is.finite(mo2_mat) & (mo2_mat > 0))
  }

  # ---- Compute start times for each cycle ----
  cycles <- lapply(cycle_levels, function(c) read_cycle(c, path))
  start_times <- vapply(cycles, .cycle_start_dt, FUN.VALUE = as.POSIXct(NA))
  start_times_fmt <- format(as.POSIXct(start_times), "%Y-%m-%d %H:%M:%S")

  # ---- Call pretty-printer with extra column ----
  print_matrix_exp(
    x            = mo2_mat,
    r2_mat       = r2_mat,
    r2_threshold = r2_threshold,
    missing_mat  = missing_mat,
    extra_col    = start_times_fmt,
    extra_col_label = "  Start Time",
    ...
  )

  invisible(NULL)
}

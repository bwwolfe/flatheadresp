#' Plot MO2 values over time for an AquaResp experiment
#'
#' Generates a ggplot visualization of \eqn{MO_{2}} values over time for each chamber in an experiment.
#' Points are colored and shaped by chamber, with cycle numbers optionally annotated above the highest
#' MO2 in each cycle.
#'
#' @param path Character. Experiment file path on disk.
#' @param chambers Integer vector. Optional selection of chambers:
#'   - Positive values = include only these chambers.
#'   - Negative values = exclude these chambers.
#'   - NULL = include all chambers.
#' @param time_col Character. Column name for time values (default = "TIME.UNIX").
#' @param mo2_col Character. Column name for MO2 values (default = "MO2").
#' @param chamber_col Character. Column name for chamber identifiers (default = "chamber").
#' @param show_cycle_labels Logical. If TRUE (default), annotate cycle numbers at the max MO2 per cycle.
#'
#' @return A ggplot object showing MO2 vs time for selected chambers.
#' @details The output is a ggplot object, so additional layers can be added (see examples).
#' @export
#' @importFrom ggplot2 ggplot theme aes
#' @importFrom rlang .data
#' @examples
#' library(ggplot2)
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' plot_exp(path = exp_dir_path)
#'
#' # include chambers 1 and 3, change y axis range
#' # alternatively can use + ylim(0,1200) but this will remove data beyond these limits (e.g.
#' # lines connecting dots beyond this MO2 range would not appear correctly
#' plot_exp(path = exp_dir_path, chambers = c(1, 3)) + coord_cartesian(ylim = c(0,1200))
#'
#' # exclude chamber 4 (the blank (for background resp) with extreme values)
#' plot_exp(path = exp_dir_path, chambers = c(-4))
#'
#' # remove the cycle labels
#' plot_exp(path = exp_dir_path, show_cycle_labels = FALSE)
plot_exp <- function(path,
                     chambers = NULL,
                     time_col = "TIME.UNIX",
                     mo2_col = "MO2",
                     chamber_col = "chamber",
                     show_cycle_labels = TRUE){

  # Generate data from path
  data <- get_exp_MO2s(path)
  data[[time_col]] <- as.POSIXct(data[[time_col]])
  data[[chamber_col]] <- as.integer(data[[chamber_col]])

  # Validate columns
  required_cols <- c(time_col, mo2_col, chamber_col)
  if (!all(required_cols %in% names(data))) {
    stop("Data must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Resolve chamber selection
  all_chambers <- as.integer(unique(data[[chamber_col]]))
  sel_chambers <- .resolve_chambers(all_chambers, chambers)
  data <- data[data[[chamber_col]] %in% sel_chambers, ]
  data[[chamber_col]] <- as.factor(data[[chamber_col]])

  # Build base plot
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x     = .data[[time_col]],
      y     = .data[[mo2_col]],
      color = .data[[chamber_col]],
      shape = .data[[chamber_col]],
      group = .data[[chamber_col]]
    )
  ) +
    ggplot2::geom_line(size = 0.4) +
    ggplot2::geom_point(size = 3) +
    cowplot::theme_cowplot(16) +
    ggplot2::labs(color = "Chamber", shape = "Chamber") +
    ggplot2::scale_x_datetime(name = "Time") +
    ggplot2::scale_y_continuous(
      name = expression(MO[2]~(mg~O[2]~kg^-1~h^-1)),
      n.breaks = 10
    ) +
    ggplot2::theme(
      legend.position = "right",
      legend.box.background = ggplot2::element_rect(colour = "black"),
      legend.box.margin = ggplot2::margin(6, 6, 6, 8)
    ) +
    ggplot2::scale_colour_manual(
      values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728")
    )

  # Optionally add cycle labels
  if (isTRUE(show_cycle_labels)) {
    if (!("cycle" %in% names(data))) {
      warning("`cycle` column not found; cycle labels skipped.")
    } else {
      label_data <- data |>
        dplyr::group_by(cycle) |>
        dplyr::slice_max(order_by = .data[[mo2_col]], n = 1, with_ties = FALSE) |>
        dplyr::ungroup()

      p <- p +
        ggplot2::geom_label(
          data = label_data,
          ggplot2::aes(label = cycle),
          nudge_y = 0.05 * max(data[[mo2_col]], na.rm = TRUE),
          size = 4,
          color = "black"
        )
    }
  }

  return(p)
}

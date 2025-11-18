#' Read a single AquaResp cycle file
#'
#' AquaResp stores each experimental cycle (i.e., the measurement period when MO₂ is calculated)
#' as a separate file. This function reads the data from a cycle file in the "All slopes" folder
#' of an AquaResp experiment.
#'
#' @param cycle_number Numeric. The cycle number to read.
#' @param path Character. The path to the AquaResp experiment directory.
#'
#' @return A dataframe containing the data from the specified cycle.
#' @export
#'
#' @examples
#' # Locate the example AquaResp experiment directory shipped with the package
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' # Read cycle 1 from the example experiment
#' read_cycle(cycle_number = 1, path = exp_dir_path)

read_cycle <- function(cycle_number, path) {
  file_path <-
    file.path(path, "All slopes", paste0("Cycle_", cycle_number, ".txt"))

  if (!file.exists(file_path)) {
    warning(paste("Cycle file", file_path, "not found. Returning NULL."))
    return(NULL)
  }

  res <- utils::read.table(file_path, sep = ";", header = TRUE)

  if ("X" %in% names(res) && all(is.na(res$X))) {
    res <- res[, setdiff(names(res), "X"), drop = FALSE]
  }

  return(res)
}

#' Calculate MO2 values for a specific AquaResp cycle
#'
#' Calculates mass-specific oxygen consumption (MO₂) for each chamber in a given cycle using linear regression on PO₂ data and experimental metadata.
#'
#' @param cycle_number Numeric. The cycle number to calculate MO₂ for.
#' @param path Character. The path to the AquaResp experiment directory.
#'
#' @return A named numeric vector of MO₂ values for each chamber.
#' @export
#'
#' @examples
#' # Example data directory shipped with the package
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' calc_cycle_mo2s(cycle_number = 1, path = exp_dir_path)
calc_cycle_mo2s <-
  function(cycle_number, path) {

    chambers <- get_chambers(path)

    cycle <- read_cycle(cycle_number, path)

    po2_cols <-
      grep("po2", names(cycle))

    cycle_time_x <- cycle$Unix.Time - min(cycle$Unix.Time) + 1

    cycle_slope_list <-
      lapply(po2_cols, function(i) {
        slope_vec <-
          stats::lm(cycle[[i]] ~ cycle_time_x) |>
          summary() |>
          (\(x)x$coefficients[2,])()
        return(slope_vec)
      }
      )

    meta <- get_exp_metadata(path)

    beta <- meta$`Oxygen.solubilty..mg.O2...L`

    rRespFish <-
      meta$`Real.volume..vresp...vfish...neutrally.bouyant...L` /
      meta$`Mass.of.fish..kg`

    sapply(chambers, \(i){
      slope <- cycle_slope_list[[i]][1]
      MO2 <- -1 * (slope / 100) * beta[i] * rRespFish[i] * 3600.0
      return(MO2)
    }) |> stats::setNames(chambers)
  }


#' Get min and max PO2 values for each chamber in a cycle
#'
#' Extracts the minimum and maximum non-zero PO₂ values for each chamber from a cycle dataframe.
#'
#' @param cycle Dataframe. A cycle dataframe as returned by `read_cycle()`.
#'
#' @return A named numeric vector of min and max PO₂ values for each chamber.
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' cycle <- read_cycle(1, exp_dir_path)
#' get_cycle_min_max(cycle)
get_cycle_min_max <-
  function(cycle) {
    po2_cols <-
      grep("po2", names(cycle))

    unlist(lapply(cycle[po2_cols], \(x) {
      c(min = min(x[x != 0]), max = max(x[x != 0]))
    }))
  }


#' Calculate percentage of zero PO2 values in a cycle
#'
#' Computes the percentage of zero values in each PO₂ column of a cycle dataframe.
#' Values greater than zero would have an incorrect R² as calculated by AquaResp.
#'
#' @param cycle Dataframe. A cycle dataframe as returned by `read_cycle()`.
#'
#' @return A named list of percentages (0–100) indicating missingness per chamber.
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' cycle <- read_cycle(1, exp_dir_path)
#' get_cycle_missingness(cycle)
get_cycle_missingness <-
  function(cycle) {
    po2_cols <-
      grep("po2", names(cycle))

    lapply(cycle[po2_cols], function(x)
      round(length(which(x == 0)) / length(x), 4) * 100) |>
      stats::setNames(paste0(names(cycle)[po2_cols], ".pct0"))
  }


#' Calculate correlation between PO2 and time for each chamber
#'
#' Performs Pearson correlation (r) tests between PO₂ values and Unix time for each chamber in a cycle, retyrns r ^ 2
#'
#' @param cycle Dataframe. A cycle dataframe as returned by `read_cycle()`.
#'
#' @return A named numeric vector of correlation coefficients and p-values for each chamber.
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' cycle <- read_cycle(1, exp_dir_path)
#' get_cycle_R2s(cycle)
get_cycle_R2s <-
  function(cycle) {
    po2_cols <-
      grep("po2", names(cycle))

    cors <-
      lapply(cycle[po2_cols], function(x, secs) {
        # was using 2nd column "Seconds..." renamed 'secs' above but at midnight
        # there's a bug
        secs.no0 <- secs[x != 0]
        x.no0 <- x[x != 0]
        cor_result <-
          stats::cor.test(x = x.no0, y = secs.no0, method = "pearson")
      }, secs = cycle$Unix.Time)

    r2s <- lapply(cors, function(x)
      x$estimate^2) |> #returns r ^ 2
      stats::setNames(paste0(names(cycle[po2_cols]), ".r2"))
    return(r2s)
  }


#' Summarize PO2 metrics across all cycles in an experiment
#'
#' Generates a summary dataframe with PO₂ statistics (min, max, correlation, missingness)
#' for each cycle and chamber in the experiment. The correlations are calculated with
#' missing values (PO₂ = 0) omitted to get around the bug in AquaResp.
#'
#' @param path Character. The path to the AquaResp experiment directory.
#'
#' @return A dataframe with one row per cycle and columns for each PO₂ metric for each chamber.
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' get_exp_cycle_summary(path = exp_dir_path)
get_exp_cycle_summary <- function(path){

  cycles <-
    gsub("Cycle_|.txt", "", list.files(file.path(path, "All slopes"))) |>
    as.numeric() |> sort()

  sapply(cycles, function(cycle, path) {
    cycle <- read_cycle(cycle, path)

    cycle_summary <-
      c(
        get_cycle_min_max(cycle),
        get_cycle_R2s(cycle),
        unlist(get_cycle_missingness(cycle))
      )

    return(cycle_summary)
  },
  path = path) |> t() |> as.data.frame()
}


#' Plot PO2 values for a specific cycle
#'
#' Generates a scatter plot of PO₂ values over time for each chamber in a given cycle.
#'
#' @param cycle_number Numeric. The cycle number to plot.
#' @param path Character. The path to the AquaResp experiment directory.
#'
#' @return A base-R plot of time vs PO₂.
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' plot_cycle_po2(cycle_number = 1, path = exp_dir_path)
plot_cycle_po2 <- function(cycle_number, path) {

  the_cycle <- read_cycle(cycle_number, path)
  color_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728")
  po2cols <- grep("^ch[0-9]+\\.po2$", names(the_cycle))
  plot(the_cycle$Unix.Time, the_cycle[[po2cols[1]]], col = color_palette[1], ylim = c(0,110))
  for (i in 2:length(po2cols)) {
    graphics::points(the_cycle$Unix.Time, the_cycle[[po2cols[i]]], col = color_palette[i], ylim = c(0,110))
  }
}

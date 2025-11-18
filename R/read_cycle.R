#' Read a single AquaResp cycle file
#'
#' AquaResp stores each experimental cycle (i.e., the measurement period when \eqn{MO_{2}} is calculated)
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
#' Calculates mass-specific oxygen consumption (\eqn{MO_{2}}) for each chamber in a given cycle using linear regression on \eqn{PO_{2}} data and experimental metadata.
#'
#' @param cycle_number Numeric. The cycle number to calculate \eqn{MO_{2}} for.
#' @param path Character. The path to the AquaResp experiment directory.
#'
#' @return A named numeric vector of \eqn{MO_{2}} values for each chamber.
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
#' Extracts the minimum and maximum non-zero \eqn{PO_{2}} values for each chamber from a cycle dataframe.
#'
#' @param cycle Dataframe. A cycle dataframe as returned by `read_cycle()`.
#'
#' @return A named numeric vector of min and max \eqn{PO_{2}} values for each chamber.
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
#' Computes the percentage of zero values in each \eqn{PO_{2}} column of a cycle dataframe.
#' Values greater than zero would have an incorrect \eqn{R^{2}} as calculated by AquaResp.
#'
#' @param cycle Dataframe. A cycle dataframe as returned by `read_cycle()`.
#'
#' @encoding UTF-8
#'
#' @return A named vector of percentages (0 to 100%) indicating missingness per chamber.
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
      stats::setNames(gsub(".po2", "", paste0(names(cycle)[po2_cols], ".pct0"))) |>
      unlist()
  }

#' Calculate correlation between PO2 and time for each chamber
#'
#' Performs Pearson correlation (r) tests between \eqn{PO_{2}} values and Unix time for each chamber in a cycle, retyrns r ^ 2
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

    r2s <- lapply(cors, function(x) {
      x$estimate^2 |> unname()
      }) |> #returns r ^ 2
      stats::setNames(gsub(".po2","",paste0(names(cycle[po2_cols]), ".r2")))

    r2_pvals <- lapply(cors, function(x) {
    x$p.value
    }) |>
      stats::setNames(gsub(".po2","",paste0(names(cycle[po2_cols]), ".p")))

  unlist(c(r2s, r2_pvals))
  }


#' Summarize PO2 metrics across all cycles in an experiment
#'
#' Generates a summary dataframe with \eqn{PO_{2}} statistics (min, max, correlation, missingness)
#' for each cycle and chamber in the experiment. The correlations are calculated with
#' missing values (\eqn{PO_{2}} = 0) omitted to get around the bug in AquaResp.
#'
#' @encoding UTF-8
#'
#' @param path Character. The path to the AquaResp experiment directory.
#'
#' @return A dataframe with one row per cycle and columns for each \eqn{PO_{2}} metric for each chamber.
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
        get_cycle_missingness(cycle)
      )

    return(cycle_summary)
  },
  path = path) |> t() |> as.data.frame()
}

#' Pivot output from get_exp_cycle_summary() to long format
#'
#' This is a helper function
#'
#' @encoding UTF-8
#'
#' @param df Dataframe. Output from get_exp_cycle_summary()
#'
#' @return long format summary with columns for chamber and cycle and metrics - PO_2 minimum, PO_2 maximum, cycle correlation R^2, R^2 p value, and % missing PO_2 data.
#' @export
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' x <- get_exp_cycle_summary(path = exp_dir_path)
#' y <- pivot_cycle_summary_long(x)
#' rbind(head(y),tail(y))
pivot_cycle_summary_long <- function(df) {
  df$cycle <- 1:nrow(df)

  long_df <-df |>
    tidyr::pivot_longer(cols = -cycle,
                        names_to = c("chamber", "metric"),
                        names_pattern = "ch(\\d+)\\.(.+)") |>
    dplyr::mutate(chamber = as.integer(chamber)) |>
    tidyr::pivot_wider(names_from = metric, values_from = value)
  return(long_df)
}

#' Summarise experiment metadata
#'
#' Prints a summary of key experimental details including chambers, fish masses,
#' respirometer volumes, cycle timing, and environmental conditions.
#'
#' @param path Character. The path to the AquaResp experiment directory.
#'
#' @return Invisibly returns NULL. Called for its side effect of printing to the console.
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                              package = "flatheadresp")
#' summarise_experiment(path = exp_dir_path)
summarise_experiment <- function(path) {
  chambers <- get_chambers(path)
  meta <- get_exp_metadata(path)

  cat("Experiment started", as.character(meta$exp_start[1]), "\n\n")
  cat(length(chambers), "chambers\n\n")

  # Fish masses
  cat("Mass of fish by chamber:\n")
  for (i in seq_along(chambers)) {
    chamber <- chambers[i]
    mass <- meta$`Mass.of.fish..kg`[i]
    cat("  Chamber", chamber, ":", mass, "kg\n")
  }

  # Respirometer volumes
  volumes <- meta$`Volume.respirometer..L`
  if (length(unique(volumes)) == 1) {
    cat("\nRespirometer volume is", volumes[1], "L for all chambers.\n")
  } else {
    cat("\nRespirometer volume by chamber:\n")
    for (i in seq_along(chambers)) {
      chamber <- chambers[i]
      vol <- volumes[i]
      cat("  Chamber", chamber, ":", vol, "L\n")
    }
  }

  # Cycle timing
  flush_time <- meta$`Flush.time..s`[1]
  wait_time <- meta$`Wait.time..s`[1]
  measure_time <- meta$`Measurement.time..s`[1]

  cat("\nCycle timing:\n")
  cat("  Flush time:", flush_time, "seconds\n")
  cat("  Wait time:", wait_time, "seconds\n")
  cat("  Measurement time:", measure_time, "seconds\n")

  # Salinity and temperature
  salinity <- meta$Salinity[1]
  temperature <- meta$Temperature[1]

  cat("\nEnvironmental conditions:\n")
  cat("  Salinity:", salinity, "ppt\n")
  cat("  Temperature:", temperature, "\u00B0C\n")

  # data QA/QC metrics
  cycle_summary <- get_exp_cycle_summary(path) |> pivot_cycle_summary_long()

}

#' Plot PO2 values for a specific cycle
#'
#' Generates a scatter plot of \eqn{PO_{2}} values over time for each chamber in a given cycle.
#'
#' @param cycle_number Numeric. The cycle number to plot.
#' @param path Character. The path to the AquaResp experiment directory.
#' @param ylim Numeric vector. The y axis plot limits, if NULL (default) it will be the range of po2 values across the cycle. May need to change to c(0, 100) or e.g. c(80,100) to zoom enough to see if there were any 0s or missing data (-300) during the cycle.
#'
#' @return A base-R plot of time vs \eqn{PO_{2}}.
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' plot_cycle_po2(cycle_number = 1, path = exp_dir_path)

plot_cycle_po2 <- function(cycle_number, path, ylim = NULL) {

  the_cycle <- read_cycle(cycle_number, path)
  color_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728")
  po2cols <- grep("^ch[0-9]+\\.po2$", names(the_cycle))
  if(is.null(ylim))
    ylim <- range(the_cycle[po2cols]) +
              c(-0.05, 0.05) * diff(range(the_cycle[po2cols]))
  plot(the_cycle$Unix.Time, the_cycle[[po2cols[1]]],
       col = color_palette[1], ylim = ylim)
  for (i in 2:length(po2cols)) {
    graphics::points(the_cycle$Unix.Time, the_cycle[[po2cols[i]]],
                     col = color_palette[i])
  }
}

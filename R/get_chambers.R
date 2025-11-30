#' Get the active chambers from an AquaResp experiment directory
#'
#' This function scans the experiment directory and returns the chamber numbers
#' for which summary data files are available. These correspond to the active chambers
#' used during the AquaResp experiment.
#'
#' @param path Character. The directory location of the AquaResp experiment on disk.
#' This should be the folder created by AquaResp that contains folders like "All slopes",
#' "Experimental information", and "Oxygen data raw".
#'
#' @return A numeric vector of active chamber numbers.
#' @export
#'
#' @examples
#' # Locate the example AquaResp experiment directory shipped with the package
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' print(exp_dir_path)
#' get_chambers(path = exp_dir_path)
get_chambers <- function(path) {
  gsub("Summary data resp |.txt",
       "",
       list.files(path, pattern = "Summary data resp .*\\.txt")) |>
    as.numeric() |>
    sort()

}

#' Get metadata for a specific chamber from an AquaResp experiment
#'
#' This function reads the first lines of the summary data file for a given chamber
#' and extracts experimental metadata such as fish mass, respirometer volume, and start time.
#'
#' @param path Character. The directory location of the AquaResp experiment on disk.
#' @param chamber Numeric. The chamber number to retrieve metadata for. Use `get_chambers()` to list available chambers.
#'
#' @return A named list containing metadata for the specified chamber.
#' @export
#'
#' @examples
#' # Locate the example AquaResp experiment directory shipped with the package
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' print(exp_dir_path)
#'
#' get_chamber_metadata(path = exp_dir_path, chamber = 1)
get_chamber_metadata <- function(path, chamber) {

  sum_metadata <- list()

  sum_path <- file.path(path, paste0("Summary data resp ", chamber, ".txt"))

  con <- file(sum_path, "r")
  for (j in 1:15) {
    line <- readLines(con, n = 1)
    if (grepl(":", line)) {
      parts <- strsplit(line, ":")[[1]]
      key <- trimws(parts[1])
      value <- trimws(parts[2]) |> as.numeric() |> round(4)
      sum_metadata[[key]] <- value
    }
  }
  close(con)

  sum_metadata$exp_start <-
    sum_metadata$`Experiment start, UNIX time` |> as.POSIXct()

  sum_metadata <- c(list(chamber = chamber), sum_metadata)
  return(sum_metadata)
}

#' Get metadata for all chambers in an AquaResp experiment
#'
#' This function retrieves metadata for all active chambers in the experiment directory
#' and returns it as a dataframe, with one row per chamber.
#'
#' @param path Character. The directory location of the AquaResp experiment on disk.
#'
#' @return A dataframe where each row contains metadata for one chamber.
#' @export
#'
#' @examples
#' # Locate the example AquaResp experiment directory shipped with the package
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' print(exp_dir_path)
#' get_exp_metadata(path = exp_dir_path)

get_exp_metadata <- function(path) {
  chambers <- get_chambers(path)

  metadata_list <-
    lapply(chambers, function(i) {
      get_chamber_metadata(path, chamber = i)

    })

  do.call(rbind, lapply(metadata_list, as.data.frame))
}


#' Get \eqn{MO_{2}} values from an AquaResp experiment
#'
#' This function reads the summary data files and returns mass-specific oxygen consumption
#' (\eqn{MO_{2}}) values for one or more chambers. If no chambers are specified, data from all
#' chambers is returned.
#'
#' @param path Character. The directory location of the AquaResp experiment on disk.
#' @param chambers Optional integer vector. Chamber numbers to include or exclude.
#'   Use negative values to exclude chambers (e.g., \code{c(-2, -3)} excludes chambers 2 and 3).
#'
#' @return A \code{data.frame} with columns for chamber number, cycle number, and \eqn{MO_{2}} values.
#' @export
#'
#' @examples
#' # Locate the example AquaResp experiment directory shipped with the package
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment", package = "flatheadresp")
#' print(exp_dir_path)
#' get_exp_mo2s(path = exp_dir_path) # All chambers
#' get_exp_mo2s(path = exp_dir_path, chambers = 1) # Chamber 1 only
#' get_exp_mo2s(path = exp_dir_path, chambers = c(-2, -3)) # Exclude chambers 2 and 3

get_exp_mo2s <- function(path, chambers = NULL) {

  # Get all available chambers
  all_chambers <- get_chambers(path)

  # Resolve selection
  if (is.null(chambers)) {
    sel_chambers <- all_chambers
  } else {
    chambers <- as.integer(chambers)

    # Validate requested chambers
    missing <- setdiff(abs(chambers), all_chambers)
    if (length(missing) > 0) {
      cli::cli_abort("Invalid chamber(s): {missing}. Available: {all_chambers}")
    }

    # Apply inclusion/exclusion logic
    if (any(chambers < 0)) {
      exclude <- abs(chambers[chambers < 0])
      sel_chambers <- setdiff(all_chambers, exclude)
    } else {
      sel_chambers <- chambers
    }
  }

  # Read summary data for selected chambers
  summary_data_list <- lapply(sel_chambers, function(i) {
    sum_path <- file.path(path, paste0("Summary data resp ", i, ".txt"))

    sum_data <- utils::read.table(sum_path,
                                  skip = 15,
                                  sep = ";",
                                  header = TRUE) |>
      (\(res) res[, setdiff(names(res), "X"), drop = FALSE])()

    data.frame(chamber = i, cycle = seq_len(nrow(sum_data)), sum_data)
  })

  # Combine and reorder columns
  summary_dataframe <- do.call("rbind", summary_data_list)
  summary_dataframe <- summary_dataframe[, c("cycle", "chamber", setdiff(names(summary_dataframe), c("cycle", "chamber")))]

  # Sort by cycle then chamber
  summary_dataframe[order(summary_dataframe$cycle, summary_dataframe$chamber), ]
}

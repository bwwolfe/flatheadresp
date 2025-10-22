
#' Adjust MO2 values using a new fish mass
#'
#' This function recalculates mass-specific oxygen consumption (MO2) values from an AquaResp experiment
#' using a new fish mass value. Like the values calculated by AquaResp, it assumes the fish has neutral buoyancy (i.e., density = 1 g/mL),
#' and corrects the MO2 values based on the change in fish mass and the resulting change in effective
#' respirometer volume.
#'
#' @param path Character. The directory location of the AquaResp experiment on disk.
#' This should be the folder created by AquaResp that contains the summary files.
#' @param chamber Numeric. The chamber number to correct MO2 values for.
#' @param new_mass Numeric. The corrected fish mass in kilograms (kg).
#'
#' @return A numeric vector of corrected MO2 values for the specified chamber.
#' @export
#'
#' @examples
#' \dontrun{
#' # Correct MO2 values for chamber 2 with a new fish mass of 0.345 kg and a new buoyancy of 0.91 g/
#' corrected_MO2s <-
#'    fix_exp_MO2s(path = "~/data/Resp/Experiment_1", chamber = 2, new_mass = 0.345)
#' }

fix_exp_MO2s <- function(path, chamber, new_mass) {

  i <- chamber

  sum_path <- file.path(path, paste0("Summary data resp ", i, ".txt"))

  sum_data <-
    utils::read.table(sum_path,
               skip = 15,
               sep = ";",
               header = TRUE)

  sum_metadata <- list()

  con <- file(sum_path, "r")
  for (i in 1:15) {
    line <- readLines(con, n = 1)
    if (grepl(":", line)) {
      parts <- strsplit(line, ":")[[1]]
      key <- trimws(parts[1])
      value <- trimws(parts[2]) |> as.numeric() |> round(4)
      sum_metadata[[key]] <- value
    }
  }
  close(con)

  mfish_corrected <- new_mass
  MO2 <- sum_data$MO2
  mfish <- sum_metadata$`Mass of fish, kg`
  vresp <- sum_metadata$`Volume respirometer, L`
  vreal <- vresp - mfish
  rRespFish <- vreal / mfish

  # below is cribbed from the aquaresp python code in AquaAnalyse.py
  # from:
  # https://github.com/bigb8/AquaResp/blob/master/Aquaresp%203/lib/AquaAnalyse.py
  #  pO2max, pO2maxkpa = ao.partialpressureoxygen(temp, patm, "mmhg")
  #  oxysolmmhg,oxysolkpa = ao.oxygensolubility(temp,salinity)
  #  vreal = vresp - mfish
  #  rRespFish = float(vreal)/mfish
  #  MO2 = -1*(slope/100)*beta*rRespFish*3600.0
  #  MO2_TOT = -1*(slope/100)*beta*vresp*3600.0
  #These are equal - for fixing bad fish masses
  #cbind(sum_data$MO2 / rRespFish * vresp, sum_abs_data$MO2)

  vreal_corrected <- vresp - mfish_corrected

  rRespFish_corrected <- vreal_corrected / mfish_corrected

  MO2_corrected <- MO2 / rRespFish * rRespFish_corrected
  cat('\n\t',paste("fish mass:", mfish, "->", new_mass))
  cat('\n\t',paste("resp real volume:", vreal, "->", vreal_corrected), '\n\n')

  return(MO2_corrected)
}

#' Test version - Adjust the MO2 data from an AquaResp chamber with a new animal mass value and density
#'
#' This function recalculates mass-specific oxygen consumption (MO2) values from an AquaResp experiment
#' using a new fish mass and an optional density value. It assumes the fish is neutrally buoyant unless
#' a specific density (in g/mL) is provided. The function uses existing metadata and MO2 data from the
#' experiment directory to perform the correction.
#'
#' @param path Character. The directory location of the AquaResp experiment on disk. This should be the folder created by AquaResp that contains the summary files.
#' @param chamber Numeric. The chamber number to correct MO2 values for. Use `get_chambers()` to list available chambers.
#' @param new_mass Numeric. The corrected fish mass in kilograms (kg).
#' @param density Numeric. The density of the fish in g/mL. Defaults to 1.0 (neutral buoyancy).
#'
#' @return A numeric vector of corrected MO2 values for the specified chamber.
#' @export
#'
#' @examples
#' \dontrun{
#' # Correct MO2 values for chamber 2 with a new fish mass of 0.345 kg and density of 1.05 g/mL
#' corrected_MO2s <-
#'  fix_exp_MO2s(path = "~/data/Resp/Experiment_1", chamber = 2,
#'               new_mass = 0.345, density = 1.05)
#' }

fix_exp_MO2s_2 <- function(path, chamber, new_mass, density = 1.0) {

  # Get metadata and MO2 data using existing functions
  sum_metadata <- get_chamber_metadata(path, chamber)
  sum_data <- get_exp_MO2s(path, chamber)

  MO2 <- sum_data$MO2
  mfish <- sum_metadata$`Mass of fish, kg`
  vresp <- sum_metadata$`Volume respirometer, L`

  # Calculate volumes using density
  vfish <- round((mfish * 1000) / density / 1000,4)  # kg -> g -> mL -> L
  vreal <- vresp - vfish
  rRespFish <- vreal / mfish

  vfish_corrected <- round((new_mass * 1000) / density / 1000,4)
  vreal_corrected <- vresp - vfish_corrected
  rRespFish_corrected <- vreal_corrected / new_mass

  MO2_corrected <- MO2 / rRespFish * rRespFish_corrected

  cat('\n\t',paste("fish mass:", mfish, "->", new_mass))
  cat('\n\t',paste("resp real volume:", vreal, "->", vreal_corrected), '\n\n')

  return(MO2_corrected)
}


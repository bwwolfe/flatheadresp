# data-raw/temp_correct_data.R
# Purpose: Build the public dataset 'temp_correct_data' and save to data/my_dataset.rda
# Note: This script is NOT included in the built package (data-raw is .Rbuildignore'd).

## created with usethis::use_data_raw("temp_correct_data")

# water temperature correction for flathead experiments
# author: barrett
# this creates dataframe temp_correct_data called by correct_for_temperature()
# that ships with the package

# load temp logger data

library(ggplot2)
library(plotly)
library(mgcv)

hobo <-
  read.csv("data-raw/Resp_water_jacket_25_apr_2025.csv", skip = 1, row.names = 1) |>
  setNames(c("datetime", "temp", "light")) |>
  transform(datetime = as.POSIXct(datetime, format = "%m/%d/%y %I:%M:%S %p"))

p <- ggplot(hobo) +
  #geom_point(aes(y = light), colour = "blue") +
  geom_line(aes(x = datetime,y = temp),colour = "red4") +
  scale_x_datetime(date_breaks = "day", date_labels = "%b %d") +
  coord_cartesian(ylim = c(15, 19),
                  xlim = as.POSIXct(paste0("2025", c("-04-25", "-05-03")))) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = -1, angle = 45))


ggplotly(p)

# smooth temp logger temp data -
# two reasons - one is its a bit jumpy due to logger resolution.
#   the other is that some trials might be longer than the ones with data
#   so need to extend the series
#   first i did a LOESS but that doesn't allow extrapolating beyond the original data
#   second is GAM which will allow for combining the two trials for which there is un-controlled temp data.

#   trial 12 april 29 exp start 15:40 end 13:55 next day

trial12 <-
  hobo[which(hobo$datetime == as.POSIXct("2025-04-29 15:40:00")):which(hobo$datetime == as.POSIXct("2025-04-30 13:55:00")),]

trial12_ext <-
  hobo[which(hobo$datetime == as.POSIXct("2025-04-29 15:40:00")):which(hobo$datetime == as.POSIXct("2025-04-30 15:00:00")),]

trial13 <-
  hobo[which(hobo$datetime == as.POSIXct("2025-04-30 15:25:00")):which(hobo$datetime == as.POSIXct("2025-05-01 14:25:00")),]

trial13_ext <-
  hobo[which(hobo$datetime == as.POSIXct("2025-04-30 15:30:00")):which(hobo$datetime == as.POSIXct("2025-05-01 16:00:00")),]

# LOESS

loess(temp ~ as.integer(datetime), data = trial12)  |> predict() |> plot(x = trial12$datetime, y = _)
points(temp ~ datetime, trial12, col = "pink4")

# didn't end up using since can't extrapolate longer

#GAM approach

gamdata <-
  rbind(
    transform(trial12,
              trial = "twelve",
              x = as.numeric(datetime-min(datetime))),
    transform(trial13,
              trial = "thirteen",
              x = as.numeric(datetime-min(datetime)))
  ) |>
  transform(trial = factor(trial))

temp_correct_data <-
  data.frame(x = seq(from = 0, by = 300, length.out = 300))

temp_correct_data$est_temp <-
  gam(temp ~ s(x) +
        s(x, trial, bs = "fs"),
      data = gamdata, method = "REML") |>
  predict(exclude = "s(x,trial)",
          newdata = temp_correct_data,
          newdata.guaranteed = T)
temp_correct_data$x <- temp_correct_data$x/60/60
names(temp_correct_data)[1] <- "TIME.HOURS"
plot(temp_correct_data)

usethis::use_data(temp_correct_data, overwrite = TRUE)

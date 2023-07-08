## 2000.csv.bz is as downloaded from
## https://doi.org/10.7910/DVN/HG7NV7/YGU3TD
## on 11 May 2023
## airports.csv is as downloaded from
## https://doi.org/10.7910/DVN/HG7NV7/XTPZZY
## on 11 May 2023

library("dplyr")

data_path <- "~/Repositories/bigbr/manuscript/bigbr_supplementary_material/diverted-flights/data"

air <- read.csv(file.path(data_path, "2000.csv.bz2"))
airports <- read.csv(file.path(data_path, "airports.csv"))

## See
## https://datascience.stackexchange.com/questions/13567/ways-to-deal-with-longitude-latitude-feature
## for lon lat conversion to x, y, z
air <- air |>
    mutate(Month = factor(Month),
           DayOfWeek = factor(DayOfWeek),
           UniqueCarrier = factor(UniqueCarrier)) |>
    left_join(airports[c("long", "lat", "iata")], by = c("Origin" = "iata")) |>
    rename(Orig_lat = lat, Orig_lon = long) |>
    left_join(airports[c("long", "lat", "iata")], by = c("Dest" = "iata")) |>
    rename(Dest_lat = lat, Dest_lon = long) |>
    mutate(Orig_lat = as.numeric(Orig_lat), Orig_lon = as.numeric(Orig_lon),
           Dest_lat = as.numeric(Dest_lat), Dest_lon = as.numeric(Dest_lon),
           Orig_x = cos(Orig_lat) * cos(Orig_lon),
           Orig_y = cos(Orig_lat) * sin(Orig_lon),
           Orig_z = sin(Orig_lat),
           Dest_x = cos(Dest_lat) * cos(Dest_lon),
           Dest_y = cos(Dest_lat) * sin(Dest_lon),
           Dest_z = sin(Dest_lat),
           Delay = ActualElapsedTime - CRSElapsedTime,
           Delay = ifelse(Cancelled == 1 | Diverted == 1, Inf, Delay),
           Delayed = Delay > 0,
           ## Set AQ as the reference level for UniqueCarrier
           UniqueCarrier = relevel(UniqueCarrier, ref = "AQ"))

saveRDS(air, file = file.path(data_path, "air2000_combined.rds"))

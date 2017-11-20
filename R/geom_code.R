
library(readr)
library(lubridate)
library(dplyr)
library(tidyr)
library(purrr)
library(geosphere)
library(ggplot2)
library(grid)
library(ggmap)

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#' Load raw hurricane data from ebtrk_atlc_1988_2015.txt file
#'
#' This function was provided by Coursera to load data from the provided file ebtrk_atlc_1988_2015.txt
#'
#' @return Returns a data frame
#' @export
#'
#' @examples
#' \dontrun{
#'   df_raw <- load_data()
#'   head(df)
#' }
#'
#' @import readr
#'
load_data <- function(){
  # note that this section of the code to read our data was provided on coursera
  ext_tracks_widths <- c(7, 10, 2, 2, 3, 5, 5, 6, 4, 5, 4, 4, 5, 3, 4, 3, 3, 3,
                         4, 3, 3, 3, 4, 3, 3, 3, 2, 6, 1)
  ext_tracks_colnames <- c("storm_id", "storm_name", "month", "day",
                           "hour", "year", "latitude", "longitude",
                           "max_wind", "min_pressure", "rad_max_wind",
                           "eye_diameter", "pressure_1", "pressure_2",
                           paste("radius_34", c("ne", "se", "sw", "nw"), sep = "_"),
                           paste("radius_50", c("ne", "se", "sw", "nw"), sep = "_"),
                           paste("radius_64", c("ne", "se", "sw", "nw"), sep = "_"),
                           "storm_type", "distance_to_land", "final")

  file_path <- file.path("data", "ebtrk_atlc_1988_2015.txt")
  ext_tracks <- read_fwf(file_path,
                         fwf_widths(ext_tracks_widths, ext_tracks_colnames),
                         na = "-99")

  # format the longitude to ensure it is numeric and has negative values for locations in the Western hemisphere
  return(ext_tracks)
}


#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#' Clean raw data collected with load_data function
#'
#' This will take our raw data, filter it for a particular storm observation and
#' modify it to the requested assignment format.
#'
#' @return Returns a data frame
#' @export
#'
#' @param df The raw data frame collected with load_data
#' @param storm The storm_id that we are filtering for
#' @param dt A datetime formatted in lubridate's ymd_h
#' @param westernHemisphere Format longitude to ensure it has negative values for locations in the Western hemisphere
#'
#' @examples
#' \dontrun{
#'   df_raw <- load_data()
#'   katrina_observation <- clean_data(df_raw)
#'   print(katrina_observation)
#'
#'   ike_observation <- clean_data(df_raw, storm="IKE-2008", dt=ymd_h("2008-09-13 6"))
#'   print(ike_observation)
#' }
#'
#' @import readr, dplyr, tidyr, lubridate, purrr
#'
clean_data <- function(df, storm="KATRINA-2005", dt=ymd_h("2005-08-29 12"),
                       westernHemisphere=TRUE){
  # get the storm_id and date columns in the format requested
  df$date <- as.Date(with(df, paste(year, month, day,sep="-")), "%Y-%m-%d")
  df$date <- ymd_h(paste(df$date, df$hour)) # they really want the datetime rather than date
  df$storm_id =paste(df$storm_name, df$year, sep="-")

  # choose our columns of interest
  cols <- c("storm_id", "date", "latitude", "longitude",
            "radius_34_ne", "radius_34_se", "radius_34_sw", "radius_34_nw",
            "radius_50_ne", "radius_50_se", "radius_50_sw", "radius_50_nw",
            "radius_64_ne", "radius_64_se", "radius_64_sw", "radius_64_nw")
  df <- df %>% select(cols)
  df <- df %>% filter(storm_id==storm, date==dt)
  #print(df)
  #print("***********")

  # now we need to break up these columns such that they are in the requested format of
  # wind_speed (knots), ne, nw, se, sw
  df <- df %>% gather(variable, value, -storm_id, -date,-latitude, -longitude, -storm_id) %>%
    mutate(wind_speed = purrr::map_chr(strsplit(variable, "_"), function(x){x[2]})) %>%
    mutate(quadrant = purrr::map_chr(strsplit(variable, "_"), function(x){x[3]})) %>%
    select(-variable) %>% # we don't need this column anymore
    spread(quadrant, value)

  if(westernHemisphere){df$longitude <- df$longitude * -1}
  return(df)
}

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

#' Defines Our GeomHurricane class
#'
#' As per the course materials:
#' The ggproto() function is used to construct a new class corresponding to your new geom.
#' The geom_* function is constructed as a regular function.
#' As a hint, notice that the wind radii geom essentially shows a polygon for each of the wind levels.
#' So we are interested in the polygonGrob.
#' As per coursera, we use the package
#' \href{https://cran.r-project.org/web/packages/geosphere/vignettes/geosphere.pdf}{geosphere}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   geom_hurricane(data = storm_observation,
#'   aes(x = longitude, y = latitude, r_ne = ne, r_se = se, r_nw = nw, r_sw = sw,
#'   fill = wind_speed, colour = wind_speed), scale_radii = scale_radii)
#' }
#'
#' @import ggplot2, geosphere, grid
#'
GeomHurricane <- ggplot2::ggproto("GeomHurricane", Geom,
                         required_aes = c("x", "y", "r_ne", "r_se", "r_nw", "r_sw"),
                         default_aes = aes(alpha = 0.75, fill="red", colour = "red", scale_radii=1),
                         draw_key = draw_key_polygon,
                         draw_group = function(data, panel_scales, coord) {
                           # draw_group for three polygon[GRID.polygon] objects in this case
                           # as per the book... The draw_panel/group function has three arguments to it
                           # The data element is a data frame containing one column for each aesthetic specified
                           # panel_scales is a list containing information about the x and y scales
                           # coord is an object that describes the coordinate system of your plot

                           # https://www.rdocumentation.org/packages/geosphere/versions/1.5-5/topics/destPoint
                           # geosphere docs we have: destPoint(p, b, d, a=6378137, f=1/298.257223563, ...)
                           # p = Longitude and Latitude of point(s), in degrees. Can be a vector of two numbers,
                           # a matrix of 2 columns (first one is longitude, second is latitude)
                           # b = numeric. Bearing (direction) in degrees
                           # d = numeric. Distance in meters (A nautical mile is a unit of length equal to exactly 1,852 meters)

                           #print(data) # data has 1 row (when using draw_group) and columns for aesthetics
                           #     colour  fill     x    y    r_ne r_se r_nw r_sw PANEL group alpha scale_radii
                           #1    red    red    -89.6 29.5  200  200  100  150     1     1     1           1
                           #draw_panel would have returned 3 rows in this case along with columns for aesthetics
                           #print("**********")

                           sections <- data.frame(section = c("r_ne", "r_se", "r_sw", "r_nw"),
                                                  angleFrom = c(0, 90, 180, 270),
                                                  angleTo = c(90, 180, 270, 360))

                           df_all <- data.frame()
                           # now for each section in the each row of data: ne, se, nw, sw
                           for(i in 1:nrow(sections)){
                             p <- c(data$x, data$y)
                             b <- sections[i, "angleFrom"]:sections[i, "angleTo"]
                             direction <- as.character(sections[i, "section"])
                             # nautical mile is a unit of length equal to exactly 1,852 meters.
                             d = data[1, direction] * 1852 * data$scale_radii
                             geoPoints <- geosphere::destPoint(p=p, b=b, d=d)
                             #print(geoPoints) # columns called lon, lat
                             df_geo <- data.frame(x = geoPoints[,"lon"], # collect points provided by geosphere
                                                  y = geoPoints[,"lat"],
                                                  colour = data$colour, # and the data aesthetics
                                                  fill = data$fill,
                                                  PANEL = data$PANEL,
                                                  group = data$group,
                                                  alpha = data$alpha
                                                  )
                             if(nrow(df_all) == 0){
                               df_all <- df_geo
                             }else{
                               df_all <- rbind(df_all, df_geo)
                             }
                           }

                           coords <- coord$transform(df_all, panel_scales)
                           #print(coords)
                           #   x         y         colour   fill PANEL group alpha
                           #1  0.6271408 0.3785001 yellow yellow     1     3     1
                           #2  0.6271435 0.3792029 yellow yellow     1     3     1
                           #3  0.6271585 0.3799057 yellow yellow     1     3     1

                           # Construct a polygon grob
                           pgrob <- polygonGrob(
                             x = coords$x,
                             y = coords$y,                                # needs 'col =' not 'colour ='
                             gp = gpar(col = as.character(coords$colour), # will be black and grey
                                       fill = as.character(coords$fill),  # unless you do explicit conversion to char
                                       alpha = coords$alpha)
                           )
                   }
)

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

#' The geom_hurricane function
#'
#' Used in conjunction with GeomHurricane.
#' See the layout of course materials p.g. 452
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   geom_hurricane(data = storm_observation,
#'   aes(x = longitude, y = latitude, r_ne = ne, r_se = se, r_nw = nw, r_sw = sw,
#'   fill = wind_speed, colour = wind_speed), scale_radii = scale_radii)
#' }
#'
#' @import ggplot2
#'
geom_hurricane <- function(mapping = NULL, data = NULL,
                           stat="identity", position = "identity", na.rm = FALSE,
                           show.legend = NA,
                           inherit.aes = TRUE, ...) {
  ggplot2::layer(
    geom = GeomHurricane,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
#' Map a particular storm observation
#'
#' This will take a storm observation gathered with the clean_data function
#' and map this observation.
#'
#' @param storm_observation A dataframe collected with clean_data representing a storm observation
#' @param location A character string passed to get_map for the location we are looking to map
#' @param zoom The zoom used in ggmap::get_map. Default is 6.
#' @param scale_radii Numeric value which allows us to plot wind radii charts with radii scaled back.
#' @param save_png Boolean value specifying whether to save the resulting map into your working directory. Default is FALSE.
#'
#' @examples
#' \dontrun{
#' df_raw <- load_data()
#' katrina_observation <- clean_data(df_raw) # get the default storm data
#' map_observation(katrina_observation, scale_radii = 1, save_png=TRUE)
#' }
#'
#' @export
#'
map_observation <- function(storm_observation, location="Louisiana", zoom=6, scale_radii=1, save_png=FALSE){
  # chart this data as suggested by coursera
  # with ggmap and testing keep getting: geocode failed with status OVER_QUERY_LIMIT, location = "Louisiana"
  # and have to wait

  map <- get_map(location, zoom=zoom, maptype = "toner-background", source = "stamen") %>%
  ggmap(extent = "device") +
  geom_hurricane(data = storm_observation,
                 aes(x = longitude, y = latitude, r_ne = ne, r_se = se, r_nw = nw, r_sw = sw,
                 fill = wind_speed, colour = wind_speed), scale_radii = scale_radii) +
  scale_color_manual(name = "Wind speed (kts)",
                     values = c("red", "orange", "yellow")) +
  scale_fill_manual(name = "Wind speed (kts)",
                    values = c("red", "orange", "yellow"))

  if(save_png){# save as a png if the file path is specified
    storm_name <- storm_observation[1, "storm_id"]
    png_name <- paste0(storm_name," scale_radii_", scale_radii,".png")
    png(png_name)
    print(map)
    dev.off()
  }
}

#------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------

df_raw <- load_data()
write.csv(df_raw, file = "df_raw.csv") # writing down the raw data

# storm="KATRINA-2005", dt=ymd_h("2005-08-29 12")
katrina_observation <- clean_data(df_raw) # get the default storm data
map_observation(katrina_observation, scale_radii = 1, save_png=TRUE)
map_observation(katrina_observation, scale_radii = 0.5, save_png=TRUE)

# we have plotted Katrina as a test
# now we'll look to plot hurricane Ike
# in the raw data there is data for Ike in 2008 from Sep 1 -> 15
# we are required to use an observation time when the storm was near or over the United States
ike_observation <- clean_data(df_raw, storm="IKE-2008", dt=ymd_h("2008-09-13 6"))
map_observation(ike_observation, location ="Louisiana", zoom=6, scale_radii = 1, save_png=TRUE)
map_observation(ike_observation, location ="Louisiana", zoom=6, scale_radii = 0.5, save_png=TRUE)

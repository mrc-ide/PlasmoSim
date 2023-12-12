#------------------------------------------------
#' @title Import file
#'
#' @description Import file from the inst/extdata folder of this package
#' 
#' @param name name of file.
#'
#' @export

plasmosim_file <- function(name) {
  
  # load file from inst/extdata folder
  name_full <- system.file("extdata/", name, package = 'PlasmoSim', mustWork = TRUE)
  ret <- readRDS(name_full)
  
  # return
  return(ret)
}

# -----------------------------------
# takes matrix as input, converts to list format for use within Rcpp code
#' @noRd
matrix_to_rcpp <- function(x) {
  return(split(x, f = 1:nrow(x)))
}

# -----------------------------------
# takes list format returned from Rcpp and converts to matrix
#' @noRd
rcpp_to_matrix <- function(x) {
  ret <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
  return(ret)
}

#------------------------------------------------
# update progress bar
#' @noRd
update_progress <- function(pb_list, name) {
  knitrProgressBar::update_progress(pb_list$pb_env[[name]])
}

#------------------------------------------------
#' @title Get great circle distance between spatial points
#'
#' @description Get great circle distance between spatial points.
#' 
#' @param lat vector of latitudes.
#' @param lon vector of longitudes.
#'
#' @export

get_spatial_distance <- function(lat, lon) {
  
  # check inputs
  assert_vector(lat)
  assert_numeric(lat)
  assert_vector(lon)
  assert_numeric(lon)
  
  # calculate distance matrix
  ret <- apply(cbind(lat, lon), 1, function(y) {lonlat_to_bearing(lat, lon, y[1], y[2])$gc_dist})
  diag(ret) <- 0
  
  return(ret)
}

#------------------------------------------------
#' @title Calculate great circle distance and bearing between coordinates
#'
#' @description Calculate great circle distance and bearing between spatial
#'   coordinates.
#'
#' @param origin_lon The origin longitude
#' @param origin_lat The origin latitude
#' @param dest_lon The destination longitude
#' @param dest_lat The destination latitude
#'
#' @export
#' @examples
#' # one degree longitude should equal approximately 111km at the equator
#' lonlat_to_bearing(0, 0, 1, 0)

lonlat_to_bearing <- function(origin_lon, origin_lat, dest_lon, dest_lat) {
  
  # check inputs
  assert_vector(origin_lon)
  assert_numeric(origin_lon)
  assert_vector(origin_lat)
  assert_numeric(origin_lat)
  assert_vector(dest_lon)
  assert_numeric(dest_lon)
  assert_vector(dest_lat)
  assert_numeric(dest_lat)
  
  # convert input arguments to radians
  origin_lon <- origin_lon*2*pi/360
  origin_lat <- origin_lat*2*pi/360
  dest_lon <- dest_lon*2*pi/360
  dest_lat <- dest_lat*2*pi/360
  
  # get change in lon
  delta_lon <- dest_lon - origin_lon
  
  # calculate bearing
  bearing <- atan2(sin(delta_lon)*cos(dest_lat), cos(origin_lat)*sin(dest_lat)-sin(origin_lat)*cos(dest_lat)*cos(delta_lon))
  
  # calculate great circle angle. Use temporary variable to avoid acos(>1) or 
  # acos(<0), which can happen due to underflow issues
  tmp <- sin(origin_lat)*sin(dest_lat) + cos(origin_lat)*cos(dest_lat)*cos(delta_lon)
  tmp <- ifelse(tmp > 1, 1, tmp)
  tmp <- ifelse(tmp < 0, 0, tmp)
  gc_angle <- acos(tmp)
  
  # convert bearing from radians to degrees measured clockwise from due north,
  # and convert gc_angle to great circle distance via radius of earth (km)
  bearing <- bearing*360/(2*pi)
  bearing <- (bearing+360)%%360
  earth_rad <- 6371
  gc_dist <- earth_rad*gc_angle
  
  # return list
  ret <-list(bearing = bearing,
             gc_dist = gc_dist)
  return(ret)
}

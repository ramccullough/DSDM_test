prepare_coords_data <- function (coordinates_file) {
  
  attr <- rgdal::readOGR(coordinates_file, 'Survey_points') %>%
    attributes()
  
  coords <- cbind(attr$coords[,1], attr$coords[,2], attr$data)[,c(1:3)] 
  colnames(coords) <- c('Longitude', 'Latitude', 'Site_point')
  coords$Site_point <- as.character(coords$Site_point)
  
  my_pts <- tidyr::separate(coords, 'Site_point', c('Site', 'Point'))
  my_pts$Site <- as.integer(my_pts$Site)
  return(my_pts)
  
}


# Distance between two points
get_dist = function(x, x_lag, y, y_lag){
  sqrt((x - x_lag)^2 + (y - y_lag)^2)
}

# Convert radians to degrees
rad2deg <- function(rad) {
  return((180 * rad) / pi)
}

get_angle = function(x, x_lag1, x_lag2, y, y_lag1, y_lag2){
      
  # create coordinates for vector AB
  ab_coord_x = x_lag1 - x_lag2
  ab_coord_y = y_lag1 - y_lag2
  # create coordinates for vector BC
  bc_coord_x = x - x_lag1
  bc_coord_y = y - y_lag1
  # create dot product and determinant
  dot = (ab_coord_x * bc_coord_x) + (ab_coord_y * bc_coord_y)
  det = (ab_coord_x * bc_coord_y) - (ab_coord_y * bc_coord_x)
  # generate angle
  angle = rad2deg(atan2(det, dot))

  return(angle)
  
}

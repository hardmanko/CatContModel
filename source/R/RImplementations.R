


circMean_OLD = function(angles, weights=1, degrees=TRUE) {
  if (degrees) {
    angles = d2r(angles)
  }
  if (length(weights) <= 1) {
    weights = rep(weights, length(angles))
  }
  weights = weights / sum(weights)
  
  cosx = sum(cos(angles) * weights)
  sinx = sum(sin(angles) * weights)
  
  rval = atan2(sinx, cosx) %% (2 * pi)
  if (degrees) {
    rval = r2d(rval)
  }
  rval
}


#' Circular Absolute Distance Between Values
#' 
#' The circular absolute distance is the smallest angular distance between two values. It is at most 180 degrees.
#' 
#' @param x Vector of values in degrees or radians.
#' @param y Vector of values in degrees or radians.
#' @param degrees If `TRUE`, x and y are treated as though they are in degrees. If `FALSE`, x and y are treated as being in radians.
#' 
#' @export
circAbsDist = function(x, y, degrees=TRUE) {
  offset = 2 * pi
  if (degrees) {
    offset = 360
  }
  
  x = x %% offset
  y = y %% offset
  
  d1 = abs(x - y)
  d2 = abs((x + offset) - y) %% offset
  d3 = abs(x - (y + offset)) %% offset
  
  pmin(d1,d2,d3)
}

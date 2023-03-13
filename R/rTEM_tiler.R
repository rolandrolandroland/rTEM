#' Tiler function
#' @param pattern input pp3 pattern
#' @param x_tiles,y_tiles,z_tiles number of tiles in each dimension (set equal to 1 for no change)
#' @param z_trim,y_trim,z_trim amount to trim off of each dimension.  Trimming will be applied equall to both sides
#'
#' @description This function takes an input pp3 pattern and tiles it in x, y, and z.  Then, it will trim it
#' according to the `trim` inputs
#' @export
tiler = function(pattern, x_tiles = 1, y_tiles = 1, z_tiles = 1,
                 x_trim = 0, y_trim = 0, z_trim = 0) {
  # get the length of each dimension
  x_length = window(pattern)$domain$xrange[2] - window(pattern)$domain$xrange[1]
  y_length = window(pattern)$domain$yrange[2] - window(pattern)$domain$yrange[1]
  z_length = window(pattern)$domain$zrange[2] - window(pattern)$domain$zrange[1]

  # get new widow
  x_dim =  window(pattern)$domain$xrange[1] + x_length * (x_tiles)
  y_dim =  window(pattern)$domain$yrange[1] + y_length * (y_tiles)
  z_dim =  window(pattern)$domain$zrange[1] + z_length * (z_tiles)
  new_window = as.box3(c(window(pattern)$domain$xrange[1], x_dim),
                       c(window(pattern)$domain$yrange[1], y_dim),
                       c(window(pattern)$domain$zrange[1], z_dim))
  # total number of tiles
  total_tiles = x_tiles * y_tiles * z_tiles
  # total number of points
  total_points = total_tiles * npoints(pattern)

  # save original data
  start_data = pattern$data

  # create matrix will full data
  full_data = pattern$data

  # only add tiles if number of tiles is greater than 1
  if (x_tiles > 1) {
    # new tile for each value x_tiles is above 1
    for (i in 2:x_tiles) {
      # new x_coords are old x_coords with a shift of x_length added
      x_coords =  pattern$data$x + x_length * (i -1)
      dat = start_data
      dat$x = x_coords
      # add the new data to the old data
      full_data = rbind(full_data, dat)
    }
  }
  ## now we tile in y, so we must include the previously tiled data
  start_data = full_data
  if (y_tiles > 1) {
    for (i in 2:y_tiles) {
      y_coords = start_data$y + y_length * (i -1)
      dat = start_data
      dat$y = y_coords
      full_data = rbind(full_data, dat)
    }
  }

  # repeat for z
  start_data = full_data
  if (z_tiles > 1) {
    for (i in 2:z_tiles) {
      z_coords = start_data$z + z_length * (i -1)
      dat = start_data
      dat$z = z_coords
      full_data = rbind(full_data, dat)
    }
  }
  new_pattern = pp3(x =full_data$x, y = full_data$y, z = full_data$z,
                    new_window, marks = full_data$marks)

  ####  get rid of any duplicates that come from tiling with points on the border
  unique_data = as.data.frame.hyperframe(new_pattern$data)
  unique_data = unique(unique_data, margin = 1)
  new_pattern = pp3(x = unique_data$x, y = unique_data$y, z = unique_data$z,
                    new_window, marks = unique_data$marks)
  ##
  ##
  # trim off sides of pattern
  x_min = window(new_pattern)$domain$xrange[1] + x_trim/2
  x_max = window(new_pattern)$domain$xrange[2] - x_trim/2

  y_min = window(new_pattern)$domain$yrange[1] + y_trim/2
  y_max = window(new_pattern)$domain$yrange[2] - y_trim/2

  z_min = window(new_pattern)$domain$zrange[1] + z_trim/2
  z_max = window(new_pattern)$domain$zrange[2] - z_trim/2

  trimmed_window = as.box3(c(x_min, x_max),
                           c(y_min, y_max),
                           c(z_min, z_max))

  trimmed_pattern = pp3(x =new_pattern$data$x, y = new_pattern$data$y,
                        z = new_pattern$data$z,
                        trimmed_window, marks = new_pattern$data$marks)
  trimmed_pattern = subset(trimmed_pattern, subset = trimmed_window)
  return(trimmed_pattern)
}

#' Tile and then relabel
#' @param seed input for `set.seed` function (for reproducibility)
#' @param pp3_full input pp3 pattern for relabeling and tiling
#' @param funcs summary functions to calculate on relabeling
#' @param host_formula,dopant_formula formula for host and dopant marks.  Summary functions
#'  will be calculated from dopant type points
#' @param x_tiles,y_tiles,z_tiles number of tiles in each dimension (set equal to 1 for no change)
#' @param z_trim,y_trim,z_trim amount to trim off of each dimension.  Trimming will be applied equall to both sides
#' @param K_cor edge correction(s) to be used for K function
#' @param G_cor edge correction(s) to be used for G function
#' @param F_cor edge correction(s) to be used for F function
#' @param GXGH_cor edge correction(s) to be used for GXGH function
#' @param GXHG_cor edge correction(s) to be used for GXHG function
#'
#' @description This function takes an input pp3 point pattern, relabels it, then
#' tiles the relabeled point pattern.  Then it calculates summary functions on it
tile_relabel = function(seed, pp3_full, funcs = c("K", "G", "F", "GXGH"), ...,
                        host_formula, dopant_formula,
                        x_tiles = 1, y_tiles = 1, z_tiles = 1,
                        x_trim = 0, y_trim = 0, z_trim = 0,
                        maxKr = 10, nKr = 200, maxGr = 5, nGr = 1000,
                        maxGXGHr = 3, maxGXHGr = 8, nGXr = 1000, vside = 0.3,
                        K_cor = "border", G_cor = "rs", F_cor = "rs",
                        GXGH_cor = "rs", GXHG_cor = "rs") {

  ## perform relabelings
  # Total number of host and dopant type points
  host_total = sum(pp3_full$data$marks== host_formula)
  dopant_total = sum(pp3_full$data$marks== dopant_formula)
  pp3_dopant = subset(pp3_full, marks == dopant_formula)
  set.seed(seed)

  # relabel the input pattern, while maintaining the same proportion of host and dopant points
  pp3_relabeled = rlabel(pp3_full,
                         labels = as.factor(c(rep(dopant_formula, dopant_total),
                                              rep(host_formula, host_total))), permute = TRUE)

  # select only the dopant type points
  #pp3_dopant_relabeled = subset(pp3_relabeled, marks == dopant_formula)

  pp3_tiled = tiler(pp3_relabeled, x_tiles = x_tiles, y_tiles = y_tiles, z_tiles = z_tiles,
                    x_trim = x_trim, y_trim = y_trim, z_trim = z_trim)

  pp3_dopant_tiled = subset(pp3_tiled, marks == dopant_formula)

  # calculate summary functions
  if ("K" %in% funcs) {
    if (K_cor == "border") {
      rrl_K = bK3est(pp3_dopant_tiled, rmax = maxKr, nrval =nKr)
    }
    else {
      rrl_K = K3est(pp3_dopant_tiled, rmax = maxKr, nrval = nKr, correction = K_cor)
    }
  }
  if ("G" %in% funcs) {
    rrl_G = G3est(pp3_dopant_tiled, rmax = maxGr, nrval = nGr, correction = G_cor)
  }
  if ("F" %in% funcs) {
    rrl_F = F3est(pp3_dopant_tiled, rmax = maxGr, nrval = nGr, correction = F_cor, vside = vside)
  }
  if ("GXGH" %in% funcs) {
    rrl_GXGH = G3cross(pp3_tiled, i = dopant_formula, j = host_formula,
                       rmax = maxGXGHr, nrval = nGXr, correction = GXGH_cor)
  }
  if ("GXHG" %in% funcs) {
    rrl_GXHG = G3cross(pp3_tiled, i = host_formula, j = dopant_formula,
                       rmax = maxGXHGr, nrval = nGXr, correction = GXHG_cor)

  }
  all_funcs = c("K", "G", "F", "GXGH", "GXHG")
  all_funcs %in% funcs
  relabs = c("rrl_K", "rrl_G", "rrl_F", "rrl_GXGH", "rrl_GXHG")
  out =lapply(relabs[all_funcs %in% funcs], get, envir = environment())
  names(out) =  relabs[all_funcs %in% funcs]
  return(out)
}







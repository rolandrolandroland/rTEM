#' Preprocess UPP data (developement version)
#'
#' @export
UPP_preprocess_dev = function(upp, guest_pattern,
                              min_thick = NA, max_thick = NA,
                              intensity_upp = 1, n_guests = NA) {
  if (is.na(max_thick)) {
    max_thick = max(upp$domain$zrange)
  }
  if (is.na(min_thick)) {
    min_thick = min(upp$domain$zrange)
  }
  # case 1 guest pattern 2d, UPP pattern is 3d,
  if (is.ppp(guest_pattern)) {
    ## how many points are needed?
    if (is.pp3(upp)) {
      upp = rescale_pattern(upp, intensity_upp)
      size = c(diff(guest_pattern$window$xrange), diff(guest_pattern$window$yrange), diff(guest_pattern$window$zrange), (max_thick - min_thick))
      upp = stitch.size(upp, boxSize = size)
      ## translate
      xdif = upp$domain$xrange[1] - guest_pattern$window$xrange[1]
      ydif = upp$domain$yrange[1] - guest_pattern$window$yrange[1]
      upp = pp3(x = upp$data$x - xdif, y = upp$data$y - ydif, z = upp$data$z,
                owin = c(guest_pattern$window$xrange, guest_pattern$window$yrange, upp$domain$zrange))
    }
    # case 2  guest pattern 2d, UPP is 2d
    else if (is.ppp(upp)) {
      upp = rescale_pattern(upp, intensity_upp)
      size = c(diff(guest_pattern$window$xrange), diff(guest_pattern$window$yrange), diff(guest_pattern$window$zrange))
      upp = stitch.size(upp, boxSize = size)
      xdif = upp$domain$xrange[1] - guest_pattern$window$xrange[1]
      ydif = upp$domain$yrange[1] - guest_pattern$window$yrange[1]
      upp = ppp(x = upp$data$x - xdif, y = upp$data$y - ydif,
                owin = c(guest_pattern$window$xrange, guest_pattern$window$yrange))
    }
  }


  # case 3: guest_pattern is 3d and upp is 3d
  if (is.pp3(guest_pattern)) {
    upp = rescale_pattern(upp, intensity_upp)
    size = c(diff(guest_pattern$domain$xrange), diff(guest_pattern$domain$yrange), diff(guest_pattern$domain$zrange))
    upp = stitch.size(upp, boxSize = size)
    ## translate
    xdif = upp$domain$xrange[1] - guest_pattern$domain$xrange[1]
    ydif = upp$domain$yrange[1] - guest_pattern$domain$yrange[1]
    upp = pp3(x = upp$data$x - xdif, y = upp$data$y - ydif, z = upp$data$z,
              as.box3(c(guest_pattern$domain$xrange, guest_pattern$domain$yrange, guest_pattern$domain$zrange)))
  }
  return(upp)
}


#' Preprocess UPP data
#'
#' @export
UPP_preprocess = function(upp, guest_pattern,
                          min_thick = NA, max_thick = NA,
                          intensity_upp = 1, n_guests = NA) {
  if (is.na(max_thick)) {
    max_thick = max(upp$domain$zrange)
  }
  if (is.na(min_thick)) {
    min_thick = min(upp$domain$zrange)
  }
  # case 1 guest pattern 2d, UPP pattern is 3d,
  if (is.ppp(guest_pattern)) {
    ## how many points are needed?
    if (is.pp3(upp)) {
      upp = rescale_pattern(upp, intensity_upp)
      size = c(diff(guest_pattern$window$xrange), diff(guest_pattern$window$yrange), diff(guest_pattern$window$zrange), (max_thick - min_thick))
      upp = stitch.size(upp, boxSize = size)
      ## translate
      xdif = upp$domain$xrange[1] - guest_pattern$window$xrange[1]
      ydif = upp$domain$yrange[1] - guest_pattern$window$yrange[1]
      upp = pp3(x = upp$data$x - xdif, y = upp$data$y - ydif, z = upp$data$z,
                owin = c(guest_pattern$window$xrange, guest_pattern$window$yrange, upp$domain$zrange))
    }
    # case 2  guest pattern 2d, UPP is 2d
    else if (is.ppp(upp)) {
      upp = rescale_pattern(upp, intensity_upp)
      size = c(diff(guest_pattern$window$xrange), diff(guest_pattern$window$yrange), diff(guest_pattern$window$zrange))
      upp = stitch.size(upp, boxSize = size)
      xdif = upp$domain$xrange[1] - guest_pattern$window$xrange[1]
      ydif = upp$domain$yrange[1] - guest_pattern$window$yrange[1]
      upp = ppp(x = upp$data$x - xdif, y = upp$data$y - ydif,
                owin = c(guest_pattern$window$xrange, guest_pattern$window$yrange))
    }
  }


  # case 3: guest_pattern is 3d and upp is 3d
  if (is.pp3(guest_pattern)) {
    upp = rescale_pattern(upp, intensity_upp)
    size = c(diff(guest_pattern$window$xrange), diff(guest_pattern$window$yrange), diff(guest_pattern$window$zrange), (max_thick - min_thick))
    upp = stitch.size(upp, boxSize = size)
    ## translate
    xdif = upp$domain$xrange[1] - guest_pattern$window$xrange[1]
    ydif = upp$domain$yrange[1] - guest_pattern$window$yrange[1]
    upp = pp3(x = upp$data$x - xdif, y = upp$data$y - ydif, z = upp$data$z,
              owin = c(guest_pattern$window$xrange, guest_pattern$window$yrange, upp$domain$zrange))
  }
  return(upp)
}

#' Preprocess UPP data
#'
#' @export
UPP_preprocess_dev = function(upp, guest_pattern,
                          min_thick = NA, max_thick = NA,
                          intensity_upp = 1, n_guests = NA) {
  if (is.na(max_thick)) {
    max_thick = max(upp$domain$zrange)
  }
  if (is.na(min_thick)) {
    min_thick = min(upp$domain$zrange)
  }
  # case 1 guest pattern 2d, UPP pattern is 3d,
  if (is.ppp(guest_pattern)) {
    ## how many points are needed?
    if (is.pp3(upp)) {
      upp = rescale_pattern(upp, intensity_upp)
      size = c(diff(guest_pattern$window$xrange), diff(guest_pattern$window$yrange), diff(guest_pattern$window$zrange), (max_thick - min_thick))
      upp = stitch.size(upp, boxSize = size)
      ## translate
      xdif = upp$domain$xrange[1] - guest_pattern$window$xrange[1]
      ydif = upp$domain$yrange[1] - guest_pattern$window$yrange[1]
      upp = pp3(x = upp$data$x - xdif, y = upp$data$y - ydif, z = upp$data$z,
                owin = c(guest_pattern$window$xrange, guest_pattern$window$yrange, upp$domain$zrange))
    }
    # case 2  guest pattern 2d, UPP is 2d
    else if (is.ppp(upp)) {
      upp = rescale_pattern(upp, intensity_upp)
      size = c(diff(guest_pattern$window$xrange), diff(guest_pattern$window$yrange), diff(guest_pattern$window$zrange))
      upp = stitch.size(upp, boxSize = size)
      xdif = upp$domain$xrange[1] - guest_pattern$window$xrange[1]
      ydif = upp$domain$yrange[1] - guest_pattern$window$yrange[1]
      upp = ppp(x = upp$data$x - xdif, y = upp$data$y - ydif,
                owin = c(guest_pattern$window$xrange, guest_pattern$window$yrange))
    }
  }


  # case 3: guest_pattern is 3d and upp is 3d
  if (is.pp3(guest_pattern)) {
    upp = rescale_pattern(upp, intensity_upp)
    size = c(diff(guest_pattern$domain$xrange), diff(guest_pattern$domain$yrange), diff(guest_pattern$domain$zrange))
    upp = stitch.size(upp, boxSize = size)
    ## translate
    xdif = upp$domain$xrange[1] - guest_pattern$domain$xrange[1]
    ydif = upp$domain$yrange[1] - guest_pattern$domain$yrange[1]
    upp = pp3(x = upp$data$x - xdif, y = upp$data$y - ydif, z = upp$data$z,
              as.box3(c(guest_pattern$domain$xrange, guest_pattern$domain$yrange, guest_pattern$domain$zrange)))
  }
  return(upp)
}


#' Get max or min
#' @export
get_maxes = function(envelope, observed, inds = NA, base_ind = 1,
                     take = "max",
                     funcs = c("G", "K", "F",  "G2"),
                     G_cor = "km", K_cor = "trans", F_cor = "km",
                     G2_cor = "km", G3_cor = "km", G4_cor = "km", sqrt = "K") {


  out = sapply(funcs, function(func) {
    cor = paste(func, "_cor", sep = "")
    rrl = paste("rrl_", func, sep = "")
    cor = get(cor)

    if ((sqrt == "all") || (sqrt == "K" && func == "K")) {
      take_root = function(vector) {
        return(sqrt(vector))
      }
    }

    # if not, just return the original vector
    else {
      take_root = function(vector) {
        return(vector)
      }
    }
    if (any(!is.na(inds))) {
      #print("ind")
      diffs = sapply(inds, function(ind) {
        df1 = take_root(observed[[func]][[cor]]) - take_root(envelope[[base_ind]][[rrl]]$mmean)
        df2 = take_root(envelope[[ind]][[rrl]]$hi) - take_root(envelope[[base_ind]][[rrl]]$mmean)
        df3 = take_root(observed[[func]][[cor]]) - take_root(envelope[[base_ind]][[rrl]]$mmean)
        df4 = take_root(envelope[[ind]][[rrl]]$lo) - take_root(envelope[[base_ind]][[rrl]]$mmean)
        c(df1, df2, df3, df4)
      })
    }

    else {
      #print("noind")
      diffs = take_root(observed[[func]][[cor]]) - take_root(envelope[[base_ind]]$mmean)


    }
    if (take == "max") {
      max(diffs)
    }
    else if (take == "min") {
      min(diffs)
    }
    # min = min(diffs)
    # max = max(diffs)
    #c("max" = max, "min" =min)
  })
  out
}

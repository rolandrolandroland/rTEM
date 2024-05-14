#' Get p value
#' @param observed a vector
#' @param expected a matrix or dataframe
#' @export
calc_p_value  = function(observed, expected, r = NULL,
                         funcs = c("K", "G", "F"),
                         obs_name = c(NA, "cor_name"), expect_name = c("rrl_", "cor_name"),
                         K_cor = "trans", G_cor = "km",
                         F_cor = "km", GXGH_cor = "km", GXHG_cor = "km",
                         G2_cor = "km") {
  # If we prefer to compute a P value instead of specifying a fixed significance level a, the test procedure is as follows.
  # .... Baddeley paper


  for (i in 2:10) {
    nm =paste("G", i, sep = "")
    if (any(funcs == nm)) {
      assign(paste("G", i, "_cor", sep=""), G2_cor)
    }
  }
  out = lapply(funcs, function(func) {
    cor = paste(func, "_cor", sep = "")
    cor = get(cor)

    obs_name[obs_name == "cor_name"] = cor
    obs_name[obs_name == "rrl_"] = paste("rrl_", func, sep = "")
    obs_name[is.na(obs_name)] = func

    expect_name[expect_name == "cor_name"] = cor
    expect_name[is.na(expect_name)] = func
    expect_name[expect_name == "rrl_"] = paste("rrl_", func, sep = "")
    ## make accepting of multiple data types
    # 1 list
    if (is.list(observed)) {
      if (is.numeric(observed[[1]])) {
        H_obs = as.data.frame(observed[[obs_name[1]]])
        H_obs = H_obs[,1]
        r = observed$r
      }
      #there's another level to it
      else {
        H_obs = as.data.frame(observed[[obs_name[1]]])
        r = H_obs$r
        H_obs = H_obs[,obs_name[2]]
      }
    }

    else if (is.numeric(observed) && !is.null(r)) {
      H_obs = observed
    }

    else if (is.numeric(observed) && is.null(r)) {
      stop("ERROR: your `observed` object does not contain an r vector.  Please set
         the input variable `r` equal to a vector of r values of same length as `observed` object")
    }

    else if (is.data.frame(observed) && obs_name) {
      H_obs = observed[,obs_name[1]]
      r = observed[,"r"]
    }
    else if (is.data.frame(observed) && !obs_name) {
      stop("ERROR: please set input variable `obs_name` equal to the column name for your observed variable")
    }

    ## make accepting of multiple data types for `expected` object
    # 1 list
    if (is.list(expected)) {
      if (is.numeric(expected[[1]])) {
        H_bar = as.data.frame(expected[[expect_name[1]]])
        H_bar = H_bar[,1]
      }
      #there's another level to it
      else {
        H_bar = sapply(expected, function(val) {
          val[[expect_name[1]]][[expect_name[2]]]
        })
        H_bar = as.data.frame(t(H_bar))
        #H_bar = as.data.frame(expected[[expect_name[1]]])
        #H_bar = H_bar[,expect_name[[2]]]
      }
    }

    else if (is.numeric(expected)) {
      H_bar = expected
    }

    else if (is.data.frame(expected) && expect_name) {
      H_bar = expected[,expect_name[1]]
    }
    else if (is.data.frame(expect_name) && !expect_name) {
      stop("ERROR: please set input variable `expect_name` equal to the column name for your expected variable")
    }

    j = sapply(1:length(H_obs), function(i) {
      #med = median(H_bar[,i])
      ##over = sum(H_bar[,i] < H_obs[i] & H_bar[i] < med)
      #under = sum(H_bar[,i] > H_obs[i] & H_bar[i] > med)
      #over + under
      over = sum(H_bar[,i] > H_obs[i])
    })

    m = nrow(H_bar)
    p1 = (j + 1) / (m + 1)
    p2 = (m + 1 - j) /(m+1)
    p = 2*apply(cbind(p1, p2), 1, min)


    out = data.frame(r = r,
                     p = p)
  })
  names(out) = funcs
  out
}


#' Get Features
#' @param expected a list of summary functions for the observed value
#' @param simulated a list of
#' @param funcs list of summary functions to get features for.  currently `G`, `K`, and `F` are supported
#' @param features a list.  Must be either length 1 or same length as `funcs` parameter
#' @export
get_features = function(expected, simulated,
                        funcs = c("G", "K", "F"),
                        features = list(c("max_diff", "min_diff",
                                          "total_norm_diff", "total_norm_diff_squared",
                                          "total_diff", "net_diff",
                                          "max_diff_envelope", "min_diff_envelope",
                                          "total_diff_envelope", "net_diff_envelope",
                                          "T_final", "T_final_ratio")),
                        simulated_outer_name = "rrl_",
                        simulated_inner_name = c("r", "mmean", "lo", "hi"),
                        expected_inner_name = c(),
                        sqrt = "K",
                        K_cor = "trans", G_cor = "km", G2_cor = "km",
                        G3_cor = "km", G4_cor = "km",
                        F_cor = "km", GXGH_cor = "km", GXHG_cor = "km")

{


  ## if features is only of length 1 and funcs is longer, then calculate all features in features[[1]] for each function
  if (length(features) == 1 ) {
    features = lapply(1:length(funcs), function(i) {
      features[[1]]
    })
  }

  ## now make sure funcs and features are the same length
  else if (length(fucs) != length(features)) {
    stop("ERROR:  `features` must be a list of either length 1 or same length as `funcs`")
  }

  out = lapply(1:length(funcs), function(i) {
    func = funcs[i]
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
    feats = features[[i]]
    cor = paste(func, "_cor", sep = "")
    cor = get(cor)
    dr = expected[[func]]$r[2] - expected[[func]]$r[1]

    rrl = paste(simulated_outer_name, func, sep = "")
    env_width =  (take_root(simulated[[rrl]][, simulated_inner_name[4]]) -  take_root(simulated[[rrl]][, simulated_inner_name[3]])) / 2

    if ("max_diff" %in% feats) {
      max_diff = max(take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][,simulated_inner_name[2]]))

    }

    if ("min_diff" %in% feats) {
      min_diff = min(take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][,simulated_inner_name[2]]))
    }

    if ("total_norm_diff" %in% feats) {
      total_norm_diff =sum(abs(take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][,simulated_inner_name[2]])) /env_width, na.rm = TRUE) *dr
    }

    if ("total_norm_diff_squared" %in% feats) {

      total_norm_diff_squared = sum(abs((take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][,simulated_inner_name[2]]))/env_width)^2, na.rm = TRUE ) * dr
    }


    if ("total_diff" %in% feats) {
      total_diff = sum(abs(take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][,simulated_inner_name[2]]))) * dr
    }
    if ("net_diff" %in% feats) {
      net_diff = sum(take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][,simulated_inner_name[2]])) * dr
    }
    if ("max_diff_envelope" %in% feats) {
      lo_diff = take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][, simulated_inner_name[3]])
      lo_diff[lo_diff > 0] = 0
      hi_diff = take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][, simulated_inner_name[4]])
      hi_diff[hi_diff < 0] = 0
      # to see how much it falls outside the envelope, we must take the smaller of the distances to each envelope end
      diffs = cbind(lo_diff, hi_diff)
      diffs = apply(diffs, 1, function(x) {
        x[which.max(abs(x))] })
      max_diff_envelope = max(diffs)
    }
    #plot(simulated[[rrl]]$r, diffs)
    if ("min_diff_envelope" %in% feats) {
      lo_diff = take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][, simulated_inner_name[3]])
      lo_diff[lo_diff > 0] = 0
      hi_diff = take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][, simulated_inner_name[4]])
      hi_diff[hi_diff < 0] = 0
      # to see how much it falls outside the envelope, we must take the smaller of the distances to each envelope end
      diffs = cbind(lo_diff, hi_diff)
      diffs = apply(diffs, 1, function(x) {
        x[which.max(abs(x))] })
      min_diff_envelope = max(diffs)
    }

    if ("total_diff_envelope" %in% feats) {
      lo_diff = take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][, simulated_inner_name[3]])
      lo_diff[lo_diff > 0] = 0
      hi_diff = take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][, simulated_inner_name[4]])
      hi_diff[hi_diff < 0] = 0
      # to see how much it falls outside the envelope, we must take the smaller of the distances to each envelope end
      diffs = cbind(lo_diff, hi_diff)
      diffs = apply(diffs, 1, function(x) {
        x[which.max(abs(x))] })
      total_diff_envelope = sum(abs(diffs))*dr
    }

    if ("net_diff_envelope" %in% feats) {
      lo_diff = take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][, simulated_inner_name[3]])
      lo_diff[lo_diff > 0] = 0
      hi_diff = take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][, simulated_inner_name[4]])
      hi_diff[hi_diff < 0] = 0
      # to see how much it falls outside the envelope, we must take the smaller of the distances to each envelope end
      diffs = cbind(lo_diff, hi_diff)
      diffs = apply(diffs, 1, function(x) {
        x[which.max(abs(x))] })
      net_diff_envelope = sum((diffs))*dr
    }
    if ("T_final" %in% feats) {
      t = calc_T_val_observed(observed = expected, expected =simulated, r = NULL,
                              func = func, rmin = 0, rmax = NA,
                              sqrt = sqrt,
                              c(NA, "cor_name"), expect_name = c("rrl_", "mmean"))
      T_final = max(t$T)

    }

    if ("T_final_ratio" %in% feats) {
      t_obs = calc_T_val_observed(observed = expected, expected =simulated, r = NULL,
                                  func = func, rmin = 0, rmax = NA,
                                  sqrt = sqrt,
                                  c(NA, "cor_name"), expect_name = c("rrl_", "mmean"))
      t_obs = max(t_obs$T)
      t_hi = calc_T_val_observed(observed = expected, expected =simulated[[rrl]]$hi, r = NULL,
                                 func = func, rmin = 0, rmax = NA,
                                 sqrt = sqrt,
                                 c(NA, "cor_name"), expect_name = c("rrl_", "mmean"))
      t_hi = max(t_hi$T)
      t_lo = calc_T_val_observed(observed = expected, expected =simulated[[rrl]]$lo, r = NULL,
                                 func = func, rmin = 0, rmax = NA,
                                 sqrt = sqrt,
                                 c(NA, "cor_name"), expect_name = c("rrl_", "mmean"))
      t_lo = max(t_lo$T)
      t_sim = max(c(t_hi, t_lo))
      T_final_ratio = t_obs/t_sim

    }
    #put them all together
    vals = sapply(feats, get, envir = environment())
    vals
  })
  names(out) = funcs
  return(out)
}

#' Perform T Test
#' @param k using kth largest and smallest values
#' @param m number of simulations/relabelings
#' @export
T_test = function(observed, expected, func = "G",
                  envelope_value = NA, k = NA,
                  m = NA, two_sided = TRUE,
                  rmax = NA, rmin = 0, sqrt = FALSE,
                  rrl_name = "rrl_",
                  K_cor = "trans", G_cor = "km",
                  F_cor = "km", GXGH_cor = "km", GXHG_cor = "km") {
  if (is.na(k) && is.na(envelope_value)) {
    stop("ERROR: You must provide either a k value or envelope value")
  }
  sides = 1
  if (two_sided) {
    sides = 2
  }
  if (is.na(k)) {
    significance = 1 - envelope_value
  }
  else {
    significance  = sides*k/( m+ 1)
  }
  rrl = paste(rrl_name, func, sep = "")
  cor = paste(func, "_cor", sep = "")
  cor = get(cor)
  T_obs = calc_T_val_observed(observed = observed, expected =expected, r = NULL,
                              func = func, rmin = rmin, rmax = rmax,
                              sqrt = "K",
                              c(NA, "cor_name"), expect_name = c("rrl_", "mmean"))
  hi = expected[[rrl]]$hi
  T_env_hi = calc_T_val_observed(observed = hi, expected =expected,
                                 r = expected[[rrl]]$r,
                                 func = func, rmin = rmin, rmax = rmax,
                                 sqrt = "K",
                                 c(NA, "cor_name"), expect_name = c("rrl_", "mmean"))
  lo = expected[[rrl]]$lo
  T_env_lo = calc_T_val_observed(observed = lo, expected =expected,
                                 r = expected[[rrl]]$r,
                                 func = func, rmin = rmin, rmax = rmax,
                                 sqrt = "K",
                                 c(NA, "cor_name"), expect_name = c("rrl_", "mmean"))
  T_expect = max(cbind(T_env_hi$T, T_env_lo$T))
  T_obs = max(T_obs$T)
  ind = length(obs)
  if (T_obs > T_expect) {
    return(print(paste("T test rejects null hypothesis at significance level", significance)))
  }
  if (T_obs <= T_expect) {
    return(print(paste("T test fails to reject null hypothesis at significance level", significance)))
  }
}

#' T Value for envelope
#' @export
calc_T_val_envs <- function(relabelings, func = "K",
                            rmin = 0, rmax = NA,
                            sqrt = FALSE, rrl_name = "rrl_",
                            K_cor = "trans", G_cor = "km",
                            F_cor = "km", GXGH_cor = "km", GXHG_cor = "km") {


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

  rrl = paste(rrl_name, func, sep = "")
  cor = paste(func, "_cor", sep = "")
  cor = get(cor)


  if (is.na(rmax)) {
    rmax = max(relabelings[[1]][[rrl]]$r)
  }
  mmean = sapply(relabelings, function(x) {
    as.data.frame(x[[rrl]])[,cor]
  })
  # get radius and median values
  r <- relabelings[[1]][[rrl]]$r

  med <- apply(mmean, 1, median)
  # get indices that fall within rmin and rmax
  ind <- which(r <= rmax & r >= rmin)
  # get r interval/step size
  interval <- rmax / (length(r) - 1)
  T <- apply(mmean, 2, function(x) {
    sapply(ind, function(i) {
      # T_1 <- ((x[1:i] - med[1:i]) *interval)^2 # difference between each simulation value and median value
      # sum(T_1^2)  # multiply the square of the value by the
      T_1 = (take_root(H_obs[1:i]) - take_root(H_bar[1:i]))^2
      #T_1 <- (x[1:i] - med[1:i])^2 # difference between each simulation value and median value
      sum(T_1) * interval # multiply the square of the value by the
    })
  })

  out = as.data.frame(cbind(r = r[ind], T))
}


#' T value for observed
#' @export
calc_T_val_observed = function(observed, expected, r = NULL,
                               func = "K",
                               rmin = 0, rmax = NA,
                               sqrt = FALSE, obs_name = c(NA, "cor_name"), expect_name = c("rrl_", "mmean"),
                               K_cor = "trans", G_cor = "km",
                               F_cor = "km", GXGH_cor = "km", GXHG_cor = "km") {




  cor = paste(func, "_cor", sep = "")
  cor = get(cor)

  obs_name[obs_name == "cor_name"] = cor
  obs_name[obs_name == "rrl_"] = paste("rrl_", func, sep = "")
  obs_name[is.na(obs_name)] = func

  expect_name[expect_name == "cor_name"] = cor
  expect_name[is.na(expect_name)] = func
  expect_name[expect_name == "rrl_"] = paste("rrl_", func, sep = "")



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

  ## make accepting of multiple data types
  # 1 list
  if (is.list(observed)) {
    if (is.numeric(observed[[1]])) {
      H_obs = as.data.frame(observed[[obs_name[1]]])
      H_obs = H_obs[,1]
      r = observed$r
    }
    #there's another level to it
    else {
      H_obs = as.data.frame(observed[[obs_name[1]]])
      r = H_obs$r
      H_obs = H_obs[,obs_name[2]]
    }
  }

  else if (is.numeric(observed) && !is.null(r)) {
    H_obs = observed
  }

  else if (is.numeric(observed) && is.null(r)) {
    stop("ERROR: your `observed` object does not contain an r vector.  Please set
         the input variable `r` equal to a vector of r values of same length as `observed` object")
  }

  else if (is.data.frame(observed) && obs_name) {
    H_obs = observed[,obs_name[1]]
    r = observed[,"r"]
  }
  else if (is.data.frame(observed) && !obs_name) {
    stop("ERROR: please set input variable `obs_name` equal to the column name for your observed variable")
  }

  ## make accepting of multiple data types for `expected` object
  # 1 list
  if (is.list(expected)) {
    if (is.numeric(expected[[1]])) {
      H_bar = as.data.frame(expected[[expect_name[1]]])
      H_bar = H_bar[,1]
    }
    #there's another level to it
    else {
      H_bar = as.data.frame(expected[[expect_name[1]]])
      H_bar = H_bar[,expect_name[[2]]]
    }
  }

  else if (is.numeric(expected)) {
    H_bar = expected
  }

  else if (is.data.frame(expected) && expect_name) {
    H_bar = expected[,expect_name[1]]
  }
  else if (is.data.frame(expect_name) && !expect_name) {
    stop("ERROR: please set input variable `expect_name` equal to the column name for your expected variable")
  }


  if (is.na(rmax)) {
    rmax = max(r)
  }
  ind <- which(r <= rmax & r >= rmin)
  interval <- r[2] - r[1]
  T = sapply(ind, function(i) {
    T_1 = (take_root(H_obs[1:i]) - take_root(H_bar[1:i]))^2
    sum(T_1) * interval
  })

  out = data.frame(r = r[ind],
                   T = T)
}




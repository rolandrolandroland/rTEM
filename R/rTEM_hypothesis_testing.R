#' Make curve set for observed and simulated data
#' @param obs_mat n x m matrix of m curves each of length n
#' @param sim_mat n x j matrix of j curves each of length n
#' @param r_vals vector of length n values corresponding to radius (or other value) for each obs_mat and sim_mat curves
#' @param ind index (column number) of which obs_mat curve will be set as observed
#' @details This function takes an observed and a simulated matrix and uses the curve_set function to turn it into a curve object to then be used in the global_envelope_test function
#' @export
make_curve_set <- function(obs_mat, sim_mat, r_vals, obs_ind, sim_inds = 1:ncol(sim_inds)){
  stopifnot(nrow(obs_mat)==length(r_vals),
            nrow(sim_mat)==length(r_vals))
  # ***pick ONE observed curve*** or an aggregate
  obs_vec <- obs_mat[,obs_ind]           # e.g. 1st replicate
  curve_set(obs = obs_vec,
            sim = sim_mat[,sim_inds],
            r   = r_vals)
}

#' Make curve set for observed data only
#' @param obs_mat n x m matrix of m curves each of length n
#' @param r_vals vector of length n values corresponding to radius (or other value) for each obs_mat curve
#' @param ind index (column number) of which obs_mat curve will be set as observed
#' @details This function takes an observed matrix and uses the curve_set function to turn it into a
#' curve object to then be used in the global_envelope_test function
#' @export
make_curve_set_expected = function(sim_mat, r_vals, obs_ind, sim_inds = 1:ncol(sim_mat)){
  stopifnot(nrow(sim_mat)==length(r_vals))
  # ***pick ONE observed curve*** or an aggregate
  obs_vec <- sim_mat[,obs_ind]           # e.g. 1st replicate
  curve_set(obs = obs_vec,
            sim = sim_mat[,sim_inds][,-obs_ind, drop = FALSE],
            r   = r_vals)
}

#'  Global envelope test for list of list of curves
#' @param expected list of matrices expected curve values.  Format: `expected[[func]]`
#' func is the summary function name.
#' @param observed list of matrices of observed curve values.  Same format as expected curves
#' @param funcs vector containing summary function names.  These must also be the names of the inner list elements
#' @param type which type of `global_envelope_test` to perform
#' @param r_vals vector of length n values corresponding to radius (or other value) for each obs_mat and sim_mat curves
#' @export
get_global_env_test = function(expected, observed,
                               funcs, type, r_vals, sim_inds) {

  relab_set = lapply(1:ncol(observed[[1]]), function(relab) {
    curves = lapply(funcs, function(func) {
      make_curve_set(observed[[func]],
                     expected[[func]],
                     r_vals,
                     relab,
                     sim_inds = sim_inds)
    })
    res = global_envelope_test(curves,
                               type = type)
  })
}

#'  Global envelope test for list of list of curves
#' @param expected list of matrices expected curve values.  Format: `expected[[func]]`,
#' func is the summary function name.
#' @param observed list of matrices of observed curve values.  Same format as expected curves
#' @param funcs vector containing summary function names.  These must also be the names of the inner list elements
#' @param type which type of `global_envelope_test` to perform
#' @param r_vals vector of length n values corresponding to radius (or other value) for each obs_mat and sim_mat curves
#' @export
get_global_env_test_expected = function(expected, funcs,
                                        type, r_vals, sim_inds) {
  relab_set = lapply(1:ncol(expected[[1]]), function(relab) {
    curves = lapply(funcs, function(func) {
      make_curve_set_expected(expected[[func]],
                              r_vals,
                              relab,
                              sim_inds = sim_inds)
    })
    res = global_envelope_test(curves,
                               type = type)
  })
}

#' Calculate squared area between two curves
#' @param x observed curve
#' @param mu expected curve
#' @param r r values for both curves
#' @export
integrated_sq_dev = function(x, mu, r) {
  sum((x - mu)^2) * diff(r[1:2])
}

#' Gett power value for two lists of lists of curves
#' @param expected list of curves
#' @param observed list of curves
#' @param thresh_val threshold, or `1-significance level` of test
#' @param funcs vector of summary functions
#' @param r_values radius values corresponding to each observed and expected matrix
#' @param n_relabs number of relabelings
#' @export
get_power_val = function(expected, observed, thresh_val = .95, funcs = c("G", "K"), r_vals, n_relabs, dat_form = "matrix") {
  sapply(funcs, function(func) {
    if (dat_form == "matrix" | !is.list(observed)) {
      train = observed[[func]]
      null = expected[[func]]
    }
    else {
      train = (sapply(observed, function(i) { i[[func]]}))
      null =  (sapply(expected, function(i) { i[[func]]}))
    }
    #if (func == "K") {
    #  r_vals = seq(0, maxKr, length.out = nKr)
    #}
    #else {
    #  r_vals = seq(0, maxGr, length.out = nGr)
    #}
    # compute row means.  Potentially use median??
    mu_null = rowMeans(null)
    med_null = apply(null, 1, median)

    # compute test statistic: how much each null and each train set deviate from the null mean (area under curve)
    T_null = apply(null, 2, integrated_sq_dev, mu = mu_null, r = r_vals)
    T_train = apply(train, 2, integrated_sq_dev, mu = mu_null, r = r_vals)

    ## define threshold from null distribution
    threshold = quantile(T_null, thresh_val)

    ##
    sum(T_train > threshold) / (n_relabs)
  })
}


#' Gett power value for two lists of lists of curves
#' @param expected list of curves
#' @param observed list of curves
#' @param thresh_val threshold, or `1-significance level` of test
#' @param funcs vector of summary functions
#' @param r_values radius values corresponding to each observed and expected matrix
#' @param n_relabs number of relabelings
get_power_val_deprecated = function(expected, observed, thresh_val = .95, funcs = c("G", "K"), r_vals, n_relabs) {
  sapply(funcs, function(func) {
    train = observed[[func]]
    null = expected[[func]]

    #if (func == "K") {
    #  r_vals = seq(0, maxKr, length.out = nKr)
    #}
    #else {
    #  r_vals = seq(0, maxGr, length.out = nGr)
    #}
    # compute row means.  Potentially use median??
    mu_null = rowMeans(null)
    med_null = apply(null, 1, median)

    # compute test statistic: how much each null and each train set deviate from the null mean (area under curve)
    T_null = apply(null, 2, integrated_sq_dev, mu = mu_null, r = r_vals)
    T_train = apply(train, 2, integrated_sq_dev, mu = mu_null, r = r_vals)

    ## define threshold from null distribution
    threshold = quantile(T_null, thresh_val)

    ##
    sum(T_train > threshold) / (n_relabs)
  })
}





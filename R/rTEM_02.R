#' nndist_subset
#'
#' @param pattern  point pattern of type \emph{ppp} or \emph{pp3} to calculate distances from
#' @param pattern2 point pattern of type \emph{ppp} or \emph{pp3} to calculate distances to.
#'  If \emph{NULL} then \emph{pattern} will be used
#' @param window object of class \emph{owin} (for ppp) or \emph{box3} (for pp3).
#' Only points inside the window are used to calculate distances from, but points outside the
#' window are still included as potential neighbors. If \emph{NULL} then it includes the
#' entire domain of \emph{pattern}
#' @param drop_isolated if \emph{TRUE} then points that are closer to a boundary than they
#' are to their \emph{kth} NN will be dropped: their distances are set to \emph{NA}
#' @param k an integer or vector of integers that determines which NN's to calculate
#' @param output outputs a list if "list", outputs a matrix otherwise
#' @description Edge corrected nearest neighbor distribution
#' @details This function calculates the distances to the \emph{kth} nearest neighbors (NN's) for
#' a subset of points, defined as those inside of \emph{window}. All points are still considered
#' when fiding the NN's.
#' @export
nndist_subset = function(pattern, pattern2 = NULL,
                         window = NULL, drop_isolated = FALSE,
                         k =1, output = "list") {

  # if window is not defined
  if (is.null(window) & is.pp3(pattern)) {
    window = domain(pattern)
  }
  if (is.null(window) & is.ppp(pattern)) {
    window = pattern$window
  }

  inside = subset(pattern, inside.boxx(coords(pattern) , w =window))
  if (is.null(pattern2)) {
    pattern2 = pattern
  }

  if (identical(pattern, pattern2)) {
    # first on will just be self - self (0 distance)
    dist = nncross(inside, pattern, k = k+1)
  }
  else {
    dist = nncross(inside, pattern2, k = k)
    ## if pattern is a subset of pattern2
    if (all(dist[,1] == rep(0, length(dist[,1])))) {
      dist = nncross(inside, pattern2, k = k+1)
    }
  }
  if (drop_isolated == TRUE) {
    border_dist = bdist.points(inside)
    head(border_dist)
    ## drop points if closer to border than NN (set to NA)
    out = lapply(1:length(k), function(i) {
      #all distances that are closer to neighbor than border
      distance = sapply(1:length(dist[,i]), function(j){
        if (border_dist[j] >= dist[j,i]) {
          dist[j,i]
        }
        else {
          NA
        }
      })
      # index of each neighbor
      which = sapply(1:length(dist[,i]), function(j){
        if (border_dist[j] >= dist[j,i]) {
          dist[j,i+length(k)]
        }
        else {
          NA
        }
      })
      df = data.frame("distance" = distance, "which" = which)
      names(df) = c(paste("dist", i, sep = "."), paste("which", i, sep = "."))
      df
    })
    names(out) = paste(k, "NN")
  }
  else {
    out = lapply(1:length(k), function(i) {
      distance = dist[,i]
      # index of each neighbor
      which = dist[i+length(k)]
      df = data.frame("distance" = distance, "which" = which)
      names(df) = c(paste("dist", i, sep = "."), paste("which", i, sep = "."))
      df
    })
  }
  if (output == "list") {
    return(out)
  }
    else {
    which_mat = out[[1]][,2]
    dist_mat = out[[1]][,1]
    for (i in 2:length(out)) {
      which_mat = cbind(which_mat, out[[i]][,2])
      dist_mat = cbind(dist_mat, out[[i]][,1])
    }
    out2 = cbind(dist_mat, which_mat)
    colnames(out2)[k] = sapply(k, function(x) {
      paste("dist", x, sep = ".")
    })
    colnames(out2)[(length(k)+1):ncol(out2)] = sapply(k, function(x) {
      paste("which", x, sep = ".")
    })
    return (out2)
  }
}

#' Neighbors in shell
#' @param pattern point pattern of class \emph{ppp} or \emph{pp3}
#' @param type mark defining point type
#' @param window object of class \emph{owin} (for ppp) or \emph{box3} (for pp3).
#' Only points inside the window are used to calculate distances from, but points outside the
#' window are still included as potential neighbors. If \emph{NULL} then it includes the
#' entire domain of \emph{pattern}
#' @param k_vec an integer or vector of integers that determines which NN's to calculate
#' @param drop_isolated if \emph{TRUE} then points that are closer to a boundary than they
#' are to their \emph{kth} NN will be dropped: their distances are set to \emph{NA}
#' @description Compare distribution of marks to binomial distribution
#' @details For each k in \emph{k_vec}, calculate how many of the points in
#' \emph{pattern} with mark \emph{type} have  0, 1, .. \emph{k} points of mark \emph{type} in their \emph{k}
#' nearest neighbors.
#' This can then be compared to the expected values of a binomial distribution with n = k and
#' p = fraction of points in \emph{pattern} that have mark \emph{type}
#' @export
neighbors = function(pattern, type = NULL,
                     window = NULL, k_vec = 1,
                     drop_isolated = TRUE,
                     ...) {
  ## k_vec must have every integer for this to work
  k_vec = 1:max(k_vec)

  # subset pattern to those with marks type
  pattern_subset = subset(pattern, marks %in% type)

  # for each k in k_vec, get the ditance from each point to its kth nearest neighbor and
  # and the identity of that NN
    # NA means that the point is closer to border than that NN (only used if drop_isolated is TRUE)
  nndist_g_all = nndist_subset(pattern = pattern_subset, pattern2 = pattern,
                               window = window,
                               drop_isolated = drop_isolated, k = k_vec)

  # now get the identity of each of those nearest neighbors
  g_all_marks = sapply(nndist_g_all, function(i) {
    marks(pattern)[i$which]
  })
  g_all_marks = g_all_marks[complete.cases(g_all_marks),]
  ## for each point, this returns the number of points within that shell that are same type
  num_neighbors = as.data.frame(sapply(k_vec, function(shell_size) {
    apply(g_all_marks, 1, function(i) {
     sum(i[1:shell_size] %in% type)
    })
  }))
  shell_long = unlist(lapply(k_vec, rep, nrow(num_neighbors)))
  num = c()
  for (j in 1:ncol(num_neighbors)) {
    num = c(num, num_neighbors[,j])
  }
  num_neighbors_long = data.frame(num = num,
                                  shell= shell_long)
  sum = sapply(k_vec, function(k) {
    vec = num_neighbors_long %>% filter(shell ==k)
    #print(vec)
    sapply(0:length(k_vec), function(inner) {
      sum(vec$num == (inner))
    })
  })
  sum = as.data.frame(sum)
  rownames(sum)= c(0:length(k_vec))
  colnames(sum) = k_vec
  sum
}

#' nn_binom
#' @param k_vec vector of integers to use as \emph{n} values for binomial distribution
#' @param p probability of success. Should equal the concentration of the point pattern that
#' this will be compared to using the \emph{neighbors()} function
#' @description calculates binomial distributions with \emph{n = k}
#' for each \emph{k} in \emph{k_vec} and probability of success \emph{p}
#' @export
nn_binom = function(k_vec, p) {
  sapply(k_vec, function(shell) {
    dbinom(0:max(k_vec), shell, pcp)
  })
}


#' pp3 downsizer
#' @export
pp3_downsizer = function(ranged_pos, downsize, side = "both") {
  if (side == "both") { ## remove from both sides
    win_new = box3(c(min(ranged_pos$x) + downsize, max(ranged_pos$x) - downsize),
                   c(min(ranged_pos$y) +downsize, max(ranged_pos$y)- downsize),
                   c(min(ranged_pos$z) + downsize, max(ranged_pos$z)- downsize))
    # one with all points
    ranged_pos_new = ranged_pos[ranged_pos$x < max(ranged_pos$x) - downsize &  ranged_pos$x >min(ranged_pos$x) + downsize &
                                  ranged_pos$y < max(ranged_pos$y) - downsize &  ranged_pos$y >min(ranged_pos$y) + downsize &
                                  ranged_pos$z < max(ranged_pos$z) - downsize &  ranged_pos$z >min(ranged_pos$z) + downsize,]
  }
  else if (side == "negative") { # just remove from min
    win_new = box3(c(min(ranged_pos$x) + downsize, max(ranged_pos$x)),
                   c(min(ranged_pos$y) +downsize, max(ranged_pos$y)),
                   c(min(ranged_pos$z) + downsize, max(ranged_pos$z)))
    # one with all points
    ranged_pos_new = ranged_pos[ranged_pos$x < max(ranged_pos$x) &  ranged_pos$x >min(ranged_pos$x) + downsize &
                                  ranged_pos$y < max(ranged_pos$y)  &  ranged_pos$y >min(ranged_pos$y) + downsize &
                                  ranged_pos$z < max(ranged_pos$z)  &  ranged_pos$z >min(ranged_pos$z) + downsize,]
  }
  else if (side == "positive") { # just remove from max
    win_new = box3(c(min(ranged_pos$x) , max(ranged_pos$x) - downsize),
                   c(min(ranged_pos$y), max(ranged_pos$y)- downsize),
                   c(min(ranged_pos$z), max(ranged_pos$z)- downsize))
    # one with all points
    ranged_pos_new = ranged_pos[ranged_pos$x < max(ranged_pos$x) - downsize &  ranged_pos$x >min(ranged_pos$x)  &
                                  ranged_pos$y < max(ranged_pos$y) - downsize &  ranged_pos$y >min(ranged_pos$y)  &
                                  ranged_pos$z < max(ranged_pos$z) - downsize &  ranged_pos$z >min(ranged_pos$z) ,]
  }
  else {
    print("side variable must be either 'both', 'positive', or 'negative'")
    return(NULL)
  }

  pp3_new = createSpat(ranged_pos_new, win = win_new)
  marks(pp3_new) = as.factor(ranged_pos_new$mark)
  return(pp3_new)
}

#' Collpase pp3 into ppp
#' @description function that collapses a pp3 into a ppp
#' @export
pp3_collapse = function(input_data, dim = "z", output = "ppp") {
  all_dims = c("x", "y", "z")
  output_dims = all_dims[!(all_dims %in% dim)]
  if (first(class(input_data)) == "pp3")
  {
    data = input_data$data
    mark_data = marks(input_data)
    # get all ranges
    xrange = input_data$domain$xrange
    yrange = input_data$domain$yrange
    zrange = input_data$domain$zrange

    ranges = c("xrange", "yrange", "zrange")

    # select the ones that we are keeping
    output_range = ranges[!(all_dims %in% dim)]
    range1 = lapply(output_range, get, envir = environment())[[1]]
    range2 = lapply(output_range, get, envir = environment())[[2]]


    win = owin(range1, range2)
  }
  ## otherwise it should be dataframe
  else
  {
    data = input_data
    mark_data = input_data$mark

    xrange = range(data$x)
    yrange = range(data$y)
    zrange = range(data$z)
    ranges = c("xrange", "yrange", "zrange")

    # select the ones that we are keeping
    output_range = ranges[!(all_dims %in% dim)]
    range1 = lapply(output_range, get, envir = environment())[[1]]
    range2 = lapply(output_range, get, envir = environment())[[2]]
    win = owin(range1, range2)

  }

  xvals = data$x
  yvals = data$y
  zvals = data$z
  all_vals = c("xvals", "yvals", "zvals")
  output_vals = all_vals[!(all_dims %in% dim)]
  vals_1 = lapply(output_vals, get, envir = environment())[[1]]
  vals_2 = lapply(output_vals, get, envir = environment())[[2]]

  if (output == "ppp")
  {

    new_data = ppp(vals_1, vals_2, window = win)
    marks(new_data) = mark_data

  }
  else
  {
    new_data = data.frame(vals_1, vals_2)
    colnames(new_data) = output_dims
    new_data$mark = mark_data

  }
  return(new_data)
}

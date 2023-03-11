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

#' Multimer Simulation
#' @param guest_pattern point pattern of class \emph{ppp} or \emph{pp3}.  The final multtimer
#' pattern will match this pattern in class, intensity, and domain.  If this is left as NULL, then
#' the domain will match that of \emph{upp}, will be of the class specified in \emph{output},
#' and have number of guests specified in \emph{n_guest}
#' @param upp point pattern of class \emph{ppp} or \emph{pp3} to use as underlying point pattern.
#' Multimers will be selected from these points.
#' @param output leave as \emph{"guest pattern type"} if specifying \emph{guest_pattern}.  Otherwise, set to \emph{"ppp"} for
#' 2D output or \emph{pp3} for 3D output
#' @param n_guests leave as \emph{NA} if specifying \emph{guest_pattern}. Otherwise, set to
#' integer for number of guests in multimer pattern
#' @param min_thick if \emph{guest_pattern} is 2d (ppp) and \emph{upp} is 3d (pp3) this
#' determines the smallest z value to keep before collapsing into 2d.
#' @param max_thick if \emph{guest_pattern} is 2d (ppp) and \emph{upp} is 3d (pp3) this
#' determines the largest z value to keep before collapsing into 2d.
#' @param ztrim a numeric.  This will trim this amount from both top and bottom (positive and negative z)
#' after multimers are generated and before pp3 pattern is collapsed into ppp.
#' Only applies if \emph{upp} is 3D (\emph{pp3})
#' @param size_fracs = a vector of numerics.  This will be the fraction of all guest points that are
#' assigned to groups of each size.  For example, if \emph{size_fracs = c(0.2, 0.3, 0.1, 0.4)},
#'  then 20\% of guests will be randomly assigned, 30\% will be in dimers, 10\% will be in trimers,
#'  and 30\% will be quadramers (group of 4)
#' @param num_neighbors An integer.  Number of nearest neighbors to select from when
#' forming dimers, trimers, etc..  Must be at least at least one less than the length of \emph{size_fracs}
#' @param sample_method if equal to \emph{"rank"}, the probability of a point of rank \emph{x}
#'  being chosen as a guest is \emph{probs[x]}.  If equal to  \emph{"exp"},
#'  the probability of a point of rank \emph{x} being chosen
#'   as a guest is \emph{probs[x] * exp(-exponent * distances[ranks]))}
#' @param exponent a numeric. If \emph{sample_method = "exp"}, then
#' this is the value of \emph{exponent} in \emph{exp(-exponent * distances[ranks])}
#' @param weights vector of length equal to number of dimensions in \emph{upp}. the weighted distance to each of \emph{num_neighbors} nearest neighbors
#' is calculated using \eqn{\sqrt{w_1 x^2 + w_2 y^2 + w_3 z^2}},
#'  where \emph{weights} = (\eqn{w_1, w_2, w_3}). Set to \emph{c(1, 1, 0)} for vertical dimers.
#' @param probs vector of probabilities.
#' For \eqn{probs = c(p_1, p_2, p_3, p_4)}, the probability of the first NN being
#' selected in \eqn{p_1}, the probability of the second is \eqn{p_2}, and so on
#' @param intensity_upp the \emph{upp} will be rescaled to have this intensity before the
#' marks are assigned. Leave as \emph{NA} to use \emph{upp} as it is
#' @description Under construction. See \code{\link{multimersim}} for stable version.
#'  Simulates multimers (groups of two/dimers, three/trimers, etc.) in
#' underlying point pattern (upp) to match domain and number of points in \emph{guest_pattern}.
#' @details Algorithm steps: Case 1: 2D Guest pattern, 3D UPP
#'  \itemize{
#'  \item{Step 1: Prechecks.  If no guest pattern is used, are both `n_guests` and `output`
#'  explicitly defined?  Is `n_guests` smaller than number of points in the `upp`?  If there are fewer
#'  probabilities than the number of neighbors to be cosidered (`num_neighbors`) then append 0's
#'  onto the  `probs` vector until it is long enough.  Normalize the `size_fracs` vector so that
#'  it sums to 1. Rescale `upp` so that it has an intensity of `intensity_upp`, unless such is left as `NA`}
#'  \item{Step 2:} {Select points in the scaled \emph{upp} that are
#'  inside the x-y limits of \emph{guest_pattern} and the z limits difined by `min_thick` and `max_thick`
#'  (default values are z limits of `upp`) }
#'  \item{Step 3:} {Determine total number of guests to be assigned.
#'  We will scale the number of guests either in the given `guest_pattern` or `n_guests`
#'  by the difference in volumes between the trimmed and untrimmed patterns.  Therefore, say `guest_pattern`
#'  has 1,000 points, `min_thick = 0`, `max_thick = 40`, and `ztrim = 10`, `n_points` will be
#'  increased to 1,000 * (40-0)/(30-10) = 2,000.  This way, after the trimming, about 1,000 guest points
#'  will remain
#'  (for dimers, this is number of guests / 2)
#'  and assign this many points in the scaled subsetted UPP to be guests.
#'   These are the "centroids"}
#'   \item{Step 4:} {Determine the number of groups of each size.  The number of points in groups of
#'   each size is simply the input variable `size_fracs`.  So We can divide this by the number of points
#'   in such a group size (1 for monomers, 2 for dimers, 3 for trimers, etc), and multiply by the
#'   total number of guest points}
#'   \item{Step 5:} {Randomly relabel the scaled, subsetted UPP
#'   so that there is one point labeled guest type for each group and the rest are labeled host type.}
#'   \item{Step 6:} {Randomly sample from the guest type points chosen in Step 5
#'   a number of points equal to the number of groups that need to have two or more points}
#'   \item{Step 7:} {Use \code{\link{create_groups}} to make one NN of each point chosen in Step 6
#'   guest type.  NN's are selected only in points from points assigned as hosts in Step 5.}
#'   \item{Step 8:} {Repeat Steps 6 and 7 for trimers, quadramers, etc.}
#'   \item{Step 9:} {Filter out all guest type points chosen in Steps 5-8 that all within `ztrim` units
#'   of the top (+z) or bottom (-z)}
#'   \item{Step 10:} {Create a 2D point pattern (class `ppp`) by ignoring the z coordinates
#'   of each point.  Label guest type points "G" and host type points "H"}}
#' Case 2: 2D Guest pattern, 2D UPP
#' \itemize{
#' \item{} {Same as steps for Case 1, except ignore Step 2, Step 3, and Step 9.}
#' }
#' Case 3: 3D Guest pattern, 3D UPP
#' \item{} {This will be the same as Case 1, except in Step 10 create a 3D point pattern (class `pp3`)
#' instead of a 2D point pattern}
#' @export
multimersim = function(guest_pattern = NULL, upp, output = "guest pattern type", n_guests = NA,
                           min_thick = NA, max_thick = NA, ztrim = 0,
                           size_fracs = c(0, 1),
                           num_neighbors = 6,
                           sample_method = "rank", exponent = 1,
                           weights = c(1, 1, 1),
                           probs = c(1, 0, 0, 0),
                           intensity_upp = NA) {

  ## check input parameters
  if(is.null(guest_pattern)) {
    if (output == "guest pattern type" || is.na(n_guests)) {
      stop("guest_pattern is not defined: output must be either `ppp` or `pp3` and
         n_guests must be a numeric")
    }
    ## check that guest has fewer points than UPP
    if (n_guests >= spatstat.geom::npoints(upp)) {
      stop("n_guests must be smaller than number of points in upp")
    }
  }
  else if (spatstat.geom::npoints(guest_pattern) >= spatstat.geom::npoints(upp)) {
    stop("guest pattern must have fewer points than upp")
  }

  ## make probability vector same length as num_neighbors
  if (num_neighbors > length(probs)) {
    probs = c(probs, rep(0, num_neighbors - length(probs)))
  }

  if (is.na(intensity_upp)) {
    intensity_upp = sum(spatstat.geom::intensity(upp))
  }

  ## normalize size_fracs so that they add to 1
  size_fracs = size_fracs / sum(size_fracs)


  # rescale UPP pattern to match physical system (1 point per nm)
  upp_scaled = rTEM::rescaler(upp, intensity_upp)

  # case 1 and 2:  guest_pattern is 2d
  if (spatstat.geom::is.ppp(guest_pattern) || output == "ppp") {

    if (is.null(guest_pattern) && spatstat.geom::is.pp3(upp)) {
      box_2d = as.owin(c(upp$domain$xrange,
                         upp$domain$yrange))

    }
    else if (is.null(guest_pattern) && spatstat.geom::is.ppp(upp)) {
      box_2d = upp$window
    }

    else {
      box_2d = guest_pattern$window
      n_guests = npoints(guest_pattern)
    }

    box_area = spatstat.geom::area(box_2d)

    # case 1 UPP pattern is 3d
    if (spatstat.geom::is.pp3(upp)) {
      print("Case 1: Multimers will be generated in 3D UPP and then collapsed into 2D")
      box_3d = domain(upp)

      if (is.na(max_thick)) {
        max_thick = max(upp$domain$zrange)
      }
      if (is.na(min_thick)) {
        min_thick = min(upp$domain$zrange)
      }


      # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
      upp_box = subset(upp_scaled, x >= box_2d$xrange[1] & x <= box_2d$xrange[2] &
                         y >= box_2d$yrange[1] & y <= box_2d$yrange[2] &
                         z >=min_thick & z <= max_thick)

      ## If using ztrim, then must adjust so that trimmed box has correct number of points
      # if not using ztrim, it will be zero and full/final will just be 1
      full_volume = abs(box_2d$xrange[2] - box_2d$xrange[1]) *
        abs(box_2d$yrange[2] - box_2d$yrange[1]) *
        (max_thick - min_thick)
      final_volume = abs(box_2d$xrange[2] - box_2d$xrange[1]) *
        abs(box_2d$yrange[2] - box_2d$yrange[1]) *
        (max_thick - min_thick - (ztrim *2))

      # npoints will be scaled by volume ratios
      ## label n/2 points as guest and the rest as host
      n_guest = n_guests * full_volume/final_volume
      n_total = npoints(upp_box) #* full_volume/final_volume
      n_host = n_total - n_guest


      # determine the number of groups of molecules to be present (if dimers, then n_guest/2)
      group_sizes = 1:length(size_fracs)
      n_groups = round(sum(n_guest * size_fracs / group_sizes), 0) # number of groups of each size
      size_dist = round(n_guest *size_fracs / group_sizes, 0)   # this is the number of groups that are each size

      # A points will be dimers, trimers, etc, B points are isolated
      upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups - size_dist[1]), rep("B",size_dist[1]), rep("C", n_total - n_groups)), permute = TRUE)

      # extract guest points
      #upp_background = subset(upp_labeled, marks == "B")
      # upp_multimers = subset(upp_labeled, marks == "A")
      upp_guest = subset(upp_labeled, marks == "A" | marks == "B")
      # extract host points
      upp_host = subset(upp_labeled, marks == "C")

      ## if we don't want any multimers, we are good now
      if (sum(size_fracs) == size_fracs[1]) {
        chosen = upp_guest$data[upp_guest$data$z > min_thick +ztrim & upp_guest$data$z  < max_thick - ztrim,]
        multimer_box = ppp(x = chosen$x, y = chosen$y, window = box_2d)
      }

      else {
        all_guest = upp_guest
        chosen_points = coords(upp_guest)
        new_host = upp_host
        for (i in 2:length(size_fracs)) {
          # sample the new guest points for the ones that should be trimers
          n_guest_to_add_to = nrow(chosen_points) - size_dist[i-1]
          guests_to_add_to = sample(1:nrow(chosen_points), n_guest_to_add_to)
          guests_to_add_to =pp3(x = chosen_points$x[guests_to_add_to],
                                y = chosen_points$y[guests_to_add_to],
                                z = chosen_points$z[guests_to_add_to], window = box_3d)

          chosen_points = create_groups(num_neighbors = num_neighbors, upp_guest = guests_to_add_to, upp_host = new_host,
                                        probs = probs, weights = weights, sample_method = sample_method, group_size = 2, exponent = exponent)

          all_guest = rbind(coords(all_guest), chosen_points)
          all_guest = pp3(x = all_guest$x, y = all_guest$y, all_guest$z,  window = box_3d)

          ind = match(do.call("paste", chosen_points), do.call("paste", coords(new_host)))
          new_host = pp3(new_host$data$x[-ind], new_host$data$y[-ind], new_host$data$z[-ind], window = box_3d)

        }


        marks(all_guest) = "G"
        marks(new_host) = "H"


        multimer_coords = data.frame(x = c(all_guest$data$x, new_host$data$x),
                                     y = c(all_guest$data$y, new_host$data$y), z = c(all_guest$data$z, new_host$data$z),
                                     marks = c(all_guest$data$marks, new_host$data$marks))#, window = box_3d)
        chosen = multimer_coords[multimer_coords$z > min_thick+ztrim & multimer_coords$z  < max_thick - ztrim,]
        multimer_box = ppp(x = chosen$z, y = chosen$y, marks = as.factor(chosen$marks), window = box_2d)
      }

    }
    ### END CASE 1


    # case 2 UPP is 2d
    else if (spatstat.geom::is.ppp(upp)) {
      print("Case 2: Multimers will be generated in 2D UPP for 2D output")

      # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
      upp_box = subset(upp_scaled, x >= box_2d$xrange[1] & x <= box_2d$xrange[2] &
                         y >= box_2d$yrange[1] & y <= box_2d$yrange[2])

      ## label n/2 points as guest and the rest as host
      n_guest = n_guests
      n_total = npoints(upp_box)
      n_host = n_total - n_guest
      if (n_guest >= n_total) {
        stop("Error: The underlying point pattern (UPP) has fewer points inside the window that contains
             the guest pattern than the guest pattern does.  Retry with a different UPP, or scale down your current one
             so that the intensity (point density) is greater")
      }
      weights = c(weights[1], weights[2])

      # determine the number of groups of molecules to be present (if dimers, then n_guest/2)
      group_sizes = 1:length(size_fracs)
      n_groups = round(sum(n_guest * size_fracs / group_sizes), 0) # number of groups of each size
      size_dist = round(n_guest *size_fracs / group_sizes, 0)   # this is the number of groups that are each size

      # A points will be dimers, trimers, etc, B points are isolated
      upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups - size_dist[1]), rep("B",size_dist[1]), rep("C", n_total - n_groups)), permute = TRUE)

      # extract guest points
      #upp_background = subset(upp_labeled, marks == "B")
      # upp_multimers = subset(upp_labeled, marks == "A")
      upp_guest = subset(upp_labeled, marks == "A" | marks == "B")
      # extract host points
      upp_host = subset(upp_labeled, marks == "C")
      if (sum(size_fracs) == size_fracs[1]) {
        multimer_box = ppp(x = upp_guest$data$x, y = upp_guest$data$y, window = box_2d)

      }

      else {
        chosen_points = coords(upp_guest)
        new_host = upp_host
        all_guest = upp_guest
        for (i in 2:length(size_fracs)) {
          # sample the new guest points for the ones that should be trimers
          n_guest_to_add_to = nrow(chosen_points) - size_dist[i-1]
          guests_to_add_to = sample(1:nrow(chosen_points), n_guest_to_add_to)
          guests_to_add_to =ppp(x = chosen_points$x[guests_to_add_to],
                                y = chosen_points$y[guests_to_add_to],
                                window = box_2d)


          chosen_points = create_groups(num_neighbors = num_neighbors, upp_guest = guests_to_add_to, upp_host = new_host,
                                        probs = probs, weights = weights, sample_method = sample_method, group_size = 2, exponent = exponent)
          all_guest = rbind(coords(all_guest), chosen_points)
          all_guest = ppp(x = all_guest$x, y = all_guest$y,  window = box_2d)

          ind = match(do.call("paste", chosen_points), do.call("paste", coords(new_host)))
          new_host = ppp(new_host$x[-ind], new_host$y[-ind], window = box_2d)

        }

        marks(all_guest) = "G"
        marks(new_host) = "H"
        multimer_box = ppp(x = c(all_guest$x, new_host$x), y = c(all_guest$y, new_host$y),
                           marks = c(all_guest$marks, new_host$marks), window = box_2d)

      }
    }


    return(multimer_box)
    ## END CASE 2
  }

  # case 3: guest_pattern is 3d and upp is 3d
  else if (spatstat.geom::is.pp3(guest_pattern) || output == "pp3") {
    print("Case 3:  Multimers will be generated in 3D for 3D output")
    if (is.null(guest_pattern)) {
      box_3d = domain(upp)
    }

    else {
      box_3d = domain(guest_pattern)
      n_guests = npoints(guest_pattern)
    }

    box_area = spatstat.geom::volume(box_3d)

    if (is.na(max_thick)) {
      max_thick = box_3d$zrange[2]
    }
    if (is.na(min_thick)) {
      min_thick = box_3d$zrange[1]
    }


    # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
    upp_box = subset(upp_scaled, x >= box_3d$xrange[1] & x <= box_3d$xrange[2] &
                       y >= box_3d$yrange[1] & y <= box_3d$yrange[2] &
                       z >=min_thick& z <= max_thick)


    ## If using ztrim, then must adjust so that trimmed box has correct number of points
    # if not using ztrim, it will be zero and full/final will just be 1
    full_volume = abs(box_3d$xrange[2] - box_3d$xrange[1]) *
      abs(box_3d$yrange[2] - box_3d$yrange[1]) *
      (max_thick - min_thick)
    final_volume = abs(box_3d$xrange[2] - box_3d$xrange[1]) *
      abs(box_3d$yrange[2] - box_3d$yrange[1]) *
      (max_thick - min_thick- (ztrim *2))


    # npoints will be scaled by volume ratios
    ## label n/2 points as guest and the rest as host
    n_guest = n_guests * full_volume/final_volume
    n_total = npoints(upp_box) #* full_volume/final_volume
    n_host = n_total - n_guest

    # determine the number of groups of molecules to be present (if dimers, then n_guest/2)
    group_sizes = 1:length(size_fracs)
    n_groups = round(sum(n_guest * size_fracs / group_sizes), 0) # number of groups of each size
    size_dist = round(n_guest *size_fracs / group_sizes, 0)   # this is the number of groups that are each size

    # A points will be dimers, trimers, etc, B points are isolated
    upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups - size_dist[1]), rep("B",size_dist[1]), rep("C", n_total - n_groups)), permute = TRUE)

    # extract guest points
    #upp_background = subset(upp_labeled, marks == "B")
    #upp_multimers = subset(upp_labeled, marks == "A")
    upp_guest = subset(upp_labeled, marks == "A" | marks == "B")
    # extract host points
    upp_host = subset(upp_labeled, marks == "C")
    if (sum(size_fracs) == size_fracs[1]) {
      chosen = upp_guest$data[upp_guest$data$z > min_thick + ztrim & upp_guest$data$z  < max_thick - ztrim,]
      multimer_box = ppp(x = chosen$x, y =chosen$y, window = box_2d)

    }

    else {
      chosen_points = coords(upp_guest)
      new_host = upp_host
      all_guest = upp_guest
      for (i in 2:length(size_fracs)) {
        # sample the new guest points for the ones that should be trimers
        n_guest_to_add_to = nrow(chosen_points) - size_dist[i-1]
        guests_to_add_to = sample(1:nrow(chosen_points), n_guest_to_add_to)

        guests_to_add_to =pp3(x = chosen_points$x[guests_to_add_to],
                              y = chosen_points$y[guests_to_add_to],
                              z = chosen_points$z[guests_to_add_to], window = box_3d)

        chosen_points = create_groups(num_neighbors = num_neighbors, upp_guest = guests_to_add_to, upp_host = new_host,
                                      probs = probs, weights = weights, sample_method = sample_method, group_size = 2, exponent = exponent)
        all_guest = rbind(coords(all_guest), chosen_points)
        all_guest = pp3(x = all_guest$x, y = all_guest$y, all_guest$z,  window = box_3d)

        ind = match(do.call("paste", chosen_points), do.call("paste", coords(new_host)))
        new_host = pp3(new_host$data$x[-ind], new_host$data$y[-ind], new_host$data$z[-ind], window = box_3d)

      }

      marks(all_guest) = "G"
      marks(new_host) = "H"


      multimer_coords = data.frame(x = c(all_guest$data$x, new_host$data$x), y = c(all_guest$data$y, new_host$data$y), z = c(all_guest$data$z, new_host$data$z),
                                   marks = c(all_guest$data$marks, new_host$data$marks))
      chosen = multimer_coords[multimer_coords$z > min_thick + ztrim & multimer_coords$z  < max_thick - ztrim,]
      box_3d_trimmed = box_3d
      box_3d_trimmed$zrange = c(min_thick + ztrim, max_thick - ztrim)
      multimer_box = pp3(x = chosen$z, y = chosen$y, chosen$z, marks = as.factor(chosen$marks), window = box_3d_trimmed)

    }

  }
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

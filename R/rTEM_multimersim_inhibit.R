#' @param OVER overlying point pattern
#' @export
findNearestPoints <- function(OVER, UNDER, n = 1) {
  # Ensure OVER and UNDER are pp3 objects and n is a positive integer
  if (!inherits(OVER, "pp3") || !inherits(UNDER, "pp3")) {
    stop("Both OVER and UNDER must be pp3 objects")
  }
  if (!is.numeric(n) || n <= 0 || n != round(n)) {
    stop("n must be a positive integer")
  }

  i = 0
  nn_same_under = TRUE
  OVER_sub = OVER
  UNDER_sub = UNDER
  chosen_under = c()
  while (sum(nn_same_under)) {
    # Find the nearest neighbors from OVER to UNDER
    # Exclude points in UNDER that are also in OVER
    coords_over <- data.frame(x = OVER_sub$data$x, y = OVER_sub$data$y, z = OVER_sub$data$z)
    coords_under <- data.frame(x = UNDER_sub$data$x, y = UNDER_sub$data$y, z = UNDER_sub$data$z,
                               id = seq(1, nrow(UNDER_sub$data)))

    nn <- nncross(OVER_sub, UNDER_sub, k =(n+i))
    nn_df <- as.data.frame(nn)



    # Ensure no point in UNDER is selected twice
    # listed as index of nn_df
    nn_same_under = duplicated(nn_df[,2]) | nn_df[,2] %in% chosen_under$id

    # listed as index of nn_df
    nn_dupe = nn_df[,1] == 0
    redo = (nn_same_under | nn_dupe)
    nn_unique = !redo

    # Filter nn to exclude these points
    OVER_sub = subset(OVER_sub, redo)

    chosen_under = rbind(chosen_under, coords_under[nn_df[,2][nn_unique],])
    i = i + 1
  }

  # Return the unique nearest neighbors without duplicates
  return(chosen_under)
}




#' Multimer Simulation with Overlying point pattern
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
#'  @param trans_plane plane for translating
#'  @param trans_frac fraction of distance in plane `trans_plane` to reduce distance by
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
#' \itemize{
#' \item{} {This will be the same as Case 1, except in Step 10 create a 3D point pattern (class `pp3`)
#' instead of a 2D point pattern}
#' }
#' @export
multimersim_inhibit = function(guest_pattern = NULL,
                               upp_under, upp_over = NULL, perc_over = 1,
                               output = "guest pattern type", n_guests = NA,
                               min_thick = NA, max_thick = NA, ztrim = 0,
                               size_fracs = c(0, 1),
                               num_neighbors = 6,
                               sample_method = "rank", exponent = 1,
                               weights = c(1, 1, 1),
                               probs = c(1, 0, 0, 0),
                               trans_plane = c("x", "y"), trans_frac = 0.5,
                               intensity_upp_under = NA,
                               intensity_upp_over = NA,
                               preprocessed = FALSE) {

  ## make probability vector same length as num_neighbors
  if (num_neighbors > length(probs)) {
    probs = c(probs, rep(0, num_neighbors - length(probs)))
  }



  ## check input parameters
  if(is.null(guest_pattern)) {
    if (output == "guest pattern type" || is.na(n_guests)) {
      stop("guest_pattern is not defined: output must be either `ppp` or `pp3` and
         n_guests must be a numeric")
    }
    ## check that guest has fewer points than UPP
    if (n_guests >= spatstat.geom::npoints(upp_under)) {
      stop("n_guests must be smaller than number of points in upp_under")
    }
  }
  else if (spatstat.geom::npoints(guest_pattern) >= spatstat.geom::npoints(upp_under)) {
    stop("guest pattern must have fewer points than upp_under")
  }
  ## normalize size_fracs so that they add to 1
  size_fracs = size_fracs / sum(size_fracs)

  ## get rid of any 0's at the end of the size_fracs vector
  for (i in 1:length(size_fracs)) {
    if (size_fracs[length(size_fracs)] == 0) {
      size_fracs = size_fracs[1:(length(size_fracs) - 1)]
    }
  }

  if (!preprocessed) {
    if (is.na(intensity_upp_under)) {
      intensity_upp_under = sum(spatstat.geom::intensity(upp_under))
    }



    # rescale UPP pattern to match physical system (1 point per nm)
    upp_under = rTEM::UPP_preprocess(upp_under, guest_pattern = guest_pattern,
                                     min_thick = min_thick, max_thick = max_thick,
                                     intensity_upp = intensity_upp,
                                     n_guests = n_guests)
    if (!is.null(upp_over)) {
      upp_over = rTEM::UPP_preprocess(upp_over, guest_pattern = guest_pattern,
                                      min_thick = min_thick, max_thick = max_thick,
                                      intensity_upp = intensity_upp,
                                      n_guests = n_guests)
    }

  }
  else {
    # upp_scaled = upp
  }


  # case 1 and 2:  guest_pattern is 2d
  if (spatstat.geom::is.ppp(guest_pattern) || output == "ppp") {

    if (is.null(guest_pattern) && spatstat.geom::is.pp3(upp_under)) {
      box_2d = as.owin(c(upp_under$domain$xrange,
                         upp_under$domain$yrange))

    }
    else if (is.null(guest_pattern) && spatstat.geom::is.ppp(upp_under)) {
      box_2d = upp_under$window
    }

    else {
      box_2d = guest_pattern$window
      n_guests = npoints(guest_pattern)
    }

    box_area = spatstat.geom::area(box_2d)

    # case 1 UPP pattern is 3d
    if (spatstat.geom::is.pp3(upp_under)) {
      print("Case 1: Multimers will be generated in 3D UPP and then collapsed into 2D")
      box_3d = domain(upp_under)
      if (is.na(max_thick)) {
        max_thick = max(upp_under$domain$zrange)
      }
      if (is.na(min_thick)) {
        min_thick = min(upp$domain$zrange)
      }


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
      n_total = npoints(upp_under) #* full_volume/final_volume
      n_host = n_total - n_guest


      # determine the number of groups of molecules to be present (if dimers, then n_guest/2)
      group_sizes = 1:length(size_fracs)
      n_groups = round(sum(n_guest * size_fracs / group_sizes), 0) # number of groups of each size
      size_dist = round(n_guest *size_fracs / group_sizes, 0)   # this is the number of groups that are each size


      # some percentage of start points should be determined by upp_over
      size_dist_over = round(perc_over * size_dist, 0)
      # now subtract that from the n_guest
      size_dist_under = size_dist - size_dist_over

      n_guest_over = sum(size_dist_over * group_sizes)
      n_guest_under = n_guest - n_guest_over

      n_groups_over = sum(size_dist_over)
      n_groups_under = n_groups - n_groups_over

      n_total_over = npoints(upp_over)

      if (n_groups_over > n_total_over) {
        print("overlying pattern does not have enough points")
        return("overlying pattern does not have enough points")
      }

      upp_over_labeled = rlabel(upp_over, labels = c(rep("A",n_groups_over - size_dist_over[1]),
                                                     rep("B",size_dist_over[1]),
                                                     rep("C", n_total_over - n_groups_over)),
                                permute = TRUE)
      # A points will be dimers, trimers, etc, B points are isolated, C is host
      upp_under_labeled = rlabel(upp_under, labels = c(rep("A",n_groups_under - size_dist_under[1]),
                                                       rep("B",size_dist_under[1]),
                                                       rep("C", n_total - n_groups_under)), permute = TRUE)
      # extract guest points
      #upp_background = subset(upp_labeled, marks == "B")
      # upp_multimers = subset(upp_labeled, marks == "A")
      upp_guest_under = subset(upp_under_labeled, marks == "A" | marks == "B")
      # extract host points
      upp_host_under = subset(upp_under_labeled, marks == "C")

      upp_guest_over = subset(upp_over_labeled, marks == "A" | marks == "B")
      # extract host points
      upp_host_over = subset(upp_over_labeled, marks == "C")

      ## if we don't want any multimers, we are good now
      if (sum(size_fracs) == size_fracs[1]) {

        # pick the points in under that are closest to the chosen points in over
        chosen_from_over = findNearestPoints(upp_guest_over, upp_host_under, n = 1)
        upp_guest = data.frame(rbind(upp_guest_under$data[,c("x", "y", "z")],
                                     chosen_from_over[,c("x", "y", "z")]))
        upp_host = as.data.frame(setdiff(coords(upp_host_under), upp_guest[c("x", "y", "z")]))
        upp_guest$marks = "G"
        upp_host$marks = "H"
        multimer_coords = data.frame(x = c(upp_guest$x, upp_host$x),
                                     y = c(upp_guest$y, upp_host$y),
                                     z = c(upp_guest$z, upp_host$z),
                                     marks = c(upp_guest$marks, upp_host$marks))#, window = box_3d)
        chosen = multimer_coords[multimer_coords$z > min_thick+ztrim & multimer_coords$z  < max_thick - ztrim,]
        multimer_box = ppp(x = chosen$x, y = chosen$y, marks = as.factor(chosen$marks), window = box_2d)
      }

      else {
        # chosen_from over are the points in upp_host_under that are closest to
        # each point upp_guest_over
        chosen_from_over = findNearestPoints(upp_guest_over, upp_host_under, n = 1)
        if (nrow(coords(upp_guest_under))) {
          chosen_points = as.data.frame(rbind(coords(upp_guest_under),
                                              chosen_from_over[,c("x", "y", "z")]))
        }
        else {
          chosen_points = as.data.frame(chosen_from_over[,c("x", "y", "z")])
        }


        # Find indices of points in upp_host that are NOT in chosen_points
        # This is done by matching coordinates and then finding which rows in upp_host_coords
        # do NOT have a match in chosen_points_df

        to_keep = as.data.frame(setdiff(coords(upp_host_under), chosen_points[c("x", "y", "z")]))

        new_host = pp3(to_keep$x, to_keep$y, to_keep$z, box_3d)
        all_guest = pp3(x = chosen_points$x, y = chosen_points$y, z = chosen_points$z,
                        marks = chosen_points$marks, box_3d)
        for (i in 2:length(size_fracs)) {
          print(paste("pre ", i, size_fracs[i]))

          # sample the new guest points for the ones that should be trimers
          if (sum(size_fracs[i:length(size_fracs)])) {

            n_guest_to_add_to = nrow(chosen_points) - size_dist[i-1]
            guests_to_add_to = sample(1:nrow(chosen_points), n_guest_to_add_to)
            guests_to_add_to =pp3(x = chosen_points$x[guests_to_add_to],
                                  y = chosen_points$y[guests_to_add_to],
                                  z = chosen_points$z[guests_to_add_to], window = box_3d)
            chosen_points = create_groups(num_neighbors = num_neighbors, upp_guest = guests_to_add_to, upp_host = new_host,
                                          probs = probs, weights = weights,
                                          trans_plane = trans_plane, trans_frac = trans_frac,
                                          sample_method = sample_method, group_size = 2, exponent = exponent)

            all_guest = rbind(coords(all_guest), chosen_points$translated)
            all_guest = pp3(x = all_guest$x, y = all_guest$y, all_guest$z,  window = box_3d)

            ind = match(do.call("paste", chosen_points$original), do.call("paste", coords(new_host)))
            chosen_points = chosen_points$translated
            new_host = pp3(new_host$data$x[-ind], new_host$data$y[-ind], new_host$data$z[-ind], window = box_3d)
          }
        }


        marks(all_guest) = "G"
        marks(new_host) = "H"


        multimer_coords = data.frame(x = c(all_guest$data$x, new_host$data$x),
                                     y = c(all_guest$data$y, new_host$data$y), z = c(all_guest$data$z, new_host$data$z),
                                     marks = c(all_guest$data$marks, new_host$data$marks))#, window = box_3d)
        chosen = multimer_coords[multimer_coords$z > min_thick+ztrim & multimer_coords$z  < max_thick - ztrim,]
        multimer_box = ppp(x = chosen$x, y = chosen$y, marks = as.factor(chosen$marks), window = box_2d)
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
        marks(upp_guest) = "G"
        marks(upp_host) = "H"
        multimer_box = ppp(x = c(upp_guest$x, upp_host$x), y = c(upp_guest$y, upp_host$y),
                           marks = c(upp_guest$marks, upp_host$marks), window = box_2d)

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
                                        probs = probs, weights = weights,
                                        trans_plane = trans_plane, trans_frac = trans_frac,
                                        sample_method = sample_method, group_size = 2, exponent = exponent)

          all_guest = rbind(coords(all_guest), chosen_points$translated)
          all_guest = ppp(x = all_guest$x, y = all_guest$y,  window = box_2d)

          ind = match(do.call("paste", chosen_points$original), do.call("paste", coords(new_host)))
          chosen_points = chosen_points$translated

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
      box_3d = domain(upp_scaled)
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

    box_3d_trimmed = box_3d
    box_3d_trimmed$zrange = c(min_thick + ztrim, max_thick - ztrim)

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
      marks(upp_guest) = "G"
      marks(upp_host) = "H"
      multimer_coords = data.frame(x = c(upp_guest$data$x, upp_host$data$x), y = c(upp_guest$data$y, upp_host$data$y),
                                   z = c(upp_guest$data$z, upp_host$data$z),
                                   marks = c(upp_guest$data$marks, upp_host$data$marks))
      chosen = multimer_coords[multimer_coords$z > min_thick + ztrim & multimer_coords$z  < max_thick - ztrim,]
      multimer_box = pp3(x = chosen$x, y = chosen$y, z =chosen$z, marks = as.factor(chosen$marks), window = box_3d_trimmed)

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
                                      probs = probs, weights = weights,
                                      trans_plane = trans_plane, trans_frac = trans_frac,
                                      sample_method = sample_method, group_size = 2, exponent = exponent)
        all_guest = rbind(coords(all_guest), chosen_points$translated)
        all_guest = pp3(x = all_guest$x, y = all_guest$y, all_guest$z,  window = box_3d)

        ind = match(do.call("paste", chosen_points$original), do.call("paste", coords(new_host)))
        chosen_points = chosen_points$translated
        new_host = pp3(new_host$data$x[-ind], new_host$data$y[-ind], new_host$data$z[-ind], window = box_3d)

      }

      marks(all_guest) = "G"
      marks(new_host) = "H"


      multimer_coords = data.frame(x = c(all_guest$data$x, new_host$data$x), y = c(all_guest$data$y, new_host$data$y), z = c(all_guest$data$z, new_host$data$z),
                                   marks = c(all_guest$data$marks, new_host$data$marks))
      chosen = multimer_coords[multimer_coords$z > min_thick + ztrim & multimer_coords$z  < max_thick - ztrim,]
      multimer_box = pp3(x = chosen$x, y = chosen$y, z =chosen$z, marks = as.factor(chosen$marks), window = box_3d_trimmed)

    }

  }
}

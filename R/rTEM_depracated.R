## these file are depracated:: I am keeping copies of them, but they are not to be used any more
#'
#'
#'
#'

#' Simulate 2D Multimers
#'
#' @description Simulate
multimersim_v01 = function(exp_ppp,  thickness, group_size = 2,
                       num_neighbors = 6, neighbor_number = 1, weights = c(1, 1, 1/4), probs = c(0.6, 0.2, 0.2, 0),
                       rcp_pattern, intensity_rcp = 1, maxGr, maxKr, nGr, nKr) {
  ## make probability vector same length as num_neighbors
  if (num_neighbors > length(probs)) {
    probs = c(probs, rep(0, num_neighbors - length(probs)))
  }

  # get desired intensity
  npoints = exp_ppp$n
  print(npoints)
  box_2d = exp_ppp$window

  box_area = area(box_2d)
  #vol = box_area * thickness
  # intensity_exp = npoints / vol
  # rescale RCP pattern to match physical system (1 point per nm)
  rcp_scaled = rescaler(rcp_pattern, intensity_rcp)

  # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
  rcp_box = subset(rcp_scaled, x > box_2d$xrange[1] & x < box_2d$xrange[2] &
                     y > box_2d$yrange[1] & y < box_2d$yrange[2] &
                     z >0 & z < thickness)
  ## label n/2 points as dopant and the rest as host
  n_irp = npoints
  n_total = length(rcp_box$data$x)
  n_host = n_total - n_irp

  # determine the number of groups of molecules to be present (if multimers, then n_irp/2)
  n_groups = round(n_irp/ group_size,0)
  rcp_labeled = rlabel(rcp_box, labels = c(rep("A",n_groups) , rep("C", length(rcp_box$data$x) - n_groups)), permute = TRUE)
  # extract dopant points
  rcp_dope = subset(rcp_labeled, marks == "A")
  # extract host points
  rcp_host = subset(rcp_labeled, marks == "C")
  # find 6 nearest host to each dopant
  nn_ind= lapply(1:num_neighbors, function(x) {
    nncross(rcp_dope, rcp_host, k = x)
  })

  # get nn distance x, y, and z values for each neighbor
  dist_coords = lapply(1:num_neighbors, function(x) {
    coords(rcp_dope) - coords(rcp_host[nn_ind[[x]][,2]])
  })
  # this just rearranges dist_coords to be grouped by point rather than neighbor order
  holder = c()
  neighbors = lapply(1:nrow(dist_coords[[1]]), function(x) {
    for (i in 1:length(dist_coords)) {
      holder = rbind(holder, dist_coords[[i]][x,])
    }
    return(holder)
  })

  # rank each of the nearest neighbors based upon distance and weights (distance will be mostly x-y distance)
  # rank will be from smallest to largest
  ranks = lapply(1:length(neighbors), function(outer) {
    distances = sapply(1:nrow(neighbors[[outer]]), function(inner) {
      sqrt(weights[1]*neighbors[[outer]]$x[inner]^2 + weights[2]*neighbors[[outer]]$y[inner]^2 +weights[3]*neighbors[[outer]]$z[inner]^2)
    })
    order(distances)
    rank(distances)
  })

  # semi randomly select one of the neighbors, prefering onces with lower x-y distance
  which_neighbor = t(sapply(1:length(ranks), function(i) {
    sample(1:num_neighbors, size =group_size-1, prob = probs[ranks[[i]]])
  }))
  if (group_size == 2) {
    which_neighbor = which_neighbor[1,]
  }

  which_neighbor = cbind(1:(length(which_neighbor)/(group_size-1)), which_neighbor)

  # get the coordinates of the semi randomly selected neighbor
  ind= t(sapply(1:nrow(which_neighbor), function(i) {
    sapply(2:(group_size), function(j) {
      nn_ind[[which_neighbor[i,j]]][i,2]
    })
  }))

  #duplicate = apply(ind, 2, duplicated)
  #duplicate = duplicated(c(ind))
  duplicate = duplicated(ind, MARGIN = 0)
  duplicate_ind = which(duplicate, arr.ind = TRUE)[,1]
  duplicate_ind = unique(duplicate_ind)
  not_duplicate = ind[!duplicate]
  ## points that are not duplicates
  chosen_points = rcp_host$data[not_duplicate]


  ## host pattern with previously used points removed
  host_unique_pp3 = rcp_host[-ind]
  #unique_points = rcp_host$data[-ind[duplicate],]
  print(paste("first duplicate ", sum(duplicate)))


  ### for each duplicate nearest neighbor (whenever the same host is the nearest neighbor for two dopants )
  # take the 2nd nearest neighbor
  while (sum(duplicate >0) ) {

    # find 6 nearest host to each dopant
    next_nn= lapply(1:num_neighbors, function(x) {
      nncross(rcp_dope, host_unique_pp3, k = x)
    })

    ####
    dist_coords = lapply(1:num_neighbors, function(x) {
      coords(rcp_dope) - coords(host_unique_pp3[next_nn[[x]][,2]])
    })
    ####

    # this just rearranges dist_coords to be grouped by point rather than neighbor order
    holder = c()
    neighbors = lapply(1:nrow(dist_coords[[1]]), function(x) {
      for (i in 1:length(dist_coords)) {
        holder = rbind(holder, dist_coords[[i]][x,])
      }
      return(holder)
    })

    # rank each of the nearest neighbors based upon distance and weights (distance will be mostly x-y distance)
    # rank will be from smallest to largest
    ranks = lapply(1:length(neighbors), function(outer) {
      distances = sapply(1:nrow(neighbors[[outer]]), function(inner) {
        sqrt(weights[1]*neighbors[[outer]]$x[inner]^2 + weights[2]*neighbors[[outer]]$y[inner]^2 +weights[3]*neighbors[[outer]]$z[inner]^2)
      })
      order(distances)
      rank(distances)
    })


    # semi randomly select one of the neighbors, prefering onces with lower x-y distance
    which_neighbor = t(sapply(1:length(ranks), function(i) {
      sample(1:num_neighbors, size =group_size-1, prob = probs[ranks[[i]]])
    }))
    if (group_size == 2) {
      which_neighbor = which_neighbor[1,]
    }
    which_neighbor = cbind(1:(length(which_neighbor)/(group_size-1)), which_neighbor)


    # get the coordinates of the semi randomly selected neighbor
    all_ind= t(sapply(1:nrow(which_neighbor), function(i) {
      sapply(2:(group_size), function(j) {
        next_nn[[which_neighbor[i,j]]][i,2]
      })
    }))

    #ind[duplicate] = all_ind[duplicate]
    #all_ind = c(all_ind)
    # duplicate = c(duplicate)
    # extract just the ones needed to replace duplicates
    next_ind =all_ind[duplicate]

    ## use two duplicates - one that has the index of duplicates in entire all_ind, one that has actual new duplicates
    #new = all_ind[duplicate]
    mat_1 = matrix(FALSE, nrow = nrow(all_ind), ncol = ncol(all_ind))
    mat_1[duplicate] = all_ind[duplicate]
    # zeros = test ==0
    duplicate_mat = duplicated(mat_1, MARGIN = 0)
    duplicate_mat[mat_1==0] = FALSE


    ## get duplicates and their index
    duplicate = duplicated(next_ind, MARGIN = 0)
    duplicate_ind = which(duplicate, arr.ind = TRUE)
    duplicate_ind = unique(duplicate_ind)
    not_duplicate = unique(next_ind)
    ## points that are not duplicates
    new_points = host_unique_pp3$data[not_duplicate]
    chosen_points = rbind(chosen_points, new_points)

    ## host pattern with previously used points removed
    host_unique_pp3 = host_unique_pp3[-next_ind]

    ## change duplicate back to full matrix form
    duplicate = duplicate_mat

    print(paste("loop duplicate ", sum(duplicate)))
  }
  dim(chosen_points)

  # get coordinates of all original points and their selected neighbors
  coords_multimer = rbind(coords(rcp_dope), data.frame(x = chosen_points$x, y = chosen_points$y, z =chosen_points$z))
  dim(coords_multimer)

  # make ppp and caluclate summary functions
  rcp_box = ppp(x = coords_multimer$x, y = coords_multimer$y, window = box_2d)
  rcp_G = Gest_nn(rcp_box, correction = "km", k = neighbor_number, r = seq(0, maxGr, length.out = nGr))
  rcp_K =Kest(rcp_box, r = seq(0, maxKr, length.out= nKr))
  summary = list(rrl_G = rcp_G, rrl_K = rcp_K)
  return(summary)
}

#' Relabel multimers
#'
#' multimer
#' DEPRACATED now use relabeler
multimers_relab_v01 = function (exp_ppp, num_relabs, thickness, rcp_list, maxGr, maxKr,
                            nGr, nKr,
                            group_size = 2, num_neighbors = 6, neighbor_number = 1, weights = c(1, 1, 1/4),
                            probs = c(0.6, 0.2, 0.2, 0),
                            ncores = detectCores()) {

  ## make vectors of maximum radii to take T tests to for K and G
  G_rad_vec = seq(0, maxGr, length.out = (nGr/50) +1)
  G_rad_vec = G_rad_vec[2:length(G_rad_vec)]
  K_rad_vec = seq(0, maxKr, length.out = (nKr/50)+1)
  K_rad_vec = K_rad_vec[2:length(K_rad_vec)]

  ## number of relabelings per rcp pattern
  nper = num_relabs/floor(length(rcp_list))
  # index for which rcp pattern to use
  ind = rep(1:length(rcp_list), nper)

  cl = makeForkCluster(ncores, outfile = "outfile_test")
  ##  calculate envelopes for length(ind) rcp multimer patterns
  multimers = parLapply(cl, 1:length(ind), function(x) {
    print(paste("start rep  ", x))
    val = multimersim(exp_ppp, thickness = thickness, group_size = group_size,
                      num_neighbors =num_neighbors, neighbor_number = neighbor_number, weights = weights, probs = probs,
                      rcp_list[[ind[x]]],
                      intensity_rcp = 1, maxGr = maxGr, maxKr = maxKr, nGr = nGr, nKr = nKr)
    print(paste("end rep  ", x))
    val
  })

  clusterExport(cl, "multimers", envir = environment())
  # calculate T values for K and G summaries for a variety of max radii
  multimers_T_G = parLapply(cl, G_rad_vec, function(x) {
    Tcalc(multimers, x, func = "G", rmax =maxGr)
  })
  multimers_T_K = parLapply(cl, K_rad_vec, function(x) {
    Tcalc(multimers, x, func = "K",rmax = maxKr)
  })

  # get the 95% CI envelope
  multimers = rrl_averager(multimers, envelope_value = 0.95,  K_cor = "border", G_cor = "km") # use same name to get rid of the older (and enormous) object
  stopCluster(cl)
  return(list(relabs = multimers, T_G = multimers_T_G, T_K = multimers_T_K))
}



#' Create  RCP Pattern with same intensity and number of points as experimental pattern
#' DEPRACATED now just mulimersim() funciton with group_size = 1
func_rcp_simple_v01 = function(exp_ppp, thickness, rcp_pattern, neighbor_number = 1,
                           intensity_rcp = 1, maxGr, maxKr, nGr, nKr) {
  # get desired intensity
  npoints = exp_ppp$n
  print(npoints)
  box_2d = exp_ppp$window
  box_area = spatstat.geom::area(box_2d)

  #vol = box_area * thickness
  # intensity_exp = npoints / vol
  # rescale RCP pattern to match physical system (1 point per nm)
  rcp_scaled = rescaler(rcp_pattern, intensity_rcp)

  # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
  rcp_box = subset(rcp_scaled, x > box_2d$xrange[1] & x < box_2d$xrange[2] &
                     y > box_2d$yrange[1] & y < box_2d$yrange[2] &
                     z >0 & z < thickness)
  ## label n/2 points as dopant and the rest as host
  n_irp = npoints
  n_total = length(rcp_box$data$x)
  n_host = n_total - n_irp
  rcp_labeled = rlabel(rcp_box, labels = c(rep("A", n_irp) , rep("C", n_host)), permute = TRUE)
  # extract dopant points
  rcp_dope = subset(rcp_labeled, marks == "A")
  # extract host points
  #rcp_host = subset(rcp_labeled, marks == "C")

  coords_dope = coords(rcp_dope)
  #win_box = box3(box_2d$xrange, box_2d$yrange, c(0, thickness))
  rcp_box = ppp(x = coords_dope$x, y = coords_dope$y, window = box_2d)

  rcp_G = Gest_nn(rcp_box, correction = "km", k = neighbor_number, r = seq(0, maxGr, length.out = nGr))
  rcp_K =Kest(rcp_box, r = seq(0, maxKr, length.out= nKr))
  summary = list(rrl_G = rcp_G, rrl_K = rcp_K)
  return(summary)
}

#' Perform num_relabs relabelings
#' DEPRACATED now use relabler
rcp_simple_relab_v01 = function(exp_ppp, num_relabs, thickness, neighbor_number =1, rcp_list, maxGr,
                            maxKr, nGr, nKr, ncores = detectCores()) {

  ## make vectors of maximum radii to take T tests to for K and G
  G_rad_vec = seq(0, maxGr, length.out = (nGr/50) +1)
  G_rad_vec = G_rad_vec[2:length(G_rad_vec)]
  K_rad_vec = seq(0, maxKr, length.out = (nKr/50)+1)
  K_rad_vec = K_rad_vec[2:length(K_rad_vec)]
  nper = num_relabs/floor(length(rcp_list))
  # index for which rcp pattern to use
  ind = rep(1:length(rcp_list), nper)

  ##  calculate envelopes for length(ind) rcp multimer patterns

  cl = makeForkCluster(ncores, outfile = "outfile3")
  rcp_simple = parLapply(cl, 1:length(ind), function(x) {
    func_rcp_simple(exp_ppp, thickness = thickness, rcp_list[[ind[x]]], neighbor_number = neighbor_number,
                    intensity_rcp = 1, maxGr = maxGr, maxKr = maxKr, nGr = nGr, nKr = nKr)
  })


  clusterExport(cl, "rcp_simple", envir = environment())
  # calculate T values for K and G summaries for a variety of max radii
  rcp_simple_T_G = parLapply(cl, G_rad_vec, function(x) {
    Tcalc(rcp_simple, x, func = "G", rmax =maxGr)

  })
  rcp_simple_T_K = parLapply(cl, K_rad_vec, function(x) {
    Tcalc(rcp_simple, x, func = "G", rmax = maxKr)
  })
  # get the 95% CI envelope
  rcp_simple = rrl_averager(rcp_simple, envelope_value = 0.95,  K_cor = "border", G_cor = "km") # use same name to get rid of the older (and enormous) object
  stopCluster(cl)
  return(list(relabs = rcp_simple, T_G = rcp_simple_T_G, T_K = rcp_simple_T_K))
}





#' Multimer Simulation Version  02
#' @param guest_pattern point pattern of class \emph{ppp} or \emph{pp3}.  The final multtimer
#' pattern will match this pattern in class, intensity, and domain.  If this is left as NULL, then
#' the domain will match that of \emph{upp}, will be of the class specified in \emph{output},
#' and have number of guests specified in \emph{n_guest}
#' @param upp point pattern of class \emph{ppp} or \emph{pp3} to use as underlying point pattern.
#' Multimers will be selected from these points.
#' @param output leave as \emph{"guest pattern type"} if specifying \emph{guest_pattern}.  Otherwise, set to \emph{"ppp"} for
#' 2D output or \emph{pp3} for 3D output
#' @param n_guests leave as \emph{NA} if specifying \emph{UPP}. Otherwise, set to
#' integer for number of guests in multimer pattern
#' @param min_thick if \emph{guest_pattern} is 2d (ppp) and \emph{upp} is 3d (pp3) this
#' determines the smallest z value to keep before collapsing into 2d.
#' @param max_thick if \emph{guest_pattern} is 2d (ppp) and \emph{upp} is 3d (pp3) this
#' determines the largest z value to keep before collapsing into 2d.
#' @param ztrim a numeric.  This will trim this amount from both top and bottom (positive and negative z)
#' after multimers are generated and before pp3 pattern is collapsed into ppp.
#' Only applies if \emph{upp} is 3D (\emph{pp3})
#' @param group_size size of clusters.  If equal to 1, then all points will be independent
#'  of each other
#' @param num_neighbors number of nearest neighbors to select from when
#' forming dimers, trimers, etc..
#' @param sample_method if equal to \emph{"rank"}, the probability of a point of rank \emph{x}
#'  being chosen as a guest is \emph{probs[x]}.  If equal to  \emph{"exp"},
#'  the probability of a point of rank \emph{x} being chosen
#'   as a guest is \emph{probs[x] * exp(-distances[ranks]))}
#' @param weights vector of length equal to number of dimensions in \emph{upp}. the weighted distance to each of \emph{num_neighbors} nearest neighbors
#' is calculated using \eqn{\sqrt{w_1 x^2 + w_2 y^2 + w_3 z^2}},
#'  where \emph{weights} = (\eqn{w_1, w_2, w_3}). Set to \emph{c(1, 1, 0)} for vertical dimers.
#' @param probs vector of probabilities.  Should sum to \emph{group_size}-1.
#' For \eqn{probs = c(p_1, p_2, p_3, p_4)}, the probability of the first NN being
#' selected in \eqn{p_1}, the probability of the second is \eqn{p_2}, and so on
#' @param intensity_upp the \emph{upp} will be rescaled to have this intensity before the
#' marks are assigned. Leave as \emph{NA} to use \emph{upp} as it is
#' @description Under construction. See \code{\link{multimersim}} for stable version. Simulates multimers (groups of two/dimers, three/trimers, etc.) in
#' underlying point pattern (upp) to match domain and number of points in \emph{guest_pattern}.
#' @details algorithm steps:
#'  \itemize{
#'  \item{Step 1:} {rescale \emph{upp} to match \emph{intensity_upp}}
#'  \item{Step 2:} {Select points in the scaled \emph{upp} that are
#'  inside the domain of \emph{guest_pattern}. }
#'  \item{Step 3:} {Determine number of guest groups or clusters (for dimers, this is number of guests / 2)
#'  and assign this many points in the scaled subsetted UPP to be guests.
#'   These are the "centroids"}
#'  \item{Step 4:} {Take the \emph{num_neighbors} closest point to each guest}
#'  \item{Step 5:} {Rank each neighbor by weighted distance to centroid}
#'  \item{Step 6:} {Using the probabilities in \emph{probs}, select which neighbors
#'   are to be guests (so that the cluster size is now equal to \emph{group_size})}
#'  \item{Step 7:} {For any duplicates, redo process so that correct number of guests are present}
#'  \item{Step 8:} {If \emph{guest_pattern} is 2D and \emph{UPP} is 3D, remove Z coordinate to
#'  make new pattern 2D}
#'}
#'
#' @export

multimersim_v02 = function(guest_pattern = NULL, upp, output = "guest pattern type", n_guests = NA,
                           min_thick = NA, max_thick = NA, ztrim = 0,
                           group_size = 2, num_neighbors = 6, sample_method = "rank",
                           weights = c(1, 1, 1), probs = c(1, 0, 0, 0),
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

  if (sum(probs != 0) < group_size -1) {
    stop("there must be at least 1 nonzero probability in `probs` for each member of group (length(group_size) -1")
  }
  if (is.na(intensity_upp)) {
    intensity_upp = sum(spatstat.geom::intensity(upp))
  }

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

      if (is.na(max_thick)) {
        max_thick = max(upp$domain$zrange)
      }
      if (is.na(min_thick)) {
        min_thick = min(upp$domain$zrange)
      }
      # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
      upp_box = subset(upp_scaled, x >= box_2d$xrange[1] & x <= box_2d$xrange[2] &
                         y >= box_2d$yrange[1] & y <= box_2d$yrange[2] &
                         z >=min_thick - ztrim & z <= max_thick + ztrim)

      ## If using ztrim, then must adjust so that trimmed box has correct number of points
      # if not using ztrim, it will be zero and full/final will just be 1
      full_volume = abs(box_2d$xrange[2] - box_2d$xrange[1]) *
        abs(box_2d$yrange[2] - box_2d$yrange[1]) *
        (max_thick - min_thick + ztrim *2)
      final_volume = abs(box_2d$xrange[2] - box_2d$xrange[1]) *
        abs(box_2d$yrange[2] - box_2d$yrange[1]) *
        (max_thick - min_thick)

      # npoints will be scaled by volume ratios
      ## label n/2 points as guest and the rest as host
      n_guest = n_guests * full_volume/final_volume
      n_total = npoints(upp_box) #* full_volume/final_volume
      n_host = n_total - n_guest

      # determine the number of groups of molecules to be present (if dimers, then n_guest/2)
      n_groups = round(n_guest/ group_size,0)
      upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups) , rep("C", n_total - n_groups)), permute = TRUE)
      # extract guest points
      upp_guest = subset(upp_labeled, marks == "A")
      # extract host points
      upp_host = subset(upp_labeled, marks == "C")
      if (group_size == 1) {
        multimer_box = ppp(x = upp_guest$data$x, y = upp_guest$data$y, window = box_2d)
        return(multimer_box)
        next
      }
      # find 6 nearest host to each guest
      nn_ind= lapply(1:num_neighbors, function(x) {
        nncross(upp_guest, upp_host, k = x)
      })
      # get nn distance x, y, and z values for each neighbor
      dist_coords = lapply(1:num_neighbors, function(x) {
        coords(upp_guest) - coords(upp_host[nn_ind[[x]][,2]])
      })

      # this just rearranges dist_coords to be grouped by point rather than neighbor order
      holder = c()
      neighbors = lapply(1:nrow(dist_coords[[1]]), function(x) {
        for (i in 1:length(dist_coords)) {
          holder = rbind(holder, dist_coords[[i]][x,])
        }
        return(holder)
      })

      # rank each of the nearest neighbors based upon distance and weights (distance will be mostly x-y distance)
      # rank will be from smallest to largest
      ranks = lapply(neighbors, function(i) {
        get_ranks(i, weights = weights)
      })
      distances = lapply(neighbors, function(i) {
        apply(i, 1, function(inner) {
          sqrt(sum(inner^2))
        })
      })
      # head(distances)

      # semi randomly select group_size -1 neighbords acordingly to the probabilties
      # in probs and and sample_method method `sample_method`
      which_neighbor = t(sapply(1:length(ranks), function (i) {
        pick_neighbor(ranks = ranks[[i]], probs = probs, distances = distances[[i]],
                      sample_method = sample_method, group_size = group_size)
      }))

      if (group_size == 2) {
        which_neighbor = which_neighbor[1,]
      }

      # index
      which_neighbor = cbind(1:(length(which_neighbor)/(group_size-1)), which_neighbor)

      # get the coordinates of the semi randomly selected neighbor
      ind= t(sapply(1:nrow(which_neighbor), function(i) {
        sapply(2:(group_size), function(j) {
          nn_ind[[which_neighbor[i,j]]][i,2]
        })
      }))

      #duplicate = apply(ind, 2, duplicated)
      #duplicate = duplicated(c(ind))
      duplicate = duplicated(ind, MARGIN = 0)
      duplicate_ind = which(duplicate, arr.ind = TRUE)[,1]
      duplicate_ind = unique(duplicate_ind)
      not_duplicate = ind[!duplicate]
      ## points that are not duplicates
      chosen_points = coords(upp_host)[not_duplicate,]


      ## host pattern with previously used points removed
      host_unique = upp_host[-ind]
      #unique_points = upp_host$data[-ind[duplicate],]
      print(paste("first duplicate ", sum(duplicate)))


      ### for each duplicate nearest neighbor (whenever the same host is the nearest neighbor for two guests )
      # take the 2nd nearest neighbor
      while (sum(duplicate >0) ) {

        # find 6 nearest host to each guest
        next_nn= lapply(1:num_neighbors, function(x) {
          nncross(upp_guest, host_unique, k = x)
        })

        ####
        dist_coords = lapply(1:num_neighbors, function(x) {
          coords(upp_guest) - coords(host_unique[next_nn[[x]][,2]])
        })
        ####

        # this just rearranges dist_coords to be grouped by point rather than neighbor order
        holder = c()
        neighbors = lapply(1:nrow(dist_coords[[1]]), function(x) {
          for (i in 1:length(dist_coords)) {
            holder = rbind(holder, dist_coords[[i]][x,])
          }
          return(holder)
        })



        # rank each of the nearest neighbors based upon distance and weights (distance will be mostly x-y distance)
        # rank will be from smallest to largest
        ranks = lapply(neighbors, function(i) {
          get_ranks(i, weights = weights)
        })
        distances = lapply(neighbors, function(i) {
          apply(i, 1, function(inner) {
            sqrt(sum(inner^2))
          })
        })
        # head(distances)

        # semi randomly select group_size -1 neighbords acordingly to the probabilties
        # in probs and and sample_method method `sample_method`
        which_neighbor = t(sapply(1:length(ranks), function (i) {
          pick_neighbor(ranks = ranks[[i]], probs = probs, distances = distances[[i]],
                        sample_method = sample_method, group_size = group_size)
        }))

        if (group_size == 2) {
          which_neighbor = which_neighbor[1,]
        }
        which_neighbor = cbind(1:(length(which_neighbor)/(group_size-1)), which_neighbor)


        # get the coordinates of the semi randomly selected neighbor
        all_ind= t(sapply(1:nrow(which_neighbor), function(i) {
          sapply(2:(group_size), function(j) {
            next_nn[[which_neighbor[i,j]]][i,2]
          })
        }))


        # extract just the ones needed to replace duplicates
        next_ind =all_ind[duplicate]

        ## use two duplicates - one that has the index of duplicates in entire all_ind, one that has actual new duplicates
        mat_1 = matrix(FALSE, nrow = nrow(all_ind), ncol = ncol(all_ind))
        mat_1[duplicate] = all_ind[duplicate]
        duplicate_mat = duplicated(mat_1, MARGIN = 0)
        duplicate_mat[mat_1==0] = FALSE


        ## get duplicates and their index
        duplicate = duplicated(next_ind, MARGIN = 0)
        duplicate_ind = which(duplicate, arr.ind = TRUE)
        duplicate_ind = unique(duplicate_ind)
        not_duplicate = unique(next_ind)
        ## points that are not duplicates
        new_points = coords(host_unique)[not_duplicate,]
        chosen_points = rbind(chosen_points, new_points)

        ## host pattern with previously used points removed
        host_unique = host_unique[-next_ind]

        ## change duplicate back to full matrix form
        duplicate = duplicate_mat

        print(paste("loop duplicate ", sum(duplicate)))
      }
      #dim(chosen_points)

      # get coordinates of all original points and their selected neighbors
      coords_multimer = rbind(coords(upp_guest)[,c(1,2, 3)],
                              data.frame(x = chosen_points$x, y = chosen_points$y, z = chosen_points$z))
      coords_multimer = subset(coords_multimer, z > min_thick &  z < max_thick )
      dim(coords_multimer)

      # make ppp and caluclate summary functions
      multimer_box = ppp(x = coords_multimer$x, y = coords_multimer$y, window = box_2d)
    }
    ### END CASE 1


    # case 2 UPP is 2d
    else if (spatstat.geom::is.ppp(upp)) {
      # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
      upp_box = subset(upp_scaled, x >= box_2d$xrange[1] & x <= box_2d$xrange[2] &
                         y >= box_2d$yrange[1] & y <= box_2d$yrange[2])
      ## label n/2 points as guest and the rest as host
      n_guest = n_guests
      n_total = npoints(upp_box)
      n_host = n_total - n_guest
      weights = c(weights[1], weights[2])

      # determine the number of groups of molecules to be present (if multimers, then n_guest/2)
      n_groups = round(n_guest/ group_size,0)
      upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups) , rep("C", n_total - n_groups)), permute = TRUE)
      # extract guest points
      upp_guest = subset(upp_labeled, marks == "A")
      if (group_size == 1) {
        multimer_box = ppp(x = upp_guest$x, y = upp_guest$y, window = box_2d)
        return(multimer_box)
        next
      }
      # extract host points
      upp_host = subset(upp_labeled, marks == "C")


      # find 6 nearest host to each guest
      nn_ind= lapply(1:num_neighbors, function(x) {
        nncross(upp_guest, upp_host, k = x)
      })
      # get nn distance x, y, and z values for each neighbor
      dist_coords = lapply(1:num_neighbors, function(x) {
        coords(upp_guest) - coords(upp_host[nn_ind[[x]][,2]])
      })
      # this just rearranges dist_coords to be grouped by point rather than neighbor order
      holder = c()
      neighbors = lapply(1:nrow(dist_coords[[1]]), function(x) {
        for (i in 1:length(dist_coords)) {
          holder = rbind(holder, dist_coords[[i]][x,])
        }
        return(holder)
      })




      # rank each of the nearest neighbors based upon distance and weights (distance will be mostly x-y distance)
      # rank will be from smallest to largest
      ranks = lapply(neighbors, function(i) {
        get_ranks(i, weights = weights)
      })
      distances = lapply(neighbors, function(i) {
        apply(i, 1, function(inner) {
          sqrt(sum(inner^2))
        })
      })
      # head(distances)

      # semi randomly select group_size -1 neighbords acordingly to the probabilties
      # in probs and and sample_method method `sample_method`
      which_neighbor = t(sapply(1:length(ranks), function (i) {
        pick_neighbor(ranks = ranks[[i]], probs = probs, distances = distances[[i]],
                      sample_method = sample_method, group_size = group_size)
      }))

      if (group_size == 2) {
        which_neighbor = which_neighbor[1,]
      }

      # index
      which_neighbor = cbind(1:(length(which_neighbor)/(group_size-1)), which_neighbor)

      # get the coordinates of the semi randomly selected neighbor
      ind= t(sapply(1:nrow(which_neighbor), function(i) {
        sapply(2:(group_size), function(j) {
          nn_ind[[which_neighbor[i,j]]][i,2]
        })
      }))

      #duplicate = apply(ind, 2, duplicated)
      #duplicate = duplicated(c(ind))
      duplicate = duplicated(ind, MARGIN = 0)
      duplicate_ind = which(duplicate, arr.ind = TRUE)[,1]
      duplicate_ind = unique(duplicate_ind)
      not_duplicate = ind[!duplicate]
      ## points that are not duplicates
      chosen_points = coords(upp_host)[not_duplicate,]


      ## host pattern with previously used points removed
      host_unique = upp_host[-ind]
      #unique_points = upp_host$data[-ind[duplicate],]
      print(paste("first duplicate ", sum(duplicate)))


      ### for each duplicate nearest neighbor (whenever the same host is the nearest neighbor for two guests )
      # take the 2nd nearest neighbor
      while (sum(duplicate >0) ) {

        # find 6 nearest host to each guest
        next_nn= lapply(1:num_neighbors, function(x) {
          nncross(upp_guest, host_unique, k = x)
        })

        ####
        dist_coords = lapply(1:num_neighbors, function(x) {
          coords(upp_guest) - coords(host_unique[next_nn[[x]][,2]])
        })
        ####

        # this just rearranges dist_coords to be grouped by point rather than neighbor order
        holder = c()
        neighbors = lapply(1:nrow(dist_coords[[1]]), function(x) {
          for (i in 1:length(dist_coords)) {
            holder = rbind(holder, dist_coords[[i]][x,])
          }
          return(holder)
        })

        # semi randomly select group_size -1 neighbords acordingly to the probabilties
        # in probs and and sample_method method `sample_method`
        which_neighbor = t(sapply(1:length(ranks), function (i) {
          pick_neighbor(ranks = ranks[[i]], probs = probs, distances = distances[[i]],
                        sample_method = sample_method, group_size = group_size)
        }))

        if (group_size == 2) {
          which_neighbor = which_neighbor[1,]
        }

        # index
        which_neighbor = cbind(1:(length(which_neighbor)/(group_size-1)), which_neighbor)

        # get the coordinates of the semi randomly selected neighbor
        all_ind= t(sapply(1:nrow(which_neighbor), function(i) {
          sapply(2:(group_size), function(j) {
            next_nn[[which_neighbor[i,j]]][i,2]
          })
        }))


        # extract just the ones needed to replace duplicates
        next_ind =all_ind[duplicate]

        ## use two duplicates - one that has the index of duplicates in entire all_ind, one that has actual new duplicates
        mat_1 = matrix(FALSE, nrow = nrow(all_ind), ncol = ncol(all_ind))
        mat_1[duplicate] = all_ind[duplicate]
        duplicate_mat = duplicated(mat_1, MARGIN = 0)
        duplicate_mat[mat_1==0] = FALSE


        ## get duplicates and their index
        duplicate = duplicated(next_ind, MARGIN = 0)
        duplicate_ind = which(duplicate, arr.ind = TRUE)
        duplicate_ind = unique(duplicate_ind)
        not_duplicate = unique(next_ind)
        ## points that are not duplicates
        new_points = coords(host_unique)[not_duplicate,]
        chosen_points = rbind(chosen_points, new_points)

        ## host pattern with previously used points removed
        host_unique = host_unique[-next_ind]

        ## change duplicate back to full matrix form
        duplicate = duplicate_mat

        print(paste("loop duplicate ", sum(duplicate)))
      }
      dim(chosen_points)

      # get coordinates of all original points and their selected neighbors
      coords_multimer = rbind(coords(upp_guest), data.frame(x = chosen_points$x, y = chosen_points$y))
      dim(coords_multimer)

      # make ppp and caluclate summary functions
      multimer_box = ppp(x = coords_multimer$x, y = coords_multimer$y, window = box_2d)
    }
    return(multimer_box)
    ## END CASE 2
  }

  # case 3: guest_pattern is 3d and upp is 3d
  else if (spatstat.geom::is.pp3(guest_pattern) || output == "pp3") {

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
                       z >=min_thick - ztrim& z <= max_thick + ztrim)


    ## If using ztrim, then must adjust so that trimmed box has correct number of points
    # if not using ztrim, it will be zero and full/final will just be 1
    full_volume = abs(box_3d$xrange[2] - box_3d$xrange[1]) *
      abs(box_3d$yrange[2] - box_3d$yrange[1]) *
      (max_thick - min_thick + ztrim *2)
    final_volume = abs(box_3d$xrange[2] - box_3d$xrange[1]) *
      abs(box_3d$yrange[2] - box_3d$yrange[1]) *
      (max_thick - min_thick)


    # npoints will be scaled by volume ratios
    ## label n/2 points as guest and the rest as host
    n_guest = n_guests * full_volume/final_volume
    n_total = npoints(upp_box) #* full_volume/final_volume
    n_host = n_total - n_guest


    # determine the number of groups of molecules to be present (if multimers, then n_guest/2)
    n_groups = round(n_guest/ group_size,0)
    upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups) , rep("C", n_total - n_groups)), permute = TRUE)
    # extract guest points
    upp_guest = subset(upp_labeled, marks == "A")

    if (group_size == 1) {
      multimer_box = pp3(x = upp_guest$data$x, y = upp_guest$data$y,z = upp_guest$data$z, window = box_3d)
      ind = match(do.call("paste", coords(multimer_box)), do.call("paste", coords(upp_labeled)))
      marks(upp_labeled) = "H"
      marks(upp_labeled)[ind] = "G"
      marks(upp_labeled) = as.factor(marks(upp_labeled))
      return(upp_labeled)
      next
    }
    # extract host points
    upp_host = subset(upp_labeled, marks == "C")
    # find 6 nearest host to each guest
    nn_ind= lapply(1:num_neighbors, function(x) {
      nncross(upp_guest, upp_host, k = x)
    })
    # get nn distance x, y, and z values for each neighbor
    dist_coords = lapply(1:num_neighbors, function(x) {
      coords(upp_guest) - coords(upp_host[nn_ind[[x]][,2]])
    })
    # this just rearranges dist_coords to be grouped by point rather than neighbor order
    holder = c()
    neighbors = lapply(1:nrow(dist_coords[[1]]), function(x) {
      for (i in 1:length(dist_coords)) {
        holder = rbind(holder, dist_coords[[i]][x,])
      }
      return(holder)
    })

    # rank each of the nearest neighbors based upon distance and weights (distance will be mostly x-y distance)
    # rank will be from smallest to largest
    ranks = lapply(neighbors, function(i) {
      get_ranks(i, weights = weights)
    })
    distances = lapply(neighbors, function(i) {
      apply(i, 1, function(inner) {
        sqrt(sum(inner^2))
      })
    })
    # head(distances)

    # semi randomly select group_size -1 neighbords acordingly to the probabilties
    # in probs and and sample_method method `sample_method`
    which_neighbor = t(sapply(1:length(ranks), function (i) {
      pick_neighbor(ranks = ranks[[i]], probs = probs, distances = distances[[i]],
                    sample_method = sample_method, group_size = group_size)
    }))

    if (group_size == 2) {
      which_neighbor = which_neighbor[1,]
    }

    # index
    which_neighbor = cbind(1:(length(which_neighbor)/(group_size-1)), which_neighbor)

    # get the coordinates of the semi randomly selected neighbor
    ind= t(sapply(1:nrow(which_neighbor), function(i) {
      sapply(2:(group_size), function(j) {
        nn_ind[[which_neighbor[i,j]]][i,2]
      })
    }))

    #duplicate = apply(ind, 2, duplicated)
    #duplicate = duplicated(c(ind))
    duplicate = duplicated(ind, MARGIN = 0)
    duplicate_ind = which(duplicate, arr.ind = TRUE)[,1]
    duplicate_ind = unique(duplicate_ind)
    not_duplicate = ind[!duplicate]
    ## points that are not duplicates
    chosen_points = coords(upp_host)[not_duplicate,]


    ## host pattern with previously used points removed
    host_unique = upp_host[-ind]
    #unique_points = upp_host$data[-ind[duplicate],]
    print(paste("first duplicate ", sum(duplicate)))


    ### for each duplicate nearest neighbor (whenever the same host is the nearest neighbor for two guests )
    # take the 2nd nearest neighbor
    while (sum(duplicate >0) ) {

      # find 6 nearest host to each guest
      next_nn= lapply(1:num_neighbors, function(x) {
        nncross(upp_guest, host_unique, k = x)
      })

      ####
      dist_coords = lapply(1:num_neighbors, function(x) {
        coords(upp_guest) - coords(host_unique[next_nn[[x]][,2]])
      })
      ####

      # this just rearranges dist_coords to be grouped by point rather than neighbor order
      holder = c()
      neighbors = lapply(1:nrow(dist_coords[[1]]), function(x) {
        for (i in 1:length(dist_coords)) {
          holder = rbind(holder, dist_coords[[i]][x,])
        }
        return(holder)
      })

      # rank each of the nearest neighbors based upon distance and weights (distance will be mostly x-y distance)
      # rank will be from smallest to largest
      ranks = lapply(neighbors, function(i) {
        get_ranks(i, weights = weights)
      })
      distances = lapply(neighbors, function(i) {
        apply(i, 1, function(inner) {
          sqrt(sum(inner^2))
        })
      })
      # head(distances)

      # semi randomly select group_size -1 neighbords acordingly to the probabilties
      # in probs and and sample_method method `sample_method`
      which_neighbor = t(sapply(1:length(ranks), function (i) {
        pick_neighbor(ranks = ranks[[i]], probs = probs, distances = distances[[i]],
                      sample_method = sample_method, group_size = group_size)
      }))

      if (group_size == 2) {
        which_neighbor = which_neighbor[1,]
      }

      # index
      which_neighbor = cbind(1:(length(which_neighbor)/(group_size-1)), which_neighbor)

      # get the coordinates of the semi randomly selected neighbor
      all_ind= t(sapply(1:nrow(which_neighbor), function(i) {
        sapply(2:(group_size), function(j) {
          next_nn[[which_neighbor[i,j]]][i,2]
        })
      }))


      # extract just the ones needed to replace duplicates
      next_ind =all_ind[duplicate]

      ## use two duplicates - one that has the index of duplicates in entire all_ind, one that has actual new duplicates
      mat_1 = matrix(FALSE, nrow = nrow(all_ind), ncol = ncol(all_ind))
      mat_1[duplicate] = all_ind[duplicate]
      duplicate_mat = duplicated(mat_1, MARGIN = 0)
      duplicate_mat[mat_1==0] = FALSE


      ## get duplicates and their index
      duplicate = duplicated(next_ind, MARGIN = 0)
      duplicate_ind = which(duplicate, arr.ind = TRUE)
      duplicate_ind = unique(duplicate_ind)
      not_duplicate = unique(next_ind)
      ## points that are not duplicates
      new_points = coords(host_unique)[not_duplicate,]
      chosen_points = rbind(chosen_points, new_points)

      ## host pattern with previously used points removed
      host_unique = host_unique[-next_ind]

      ## change duplicate back to full matrix form
      duplicate = duplicate_mat

      print(paste("loop duplicate ", sum(duplicate)))
    }
    dim(chosen_points)

    # get coordinates of all original points and their selected neighbors
    coords_multimer = rbind(coords(upp_guest),
                            data.frame(x = chosen_points$x, y = chosen_points$y, z = chosen_points$z))
    coords_multimer = subset(coords_multimer, z > min_thick &  z < max_thick )
    #dim(coords_multimer)

    # make ppp and caluclate summary functions
    multimer_box = pp3(x = coords_multimer$x, y = coords_multimer$y, z = coords_multimer$z, window = box_3d)
    ind = match(do.call("paste", coords(multimer_box)), do.call("paste", coords(upp_box)))
    marks(upp_box) = "H"
    marks(upp_box)[ind] = "G"
    marks(upp_box) = as.factor(marks(upp_box))
    return(upp_box)
  }
}




#' Multimer Simulation Version 03
#' @param guest_pattern point pattern of class \emph{ppp} or \emph{pp3}.  The final multtimer
#' pattern will match this pattern in class, intensity, and domain.  If this is left as NULL, then
#' the domain will match that of \emph{upp}, will be of the class specified in \emph{output},
#' and have number of guests specified in \emph{n_guest}
#' @param upp point pattern of class \emph{ppp} or \emph{pp3} to use as underlying point pattern.
#' Multimers will be selected from these points.
#' @param output leave as \emph{"guest pattern type"} if specifying \emph{guest_pattern}.  Otherwise, set to \emph{"ppp"} for
#' 2D output or \emph{pp3} for 3D output
#' @param n_guests leave as \emph{NA} if specifying \emph{UPP}. Otherwise, set to
#' integer for number of guests in multimer pattern
#' @param min_thick if \emph{guest_pattern} is 2d (ppp) and \emph{upp} is 3d (pp3) this
#' determines the smallest z value to keep before collapsing into 2d.
#' @param max_thick if \emph{guest_pattern} is 2d (ppp) and \emph{upp} is 3d (pp3) this
#' determines the largest z value to keep before collapsing into 2d.
#' @param ztrim a numeric.  This will trim this amount from both top and bottom (positive and negative z)
#' after multimers are generated and before pp3 pattern is collapsed into ppp.
#' Only applies if \emph{upp} is 3D (\emph{pp3})
#' @param group_size size of clusters.  If equal to 1, then all points will be independent
#'  of each other
#' @param num_neighbors number of nearest neighbors to select from when
#' forming dimers, trimers, etc..
#' @param sample_method if equal to \emph{"rank"}, the probability of a point of rank \emph{x}
#'  being chosen as a guest is \emph{probs[x]}.  If equal to  \emph{"exp"},
#'  the probability of a point of rank \emph{x} being chosen
#'   as a guest is \emph{probs[x] * exp(-distances[ranks]))}
#' @param weights vector of length equal to number of dimensions in \emph{upp}. the weighted distance to each of \emph{num_neighbors} nearest neighbors
#' is calculated using \eqn{\sqrt{w_1 x^2 + w_2 y^2 + w_3 z^2}},
#'  where \emph{weights} = (\eqn{w_1, w_2, w_3}). Set to \emph{c(1, 1, 0)} for vertical dimers.
#' @param probs vector of probabilities.  Should sum to \emph{group_size}-1.
#' For \eqn{probs = c(p_1, p_2, p_3, p_4)}, the probability of the first NN being
#' selected in \eqn{p_1}, the probability of the second is \eqn{p_2}, and so on
#' @param intensity_upp the \emph{upp} will be rescaled to have this intensity before the
#' marks are assigned. Leave as \emph{NA} to use \emph{upp} as it is
#' @description Under construction. See \code{\link{multimersim}} for stable version. Simulates multimers (groups of two/dimers, three/trimers, etc.) in
#' underlying point pattern (upp) to match domain and number of points in \emph{guest_pattern}.
#' @details algorithm steps:
#'  \itemize{
#'  \item{Step 1:} {rescale \emph{upp} to match \emph{intensity_upp}}
#'  \item{Step 2:} {Select points in the scaled \emph{upp} that are
#'  inside the domain of \emph{guest_pattern}. }
#'  \item{Step 3:} {Determine number of guest groups or clusters (for dimers, this is number of guests / 2)
#'  and assign this many points in the scaled subsetted UPP to be guests.
#'   These are the "centroids"}
#'  \item{Step 4:} {Take the \emph{num_neighbors} closest point to each guest}
#'  \item{Step 5:} {Rank each neighbor by weighted distance to centroid}
#'  \item{Step 6:} {Using the probabilities in \emph{probs}, select which neighbors
#'   are to be guests (so that the cluster size is now equal to \emph{group_size})}
#'  \item{Step 7:} {For any duplicates, redo process so that correct number of guests are present}
#'  \item{Step 8:} {If \emph{guest_pattern} is 2D and \emph{UPP} is 3D, remove Z coordinate to
#'  make new pattern 2D}
#'}
#'
#' @export
multimersim_v03 = function(guest_pattern = NULL, upp, output = "guest pattern type", n_guests = NA,
                           min_thick = NA, max_thick = NA, ztrim = 0,
                           group_size = 2, num_neighbors = 6, sample_method = "rank",
                           weights = c(1, 1, 1), probs = c(1, 0, 0, 0),
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

  if (sum(probs != 0) < group_size -1) {
    stop("there must be at least 1 nonzero probability in `probs` for each member of group (length(group_size) -1")
  }
  if (is.na(intensity_upp)) {
    intensity_upp = sum(spatstat.geom::intensity(upp))
  }

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

      if (is.na(max_thick)) {
        max_thick = max(upp$domain$zrange)
      }
      if (is.na(min_thick)) {
        min_thick = min(upp$domain$zrange)
      }
      # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
      upp_box = subset(upp_scaled, x >= box_2d$xrange[1] & x <= box_2d$xrange[2] &
                         y >= box_2d$yrange[1] & y <= box_2d$yrange[2] &
                         z >=min_thick - ztrim & z <= max_thick + ztrim)

      ## If using ztrim, then must adjust so that trimmed box has correct number of points
      # if not using ztrim, it will be zero and full/final will just be 1
      full_volume = abs(box_2d$xrange[2] - box_2d$xrange[1]) *
        abs(box_2d$yrange[2] - box_2d$yrange[1]) *
        (max_thick - min_thick + ztrim *2)
      final_volume = abs(box_2d$xrange[2] - box_2d$xrange[1]) *
        abs(box_2d$yrange[2] - box_2d$yrange[1]) *
        (max_thick - min_thick)

      # npoints will be scaled by volume ratios
      ## label n/2 points as guest and the rest as host
      n_guest = n_guests * full_volume/final_volume
      n_total = npoints(upp_box) #* full_volume/final_volume
      n_host = n_total - n_guest

      # determine the number of groups of molecules to be present (if dimers, then n_guest/2)
      n_groups = round(n_guest/ group_size,0)
      upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups) , rep("C", n_total - n_groups)), permute = TRUE)
      # extract guest points
      upp_guest = subset(upp_labeled, marks == "A")
      # extract host points
      upp_host = subset(upp_labeled, marks == "C")
      if (group_size == 1) {
        chosen = upp_guest[upp_guest$data$z > min_thick & upp_guest$data$z  < max_thick]
        multimer_box = ppp(x = upp_guest$data$x[chosen], y = upp_guest$data$y[chosen], window = box_2d)
        return(multimer_box)
        next
      }
      chosen_points = create_groups(num_neighbors, upp_guest, upp_host, probs, weights, sample_method, group_size)


      # get coordinates of all original points and their selected neighbors
      coords_multimer = rbind(coords(upp_guest)[,c(1,2, 3)],
                              data.frame(x = chosen_points$x, y = chosen_points$y, z = chosen_points$z))
      coords_multimer = subset(coords_multimer, z > min_thick &  z < max_thick )
      dim(coords_multimer)

      # make ppp and caluclate summary functions
      multimer_box = ppp(x = coords_multimer$x, y = coords_multimer$y, window = box_2d)
    }
    ### END CASE 1


    # case 2 UPP is 2d
    else if (spatstat.geom::is.ppp(upp)) {
      # subset so that it is now only includes the points in the xy range of the TEM points and z range of thickness
      upp_box = subset(upp_scaled, x >= box_2d$xrange[1] & x <= box_2d$xrange[2] &
                         y >= box_2d$yrange[1] & y <= box_2d$yrange[2])
      ## label n/2 points as guest and the rest as host
      n_guest = n_guests
      n_total = npoints(upp_box)
      n_host = n_total - n_guest
      weights = c(weights[1], weights[2])

      # determine the number of groups of molecules to be present (if multimers, then n_guest/2)
      n_groups = round(n_guest/ group_size,0)
      upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups) , rep("C", n_total - n_groups)), permute = TRUE)
      # extract guest points
      upp_guest = subset(upp_labeled, marks == "A")
      if (group_size == 1) {
        multimer_box = ppp(x = upp_guest$x, y = upp_guest$y, window = box_2d)
        return(multimer_box)
        next
      }
      # extract host points
      upp_host = subset(upp_labeled, marks == "C")

      chosen_points = create_groups(num_neighbors, upp_guest, upp_host, probs, weights, sample_method, group_size)



      # get coordinates of all original points and their selected neighbors
      coords_multimer = rbind(coords(upp_guest), data.frame(x = chosen_points$x, y = chosen_points$y))
      dim(coords_multimer)

      # make ppp and caluclate summary functions
      multimer_box = ppp(x = coords_multimer$x, y = coords_multimer$y, window = box_2d)
    }
    return(multimer_box)
    ## END CASE 2
  }

  # case 3: guest_pattern is 3d and upp is 3d
  else if (spatstat.geom::is.pp3(guest_pattern) || output == "pp3") {

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
                       z >=min_thick - ztrim& z <= max_thick + ztrim)


    ## If using ztrim, then must adjust so that trimmed box has correct number of points
    # if not using ztrim, it will be zero and full/final will just be 1
    full_volume = abs(box_3d$xrange[2] - box_3d$xrange[1]) *
      abs(box_3d$yrange[2] - box_3d$yrange[1]) *
      (max_thick - min_thick + ztrim *2)
    final_volume = abs(box_3d$xrange[2] - box_3d$xrange[1]) *
      abs(box_3d$yrange[2] - box_3d$yrange[1]) *
      (max_thick - min_thick)


    # npoints will be scaled by volume ratios
    ## label n/2 points as guest and the rest as host
    n_guest = n_guests * full_volume/final_volume
    n_total = npoints(upp_box) #* full_volume/final_volume
    n_host = n_total - n_guest


    # determine the number of groups of molecules to be present (if multimers, then n_guest/2)
    n_groups = round(n_guest/ group_size,0)
    upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups) , rep("C", n_total - n_groups)), permute = TRUE)
    # extract guest points
    upp_guest = subset(upp_labeled, marks == "A")

    if (group_size == 1) {
      multimer_box = pp3(x = upp_guest$data$x, y = upp_guest$data$y,z = upp_guest$data$z, window = box_3d)
      ind = match(do.call("paste", coords(multimer_box)), do.call("paste", coords(upp_labeled)))
      marks(upp_labeled) = "H"
      marks(upp_labeled)[ind] = "G"
      marks(upp_labeled) = as.factor(marks(upp_labeled))
      return(upp_labeled)
      next
    }
    # extract host points
    upp_host = subset(upp_labeled, marks == "C")

    # create groups/multimers of size group_size by using sample_method to select from the nearest num_neighbors
    # of each point in upp_guest
    chosen_points = create_groups(num_neighbors, upp_guest, upp_host, probs, weights, sample_method, group_size)

    # get coordinates of all original points and their selected neighbors
    coords_multimer = rbind(coords(upp_guest),
                            data.frame(x = chosen_points$x, y = chosen_points$y, z = chosen_points$z))
    coords_multimer = subset(coords_multimer, z > min_thick &  z < max_thick )
    #dim(coords_multimer)

    # make ppp and caluclate summary functions
    multimer_box = pp3(x = coords_multimer$x, y = coords_multimer$y, z = coords_multimer$z, window = box_3d)
    ind = match(do.call("paste", coords(multimer_box)), do.call("paste", coords(upp_box)))
    marks(upp_box) = "H"
    marks(upp_box)[ind] = "G"
    marks(upp_box) = as.factor(marks(upp_box))
    return(upp_box)
  }
}





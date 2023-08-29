#' Plot summary functions (developement version)
#' @param func Summary function to plot.
#' @param observed_values a list. observed summary function value. Must have one of the following structures
#' \itemize{
#'  \item{Option 1: List of summary functions}{observed_values[[func]][[cor]]}
#'  \item{Option 2: List of lists of summary functions}{observed_values[[pattern_num]][[func]][[cor]]}
#'  \item{Option 2: Dataframe}{observed_values[,c("r", y_value)]}
#'  }
#' @param envelopes a list envelope values found by
#' function such as   \code{\link{calc_summary_funcs}}.  Should be list of length equal
#' to the number of envelopes.  The structure should resemble:
#' `envelopes[[envelope_num]]$rrl_K[, c("r", "mmean")]`
#' @param pattern.colors colors and names for envelope lines.
#'  MUST follow some formatting rules.  Envelope names must come first.
#'  If each observed value corresponds to a different envelope, then
#'  the number of `observed_values` must match `length(envelopes)` and only
#'  the first `1:length(envelopes)` of `pattern_colors` will be used.
#'  and must set `base_value = "each"`.  If not, each envelope and each observed
#'  value needs its own color.  then the `mmean` value from
#'  `envelopes[[1]]` will be subtracted from each envelope and observed value.
#' @param fill.colors colors and names for envelope fill.  Match names to `pattern.colors`.
#' Recommended to leave as NA and it will automattically be set to match  `pattern.colors`
#' @param sqrt either `"K", "all", or "none"`. If `"K"`,
#' then only \code{expression(tilde(K)[g](r))} will be found using the differences of the
#' square roots, rather than differences.  If `"all"`, then all functions will be found
#' using such.  If `"none"` (or actually anything other than the first two options) then
#' no square roots are taken. Ignored if `raw = "TRUE"`
#' @param raw.  A logical.  If TRUE, then no envelope mmeans are subtracted from each value
#' @param base_value a character, either `"first" or "each"`.  Does each observed value
#' correspond to a different envelope?  If yes, set to `"each"`. Otherwise, set to
#' `"first"` and the envelopes[[1]] will be used for each
#' @param unit a character.  This will appear as the units in the x axis label
#' @param K_cor edge correction(s) to be used when `func = "K"`
#' @param G_cor edge correction(s) to be used when `func = "G"`
#' @param F_cor edge correction(s) to be used when `func = "F"`
#' @param GXGH_cor edge correction(s) to be used when `func = "GXGH"`
#' @param GXHG_cor edge correction(s) to be used when `func = "GXHG"`
#' @param alpha numeric between 0 and 1.  Transparency of envelopes
#' @param legend.key.size numeric.  size of legend key
#' @param legend.text.size numeric. size of legend text
#' @param legend.position.  vector of 2 numerics.  coordinates of legend
#' @param axis.title.size numeric.  Size of axis title
#' @param title a character.  Text for title
#' @param title.size numeric.  Title size
#' @param axis.text.x.size numeric.  size of text on x axis
#' @param axis.text.y.size numeric.  size of text on y axis
#' @param linewidth numeric.  Width of lines in plot
#' @param env_linewidth numeric.  Width of lines that make up envelope edges
#' @param linetype a character.  Type of lines that make up lines
#' @param env_linetype a character.  Type of lines that make up  envelope lines
#' @description Plot the observed value with envelopes for expected values for summary function
#' @details The best way to learn about this function is to read the parameter definitions.
#' @export
plot_summary_dev <- function(func = "K",
                              observed_values, envelopes, ...,
                              pattern.colors = c("Envelope 1" = "pink", "Envelope 2" = "gray", "Observed" = "blue"),
                              fill.colors =  NA,
                              sqrt = "K", raw = "FALSE",
                              base_value = "first", unit = "nm",
                              K_cor = "trans", G_cor = "km", F_cor = "km",
                              GXGH_cor = "km", GXHG_cor = "km", G2_cor = "km",
                              alpha = 0.5,
                              legend.key.size = 20, legend.text.size = 20,
                              legend.position =  c(0.75, 0.80),
                              axis.title.size = 40,
                              title = "", title.size = 40, axis.text.x.size = 40,
                              axis.text.y.size = 40, linewidth = 0.5,  env_linewidth = 0.5,
                              linetype = "solid",
                              env_linetype = "dashed") {
  if (all(is.na(fill.colors))) {
    fill.colors = pattern.colors
  }


  for (i in 2:10) {
    nm =paste("G", i, sep = "")
    if (func == nm) {
      assign(paste("G", i, "_cor", sep=""), G2_cor)
    }
  }

  cor = paste(func, "_cor", sep = "")

  cor = get(cor)

  rrl = paste("rrl_", func, sep = "")
  if (raw) {
    ylabs = list("K" = expression(K[g](r)), "G" = expression(G[g](r)),
                 "F" = expression(F[g](r)),
                 "GXGH" = expression(G[gh](r)), "GXHG" = expression(G[hg](r)),
                 "G2" = expression(G[2,g](r)), "G3" = expression(G[3,g](r)),
                 "G4" = expression(G[4,g](r)), "G5" = expression(G[5,g](r)),
                 "G6" = expression(G[6,g](r)), "G7" = expression(G[7,g](r)),
                 "G8" = expression(G[8,g](r)), "G9" = expression(G[9,g](r)),
                 "G10" = expression(G[10,g](r)))

  }
  else {
    ylabs = list("K" = expression(tilde(K)[g](r)), "G" = expression(tilde(G)[g](r)),
                 "F" = expression(tilde(F)[g](r)),
                 "GXGH" = expression(tilde(G)[gh](r)), "GXHG" = expression(tilde(G)[hg](r)),
                 "G2" = expression(tilde(G)[2,g](r)), "G3" = expression(tilde(G)[3,g](r)),
                 "G4" = expression(tilde(G)[4,g](r)), "G5" = expression(tilde(G)[5,g](r)),
                 "G6" = expression(tilde(G)[6,g](r)), "G7" = expression(tilde(G)[7,g](r)),
                 "G8" = expression(tilde(G)[8,g](r)), "G9" = expression(tilde(G)[9,g](r)),
                 "G10" = expression(tilde(G)[10,g](r)))
  }


  ## if returning the square root of the difference, rather than the difference
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
  # if there is an envelope for each observed value, then plot envelopes separately
  if ((length(envelopes) == length(observed_values)) && base_value != "first") {
    print("start each")
    long = data.frame("r" = c(),
                      "mmean" = c(),
                      "lo" = c(),
                      "hi" = c(),
                      "type" = c())
    i <- 1
    while (i <= length(envelopes)) {
      temp <- envelopes[[i]][[rrl]]
      observed <- observed_values[[i]][[func]]
      temp$type <- names(pattern.colors)[i]
      baseline = take_root(temp$mmean)

      if (raw) {
        baseline = 0
      }

      temp$lo = take_root(temp$lo) - baseline
      temp$hi = take_root(temp$hi) - baseline
      temp$mmean = take_root(observed[[cor]]) - baseline

      long <- rbind(long, temp)
      i <- i + 1
    }

    gplot <- long %>% ggplot2::ggplot(aes(
      x = r, ymin = lo,
      ymax = hi, color = type, fill = type
    )) +
      geom_ribbon(alpha = alpha, linewidth = env_linewidth, linetype = env_linetype) +
      geom_hline(yintercept = 0) +
      geom_line(aes(x= r, y = mmean, color = type), linewidth = linewidth, linetype = linetype)


    gplot <- gplot + theme(
      plot.title = element_text(hjust = 0.5, size = title.size),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = axis.text.x.size),
      axis.text.y = element_text(size = axis.text.y.size),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = axis.title.size),
      legend.key.size = unit(legend.key.size, "pt"),
      legend.text = element_text(size = legend.text.size),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position =legend.position,
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"),# ...
    ) +
      # guides(color = "none") +
      scale_color_manual(values = pattern.colors) +scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(ylabs[[func]]) +
      ggtitle(title)


  }

  #######

  else {
    print("start first")
    baseline = envelopes[[1]][[rrl]]$mmean
    if (raw) {
      baseline = 0
    }
    long = data.frame("r" = c(),
                      "mmean" = c(),
                      "lo" = c(),
                      "hi" = c(),
                      "type" = c())

    i <- 1
    while (i <= length(envelopes)) {
      temp <- envelopes[[i]][[rrl]]
      temp$type <- names(pattern.colors)[i]
      temp$lo = take_root(temp$lo) - take_root(baseline)
      temp$hi = take_root(temp$hi) - take_root(baseline)
      temp$mmean = take_root(temp$mmean) - take_root(baseline)
      long <- rbind(long, temp)
      i <- i + 1

    }

    #class(observed_funcs[[image_num]][[func]])
    if (is.data.frame(observed_values)) {
      observed <- data.frame(
        r = observed_values[["r"]],

        mmean = take_root(observed_values[[cor]]) - take_root(baseline),
        lo = take_root(observed_values[[cor]]) - take_root(baseline),
        hi = take_root(observed_values[[cor]]) - take_root(baseline),
        type = "Observed"
      )
    }


    else if (is.list(observed_values)) {
      if (is.data.frame(observed_values[[1]])) {
        observed = data.frame(row.names = c("r", "mmean", "lo", "hi", "type"))
        obs = observed_values[[func]]
        obs_01 <- data.frame(
          r = obs$r,

          mmean = take_root(obs[[cor]]) - take_root(baseline),
          lo = take_root(obs[[cor]]) - take_root(baseline),
          hi = take_root(obs[[cor]]) - take_root(baseline),
          type = names(pattern.colors)[length(envelopes) + 1])
        observed = rbind(observed, obs_01)
      }

      else if (is.list(observed_values[[1]])) {
        observed = data.frame(row.names = c("r", "mmean", "lo", "hi", "type"))
        for (i in 1:length(observed_values)) { #
          obs = (observed_values[[i]][[func]])
          obs_01 <- data.frame(
            r = obs$r,

            mmean = take_root(obs[[cor]]) - take_root(baseline),
            lo = take_root(obs[[cor]]) - take_root(baseline),
            hi = take_root(obs[[cor]]) - take_root(baseline),
            type = names(pattern.colors)[length(envelopes) + i])
          observed = rbind(observed, obs_01)
        }
      }


    }

    else {
      stop("observed_values must be either dataframe or list of dataframes with a single summary function,
             list of dataframes, each dataframe containing the same function,
             or a list of list of dataframes, the outer list being the pattern, the inner list being the different summary functions")
    }

    env_inds = 1:length(envelopes)
    obs_inds = (length(envelopes)+1):(length(pattern.colors))
    long <- rbind(long, observed)

    #long$type = as.factor(long$type)


    gplot =ggplot(long %>% filter(type %in% names(pattern.colors[env_inds])))+
      geom_ribbon(aes(x = r, ymin = (lo) ,
                      ymax = (hi) , color = type, fill = type ),
                  alpha = alpha, linewidth = env_linewidth, linetype = env_linetype) +
      geom_hline(yintercept = 0) +
      geom_line(data = long %>% filter(type %in% names(pattern.colors[obs_inds])),
                aes(x = r, y= (lo), color = type, fill =type),
                linewidth = linewidth, linetype = linetype)

    # gplot <- long %>% ggplot2::ggplot(aes(
    #  x = r, ymin = (lo) ,
    ##    ymax = (hi) , color = type, levels = c(), fill = type
    #  )) +
    #    geom_ribbon(alpha = alpha, linewidth = linewidth, linetype = linetype) +
    #    geom_hline(yintercept = 0)# +
    # geom_line(aes(x= r, y = sqrt(mmean) , color = type))



    print("start plot")
    gplot <- gplot + theme(
      plot.title = element_text(hjust = 0.5, size = title.size),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = axis.text.x.size),
      axis.text.y = element_text(size = axis.text.y.size),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = axis.title.size),
      legend.key.size = unit(legend.key.size, "pt"),
      legend.text = element_text(size = legend.text.size),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position =legend.position,
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"),# ...
    ) +
      # guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(ylabs[[func]]) +
      ggtitle(title)
  }
}



#' Mixed Multimer Simulation
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
#' @param weights list of vectors of length equal to number of dimensions in \emph{upp}. the weighted distance to each of \emph{num_neighbors} nearest neighbors
#' is calculated using \eqn{\sqrt{w_1 x^2 + w_2 y^2 + w_3 z^2}},
#'  where \emph{weights} = (\eqn{w_1, w_2, w_3}). Set to \emph{c(1, 1, 0)} for vertical dimers.
#' @param weights_fracs vector of numerics.  This will be the fraction of all guests that are assigned to the corresponding element of the `weights` list
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
#' \itemize{
#' \item{} {This will be the same as Case 1, except in Step 10 create a 3D point pattern (class `pp3`)
#' instead of a 2D point pattern}
#' }
#' @export
multimersim_mix = function(guest_pattern = NULL, upp, output = "guest pattern type", n_guests = NA,
                       min_thick = NA, max_thick = NA, ztrim = 0,
                       size_fracs = c(0, 1),
                       num_neighbors = 6,
                       sample_method = "rank", exponent = 1,
                       weights = list(c(1, 1, 1), c(1, 1, 0)),
                       weights_fracs = c(0.5, 0.5),
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

  ## get rid of any 0's at the end of the size_fracs vector
  for (i in 1:length(size_fracs)) {
    if (size_fracs[length(size_fracs)] == 0) {
      size_fracs = size_fracs[1:(length(size_fracs) - 1)]
    }
  }

  # rescale UPP pattern to match physical system (1 point per nm)
  upp_scaled = rTEM::rescale_pattern(upp, intensity_upp)

  # case 1 and 2:  guest_pattern is 2d
  if (spatstat.geom::is.ppp(guest_pattern) || output == "ppp") {

    if (is.null(guest_pattern) && spatstat.geom::is.pp3(upp)) {
      box_2d = as.owin(c(upp_scaled$domain$xrange,
                         upp_scaled$domain$yrange))

    }
    else if (is.null(guest_pattern) && spatstat.geom::is.ppp(upp)) {
      box_2d = upp_scaled$window
    }

    else {
      box_2d = guest_pattern$window
      n_guests = npoints(guest_pattern)
    }

    box_area = spatstat.geom::area(box_2d)

    # case 1 UPP pattern is 3d
    if (spatstat.geom::is.pp3(upp)) {
      print("Case 1: Multimers will be generated in 3D UPP and then collapsed into 2D")
      box_3d = domain(upp_scaled)

      if (is.na(max_thick)) {
        max_thick = max(upp_scaled$domain$zrange)
      }
      if (is.na(min_thick)) {
        min_thick = min(upp_scaled$domain$zrange)
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
      upp_labeled = rlabel(upp_box, labels = c(rep("A",n_groups - size_dist[1]),
                                               rep("B",size_dist[1]),
                                               rep("C", n_total - n_groups)), permute = TRUE)

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

        marks(upp_guest) = "G"
        marks(upp_host) = "H"


        multimer_coords = data.frame(x = c(upp_guest$data$x, upp_host$data$x),
                                     y = c(upp_guest$data$y, upp_host$data$y), z = c(upp_guest$data$z, upp_host$data$z),
                                     marks = c(upp_guest$data$marks, upp_host$data$marks))#, window = box_3d)
        chosen = multimer_coords[multimer_coords$z > min_thick+ztrim & multimer_coords$z  < max_thick - ztrim,]
        multimer_box = ppp(x = chosen$x, y = chosen$y, marks = as.factor(chosen$marks), window = box_2d)
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
      multimer_box = pp3(x = chosen$x, y = chosen$y, z =chosen$z, marks = as.factor(chosen$marks), window = box_3d_trimmed)

    }

  }
}



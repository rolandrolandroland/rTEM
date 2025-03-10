#' assemble png files into a png square deprecated Dec 20 2024
get_png_square_deprecated = function(input_path, save_path, pat, to_name, y_space = 50, x_space = 50, resolution = 70) {
  # Load PNG images

  path = paste(input_path, "mixed", sep = "/")

  mixed_files <- list.files(path = path, pattern = pat, full.names = TRUE)

  path = paste(input_path, "iso_dimer", sep = "/")
  iso_files <- list.files(path = path, pattern = pat, full.names = TRUE)

  path = paste(input_path, "vert_dimer", sep = "/")
  vert_files <- list.files(path = path, pattern = pat, full.names = TRUE)

  path = paste(input_path, "observed", sep = "/")
  observed_files = list.files(path = path, full.names = TRUE)

  mixed_images = lapply(mixed_files, readPNG)
  iso_images = lapply(iso_files, readPNG)
  vert_images = lapply(vert_files, readPNG)
  observed_images = lapply(observed_files, readPNG)

  f_ind = 1
  g_ind = 2
  g2_ind = 3
  g3_ind = 4
  k_ind = 5


  # Set space between images (in pixels)
  total_height = dim(observed_images[[k_ind]])[1] + dim(mixed_images[[k_ind]])[1] +
    dim(vert_images[[k_ind]])[1] + dim(iso_images[[k_ind]])[1] + (y_space * (ln -1))
  # Get total height and max width for the stacked image (including spaces between images)
  total_width = dim(iso_images[[k_ind]])[2] + dim(iso_images[[g_ind]])[2] +
    dim(iso_images[[g2_ind]])[2] + dim(iso_images[[f_ind]])[2] + (x_space * (3))


  # Open a PNG device with dimensions based on the total height and width
  # Adjust resolution with the `res` argument (e.g., 300 for high-resolution output)
  #png("stacked_tables.png", width = max_width, height = total_height, res = 72)
  png(paste(save_path, to_name, ".png", sep = ""), width = total_width , height = total_height, res = resolution)

  # Set up an empty plot with the right dimensions
  plot(NA, xlim = c(0, total_width), ylim = c(0, total_height), type = "n", xaxt = "n", yaxt = "n", bty = "n", xaxs = "i", yaxs = "i")

  # Initialize y_offset at the top (since plotting in R starts from bottom-left corner)
  #y_offset <- total_height
  y_offset = 0
  x_offset =0
  #x_offset = total_width
  ind_order = c(k_ind, g_ind, g2_ind, f_ind)

  images_to_use = c(observed_images[ind_order], iso_images[ind_order],
                    vert_images[ind_order], mixed_images[ind_order])
  list_images  = list(observed_images[ind_order], iso_images[ind_order],
                      vert_images[ind_order], mixed_images[ind_order])
  list_images = list(mixed_images[ind_order], vert_images[ind_order],
                     iso_images[ind_order], observed_images[ind_order])

  # Loop through the images and place them one after the other with space in between
  for (r in 1:length(list_images)) {
    x_offset = 0
    for (c in 1:length(ind_order)) {
      img_height = dim(list_images[[r]][[c]])[1]
      img_width = dim(list_images[[r]][[c]])[2]
      #y_offset <- y_offset - img_height
      #x_offset = x_offset - img_width
      rasterImage(list_images[[r]][[c]], x_offset, y_offset,
                  img_width + x_offset,
                  img_height + y_offset)
      x_offset = x_offset + x_space + img_width
    }
    y_offset = y_offset + y_space + img_height

  }
  dev.off()

}





#' Get Features deprecated Dec 19 2024
#' @param expected a list of summary functions for the observed value
#' @param simulated a list of
#' @param funcs list of summary functions to get features for.  currently `G`, `K`, and `F` are supported
#' @param features a list.  Must be either length 1 or same length as `funcs` parameter
#'  @description deprecated because the normalize total difference needed to only sum
#'  non zero env width
#'
get_features_deprecated = function(expected, simulated,
                        funcs = c("G", "K", "F"),
                        features = list(c("max_diff", "min_diff",
                                          "total_norm_diff", "total_norm_diff_squared",
                                          "total_diff", "total_squared_diff", "net_diff",
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

  for (i in 3:10) {
    nm =paste("G", i, sep = "")
    if (any(funcs == nm)) {
      assign(paste("G", i, "_cor", sep=""), G2_cor)
    }
  }
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
    current_env = environment()
    feats = features[[i]]
    cor = paste(func, "_cor", sep = "")
    cor = get(cor)#, envir= environment())
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
    if ("total_squared_diff" %in% feats) {
      total_squared_diff = sum(abs(take_root(expected[[func]][[cor]]) - take_root(simulated[[rrl]][,simulated_inner_name[2]]))^2) * dr
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


#' Plot summary functions (deprecated 240911)
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
plot_summary_240911 <- function(func = "K",
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
  ylabs = list("K" = expression(tilde(K)[g](r)), "G" = expression(tilde(G)[g](r)),
               "F" = expression(tilde(F)[g](r)),
               "GXGH" = expression(tilde(G)[gh](r)), "GXHG" = expression(tilde(G)[hg](r)),
               "G2" = expression(tilde(G)[2,g](r)), "G3" = expression(tilde(G)[3,g](r)),
               "G4" = expression(tilde(G)[4,g](r)), "G5" = expression(tilde(G)[5,g](r)),
               "G6" = expression(tilde(G)[6,g](r)), "G7" = expression(tilde(G)[7,g](r)),
               "G8" = expression(tilde(G)[8,g](r)), "G9" = expression(tilde(G)[9,g](r)),
               "G10" = expression(tilde(G)[10,g](r)))


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


    long <- rbind(long, observed)

    gplot <- long %>% ggplot2::ggplot(aes(
      x = r, ymin = (lo) ,
      ymax = (hi) , color = type, fill = type
    )) +
      geom_ribbon(alpha = alpha, linewidth = env_linewidth, linetype = env_linetype) +
      geom_line(aes(x = r, y = lo, color = type, fill = type),
                linewidth = linewidth, linetype = env_linetype,
                data = filter(long, type == "Observed") ) +
      geom_hline(yintercept = 0)# +
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


#' Find Points
#' @param OVER overlying point pattern
findNearestPoints_dep <- function(OVER, UNDER, n = 1) {
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

#' Calculate T value
#'
#' @param relabelings object from \code{\link[rTEM]{relabel_summarize}}
#' @param func Which function is used: currently on G and K are supported
#' @param rmin what value to start the calculation from
#' @param rmax max value to run calculation to
#' @param K_cor boundary correction type for K
#' @param G_cor boundary correction type for G
#'
#' @description
#' Takes input of relabeling object and calculates the T value for either K or G function
#'
#' @details This uses the method from \emph{Baddeley et al.} to calculate the T value for an envelope for either
#' K  (\emph{K3est}) G \emph{G3est} function
#' @references
#' Baddeley, A., Diggle, P. J., Hardegen, A., Lawrence, T_, Milne, R. K., & Nair, G. (2014).
#' On tests of spatial pattern based on simulation envelopes. Ecological Monographs, 84(3), 477–489. https://doi.org/10.1890/13-2042.1
#' @export
calc_T_vals_v01 <- function(relabelings, func = "K", rmin = 0, rmax, K_cor = "border", G_cor = "km") {
  if (func == "K") {
    if (K_cor == "trans") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_K$trans
      })
    } else if (K_cor == "iso") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_K$iso
      })
    } else if (K_cor == "border") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_K$bord
      })
    } else {
      print("Incorrect K edge correction")
    }

    r <- relabelings[[1]]$rrl_K$r
    med <- apply(mmean, 1, median)
    ind <- which(r <= rmax & r >= rmin)
    interval <- rmax / (length(r) - 1)
    T <- apply(mmean, 2, function(x) {
      T_1 <- sqrt(x) - sqrt(med)
      T_1 <- T_1[ind]
      sum(T_1^2) * interval
    })
    return(T)
  }
  if (func == "G") {
    if (G_cor == "km") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_G$km
      })
    } else if (G_cor == "rs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_G$rs
      })
    } else {
      print("Incorrect G edge correction")
    }
    r <- relabelings[[1]]$rrl_G$r
    med <- apply(mmean, 1, median)
    ind <- which(r <= rmax & r >= rmin)
    interval <- rmax / (length(r) - 1)
    T <- apply(mmean, 2, function(x) {
      T_1 <- x - med
      T_1 <- T_1[ind]
      sum(T_1^2) * interval
    })
  }
  return(T)

  if (func == "F") {
    if (F_cor == "km") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_F$km
      })
    } else if (F_cor == "rs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_F$rs
      })
    } else if (F_cor == "cs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_F$cs
      })
    } else {
      print("Incorrect F edge correction")
    }
    r <- relabelings[[1]]$rrl_F$r
    med <- apply(mmean, 1, median)
    ind <- which(r <= rmax & r >= rmin)
    interval <- rmax / (length(r) - 1)
    T <- apply(mmean, 2, function(x) {
      T_1 <- x - med
      T_1 <- T_1[ind]
      sum(T_1^2) * interval
    })
  }
  return(T)

  if (func == "GXGH") {
    if (GXGH_cor == "km") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXGH$km
      })
    } else if (GXGH_cor == "rs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXGH$rs
      })
    } else if (GXGH_cor == "han") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXGH$han
      })
    } else if (GXGH_cor == "none") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXGH$raw
      })
    } else {
      print("Incorrect GXGH edge correction")
    }
    r <- relabelings[[1]]$rrl_GXGH$r
    med <- apply(mmean, 1, median)
    ind <- which(r <= rmax & r >= rmin)
    interval <- rmax / (length(r) - 1)
    T <- apply(mmean, 2, function(x) {
      T_1 <- x - med
      T_1 <- T_1[ind]
      sum(T_1^2) * interval
    })
  }
  return(T)


  if (func == "GXHG") {
    if (GXHG_cor == "km") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXHG$km
      })
    } else if (GXHG_cor == "rs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXHG$rs
      })
    } else if (GXHG_cor == "han") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXHG$han
      })
    } else if (GXHG_cor == "none") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXHG$raw
      })
    } else {
      print("Incorrect GXHG edge correction")
    }
    r <- relabelings[[1]]$rrl_GXHG$r
    med <- apply(mmean, 1, median)
    ind <- which(r <= rmax & r >= rmin)
    interval <- rmax / (length(r) - 1)
    T <- apply(mmean, 2, function(x) {
      T_1 <- x - med
      T_1 <- T_1[ind]
      sum(T_1^2) * interval
    })
  }
  return(T)
}



#' Average relabelings
#' @param relabelings output of  \code{\link[rTEM]{relabel_summarize}} function
#' @param envelope.value size of envelope to compute.  Should be decimal (e.g. 0.95 = 95\%)
#' @param funcs vector of summary functions to calculate
#' @param K_cor edge correction(s) to be used for K function
#' @param G_cor edge correction(s) to be used for G function
#' @param F_cor edge correction(s) to be used for F function
#' @param GXGH_cor edge correction(s) to be used for GXGH function
#' @param GXHG_cor edge correction(s) to be used for GXHG function
#'
#' @description
#' Function take all relabelings and return averages and envelope values
#'
#' @details Take output of  \code{\link[rTEM]{relabel_summarize}} and create envelopes
#'
#' @return data frame with x value (distance), average y value, and y value envelopess
#' @export
average_relabelings_v01 <- function(relabelings, envelope.value = .95,
                                funcs = c("K", "G"),
                                K_cor = "trans", G_cor = "km", F_cor = "km",
                                GXGH_cor = "km", GXHG_cor = "km") {
  # transform envelope value to high index (0.95 envelope will be 0.025 and 0.975)
  envelope.value <- envelope.value + (1 - envelope.value) / 2

  # get index of high and low envelope values and find values at each index
  hi.ind <- round((length(relabelings) + 1) * envelope.value, 0)
  lo.ind <- round((length(relabelings) + 1) * (1 - envelope.value), 0)
  if (lo.ind == 0) {
    lo.ind <- 1
  }
  # make the relabelings their own individual objects

  # K
  # extract K(r) values
  if ("K" %in% funcs) {
    if (K_cor == "trans") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_K$trans
      })
    } else if (K_cor == "iso") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_K$iso
      })
    } else if (K_cor == "border") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_K$bord
      })
    } else {
      print("Incorrect K edge correction")
    }

    # order K(r) values by value
    ordered <- apply(mmean, 1, sort)

    lo <- ordered[lo.ind, ]
    hi <- ordered[hi.ind, ]

    # get r values
    r <- relabelings[[1]]$rrl_K$r

    # find the median at every distance
    med <- apply(mmean, 1, median)
    rrl_K <- data.frame(r = r, mmean = med, lo = lo, hi = hi)
  }

  # Repeat for G, GX, and F
  # G
  if ("G" %in% funcs) {
    if (G_cor == "km") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_G$km
      })
    } else if (G_cor == "rs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_G$rs
      })
    } else {
      print("Incorrect G edge correction")
    }

    r <- relabelings[[1]]$rrl_G$r
    ordered <- apply(mmean, 1, sort)
    lo <- ordered[lo.ind, ]
    hi <- ordered[hi.ind, ]
    med <- apply(mmean, 1, median)
    rrl_G <- data.frame(r = r, mmean = med, lo = lo, hi = hi)
  }

  # F
  if ("F" %in% funcs) {
    if (F_cor == "km") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_F$km
      })
    } else if (F_cor == "rs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_F$rs
      })
    } else if (F_cor == "cs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_F$cs
      })
    } else {
      print("Incorrect F edge correction")
    }
    r <- relabelings[[1]]$rrl_F$r
    ordered <- apply(mmean, 1, sort)
    lo <- ordered[lo.ind, ]
    hi <- ordered[hi.ind, ]
    med <- apply(mmean, 1, median)
    rrl_F <- data.frame(r = r, mmean = med, lo = lo, hi = hi)
  }

  # GXGH
  if ("GXGH" %in% funcs) {
    if (GXGH_cor == "km") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXGH$km
      })
    } else if (GXGH_cor == "rs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXGH$rs
      })
    } else if (GXGH_cor == "han") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXGH$han
      })
    } else if (GXGH_cor == "none") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXGH$raw
      })
    } else {
      print("Incorrect GXGH edge correction")
    }

    r <- relabelings[[1]]$rrl_GXGH$r
    ordered <- apply(mmean, 1, sort)
    lo <- ordered[lo.ind, ]
    hi <- ordered[hi.ind, ]
    med <- apply(mmean, 1, median)
    rrl_GXGH <- data.frame(r = r, mmean = med, lo = lo, hi = hi)
  }

  # GXGH
  if ("GXHG" %in% funcs) {
    if (GXHG_cor == "km") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXHG$km
      })
    } else if (GXHG_cor == "rs") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXHG$rs
      })
    } else if (GXHG_cor == "han") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXHG$han
      })
    } else if (GXHG_cor == "none") {
      mmean <- sapply(relabelings, function(x) {
        x$rrl_GXHG$raw
      })
    } else {
      print("Incorrect GXHG edge correction")
    }

    r <- relabelings[[1]]$rrl_GXHG$r
    ordered <- apply(mmean, 1, sort)
    lo <- ordered[lo.ind, ]
    hi <- ordered[hi.ind, ]
    med <- apply(mmean, 1, median)
    rrl_GXHG <- data.frame(r = r, mmean = med, lo = lo, hi = hi)
  }
  all_funcs <- c("K", "G", "F", "GXGH", "GXHG")
  all_funcs %in% funcs
  relabs <- c("rrl_K", "rrl_G", "rrl_F", "rrl_GXGH", "rrl_GXHG")
  out <- lapply(relabs[all_funcs %in% funcs], get, envir = environment())
  names(out) <- relabs[all_funcs %in% funcs]
  return(out)
}







## these file are deprecated:: I am keeping copies of them, but they are not to be used any more
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
  rcp_scaled = rTEM::rescale_pattern(rcp_pattern, intensity_rcp)

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
  rcp_scaled = rTEM::rescale_pattern(rcp_pattern, intensity_rcp)

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
  upp_scaled = rTEM::rescale_pattern(upp, intensity_upp)

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
  upp_scaled = rTEM::rescale_pattern(upp, intensity_upp)

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






#' 2D relabeling function
relabeler_2d_deprecated = function(relabelings, envelope_value = .95) {
  mmean =sapply(relabelings, function(x) {
    x$rrl_K$border
  })
  envelope_value =  envelope_value + (1- envelope_value)/2

  ordered = apply(mmean, 1, sort)
  hi_ind = round(length(relabelings) * envelope_value, 0)
  lo_ind = round(length(relabelings) * (1-envelope_value), 0)
  if (lo_ind == 0) {
    lo_ind = 1
  }
  lo = ordered[lo_ind,]
  min = ordered[1,]
  hi = ordered[hi_ind,]
  max = ordered[nrow(ordered),]
  r = relabelings[[1]]$rrl_K$r
  med = apply(mmean, 1, median)
  rrl_K = data.frame(r = r, mmean =med, lo = lo, hi = hi,
                     min = min, max = max,
                     lo_diff = lo - med,
                     hi_diff = hi - med,
                     min_diff = min - med,
                     max_diff = max - med,
                     lo_sqrt_diff = sqrt(lo) - sqrt(med),
                     hi_sqrt_diff = sqrt(hi) - sqrt(med),
                     min_sqrt_diff = sqrt(min) - sqrt(med),
                     max_sqrt_diff = sqrt(max) - sqrt(med))


  # G
  mmean =sapply(relabelings, function(x) {
    x$rrl_G$km
  })
  r = relabelings[[1]]$rrl_G$r
  ordered = apply(mmean, 1, sort)
  lo = ordered[lo_ind,]
  min = ordered[1,]
  hi = ordered[hi_ind,]
  max = ordered[length(relabelings),]
  med = apply(mmean, 1, median)
  rrl_G = data.frame(r = r, mmean =med, lo = lo, hi = hi,
                     min = min, max = max,
                     lo_diff = lo - med,
                     hi_diff = hi - med,
                     min_diff = min - med,
                     max_diff = max - med)

  return(list(rrl_G = rrl_G, rrl_K =rrl_K))
}

#' Plot Summary functions
#'
plot_summary_v01 <- function(func = "K",
                         observed_values, envelopes, ...,
                         pattern.colors = c("95.0% AI" = "pink", "Median" = "black", "Observed" = "blue"),
                         fill.colors = pattern.colors, unit = "nm",
                         K_cor = "trans", G_cor = "km", F_cor = "km",
                         GXGH_cor = "km", GXHG_cor = "km") {
  if (func == "K") {
    long <- as.data.frame(envelopes[[1]]$rrl_K)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_K)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }

    baseline <- sqrt(envelopes[[1]]$rrl_K$mmean)

    observed <- as.data.frame(observed_values$K)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, K_cor],
      lo = observed[, K_cor], hi = observed[, K_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = sqrt(lo) - baseline,
      ymax = sqrt(hi) - baseline, color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(K)[g](r)))
  }

  ## plot for G
  else if (func == "G") {
    long <- as.data.frame(envelopes[[1]]$rrl_G)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    head(long)
    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_G)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }

    baseline <- (envelopes[[1]]$rrl_G$mmean)

    observed <- as.data.frame(observed_values$G)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, G_cor],
      lo = observed[, G_cor], hi = observed[, G_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo) - baseline,
      ymax = (hi) - baseline, color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(G)[g](r)))
  }


  ### Plot for F
  else if (func == "F") {
    long <- as.data.frame(envelopes[[1]]$rrl_F)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_F)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }

    baseline <- (envelopes[[1]]$rrl_F$mmean)

    observed <- as.data.frame(observed_values$F)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, F_cor],
      lo = observed[, F_cor], hi = observed[, F_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo) - baseline,
      ymax = (hi) - baseline, color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(F)[g](r)))
  }

  ### Plot for GXGH
  else if (func == "GXGH") {
    long <- as.data.frame(envelopes[[1]]$rrl_GXGH)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    head(long)
    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_GXGH)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }

    baseline <- (envelopes[[1]]$rrl_GXGH$mmean)

    observed <- as.data.frame(observed_values$GXGH)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, GXGH_cor],
      lo = observed[, GXGH_cor], hi = observed[, GXGH_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo) - baseline,
      ymax = (hi) - baseline, color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(G)[gh](r)))
  }

  ### Plot for GXHG
  else if (func == "GXHG") {
    long <- as.data.frame(envelopes[[1]]$rrl_GXHG)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    head(long)
    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_GXHG)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }

    baseline <- (envelopes[[1]]$rrl_GXHG$mmean)

    observed <- as.data.frame(observed_values$GXHG)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, GXHG_cor],
      lo = observed[, GXHG_cor], hi = observed[, GXHG_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo) - baseline,
      ymax = (hi) - baseline, color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(G)[hg](r)))
  }
}

#' Plot summary functions (version 2, deprecated)
#' @param func Summary function to plot
#' @param observed_values a list. observed summary function value
#' @param envelopes a list theoretical values found by function such as `calc_summary_funcs().` Should be list
#' @param pattern.colors colors and names for envelope lines
#' @param fill.colors colors and names for envelope fill.  Match names to `pattern.colors`
#' @param K_cor edge correction(s) to be used for K function
#' @param G_cor edge correction(s) to be used for G function
#' @param F_cor edge correction(s) to be used for F function
#' @param GXGH_cor edge correction(s) to be used for GXGH function
#' @param GXHG_cor edge correction(s) to be used for GXHG function
#'
#' @description Under developement. Plot the observed value with envelopes for expected values for summary function
#' @export
plot_summary_v02 <- function(func = "K",
                             observed_values, envelopes, ...,
                             pattern.colors = c("95.0% AI" = "pink", "Median" = "black", "Observed" = "blue"),
                             fill.colors = pattern.colors,
                             base_value = "first", unit = "nm",
                             K_cor = "trans", G_cor = "km", F_cor = "km",
                             GXGH_cor = "km", GXHG_cor = "km",
                             legend.key.size = 20, legend.text.size = 20,
                             legend.position =  c(0.75, 0.80),
                             axis.title.size = 40,
                             title = "", title.size = 40, axis.text.x.size = 40, axis.text.y.size = 40,
                             alpha = 0.5, linewidth = 0.5,  env_linewidth = 0.5, linetype = "solid",
                             env_linetype = "dashed") {
  if (func == "K") {
    print("start K")


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
        temp <- as.data.frame(envelopes[[i]]$rrl_K)
        observed <- as.data.frame(observed_values[[i]]$K)
        #temp$type <- as.factor(i)
        temp$type <- names(pattern.colors)[i]
        #temp$r = observed_values[[i]]$K$r
        #temp$mmean = observed_values[[i]]$K[,K_cor]
        temp$lo = sqrt(temp$lo) - sqrt(temp$mmean)
        temp$hi = sqrt(temp$hi) - sqrt(temp$mmean)
        temp$mmean = sqrt(observed[, K_cor]) - sqrt(temp$mmean)
        long <- rbind(long, temp)
        i <- i + 1
      }

      gplot <- long %>% ggplot(aes(
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
        xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(K)[g](r))) +
        ggtitle(title)


    }

    #######

    else {
      print("start first")
      baseline = envelopes[[1]]$rrl_K$mmean
      long = data.frame("r" = c(),
                        "mmean" = c(),
                        "lo" = c(),
                        "hi" = c(),
                        "type" = c())

      i <- 1
      while (i <= length(envelopes)) {
        temp <- as.data.frame(envelopes[[i]]$rrl_K)
        temp$type <- names(pattern.colors)[i]
        temp$lo = sqrt(temp$lo) - sqrt(baseline)
        temp$hi = sqrt(temp$hi) - sqrt(baseline)
        temp$mmean = sqrt(temp$mmean) - sqrt(baseline)
        long <- rbind(long, temp)
        i <- i + 1

      }

      #
      if (is.data.frame(observed_values)) {
        observed <- as.data.frame(observed_values$K)
        observed <- data.frame(
          r = observed$r,
          mmean = sqrt(observed[, K_cor]) - sqrt(baseline),
          lo = sqrt(observed[, K_cor]) - sqrt(baseline),
          hi = sqrt(observed[, K_cor]) - sqrt(baseline),
          type = "Observed"
        )
      }
      else if (is.list(observed_values)) {
        observed = data.frame(row.names = c("r", "mmean", "lo", "hi", "type"))
        for (i in 1:length(observed_values)) { #
          obs = as.data.frame(observed_values[[i]]$K)
          obs_01 <- data.frame(
            r = obs$r,
            mmean = sqrt(obs[, K_cor]) - sqrt(baseline),
            lo = sqrt(obs[, K_cor]) - sqrt(baseline),
            hi = sqrt(obs[, K_cor]) - sqrt(baseline),
            type = names(pattern.colors)[length(envelopes) + i])
          observed = rbind(observed, obs_01)
        }
      }

      else {
        stop("observed_values must be either dataframe or list of dataframes")
      }


      long <- rbind(long, observed)

      gplot <- long %>% ggplot(aes(
        x = r, ymin = (lo) ,
        ymax = (hi) , color = type, fill = type
      )) +
        geom_ribbon(alpha = alpha, linewidth = linewidth, linetype = linetype) +
        geom_hline(yintercept = 0)# +
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
        xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(K)[g](r))) +
        ggtitle(title)
    }
  }


  ### if plotting G
  else if (func == "G") {


    # if there is an envelope for each observed value, then plot envelopes separately
    if ((length(envelopes) == length(observed_values)) && base_value != "first") {
      long = data.frame("r" = c(),
                        "mmean" = c(),
                        "lo" = c(),
                        "hi" = c(),
                        "type" = c())
      i <- 1
      while (i <= length(envelopes)) {
        temp <- as.data.frame(envelopes[[i]]$rrl_G)
        observed <- as.data.frame(observed_values[[i]]$G)
        #temp$type <- as.factor(i)
        temp$type <- names(pattern.colors)[i]
        #temp$r = observed_values[[i]]$G$r
        #temp$mmean = observed_values[[i]]$G[,G_cor]
        temp$lo = (temp$lo) - (temp$mmean)
        temp$hi = (temp$hi) - (temp$mmean)
        temp$mmean = (observed[, G_cor]) - (temp$mmean)
        long <- rbind(long, temp)
        i <- i + 1
      }

      gplot <- long %>% ggplot(aes(
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
        xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(G)[g](r))) +
        ggtitle(title)
    }

    #######

    else {
      baseline = envelopes[[1]]$rrl_G$mmean
      long = data.frame("r" = c(),
                        "mmean" = c(),
                        "lo" = c(),
                        "hi" = c(),
                        "type" = c())

      i <- 1
      while (i <= length(envelopes)) {
        temp <- as.data.frame(envelopes[[i]]$rrl_G)
        temp$type <- names(pattern.colors)[i]
        temp$lo = (temp$lo) - (baseline)
        temp$hi = (temp$hi) - (baseline)
        temp$mmean = (temp$mmean) - (baseline)
        long <- rbind(long, temp)
        i <- i + 1

      }

      #
      if (is.data.frame(observed_values)) {
        observed <- as.data.frame(observed_values$G)
        observed <- data.frame(
          r = observed$r,
          mmean = (observed[, G_cor]) - (baseline),
          lo = (observed[, G_cor]) - (baseline),
          hi = (observed[, G_cor]) - (baseline),
          type = "Observed"
        )
      }
      else if (is.list(observed_values)) {
        observed = data.frame(row.names = c("r", "mmean", "lo", "hi", "type"))
        for (i in 1:length(observed_values)) { #
          obs = as.data.frame(observed_values[[i]]$G)
          obs_01 <- data.frame(
            r = obs$r,
            mmean = (obs[, G_cor]) - (baseline),
            lo = (obs[, G_cor]) - (baseline),
            hi = (obs[, G_cor]) - (baseline),
            type = names(pattern.colors)[length(envelopes) + i])
          observed = rbind(observed, obs_01)
        }
      }

      else {
        stop("observed_values must be either dataframe or list of dataframes")
      }


      long <- rbind(long, observed)

      gplot <- long %>% ggplot(aes(
        x = r, ymin = (lo) ,
        ymax = (hi) , color = type, fill = type
      )) +
        geom_ribbon(alpha = alpha, linewidth = linewidth, linetype = linetype) +
        geom_hline(yintercept = 0)



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
        xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(G)[g](r))) +
        ggtitle(title)
    }
  }

  ### if plotting F
  else if (func == "F") {


    # if there is an envelope for each observed value, then plot envelopes separately
    if ((length(envelopes) == length(observed_values)) && base_value != "first") {
      long = data.frame("r" = c(),
                        "mmean" = c(),
                        "lo" = c(),
                        "hi" = c(),
                        "type" = c())
      i <- 1
      while (i <= length(envelopes)) {
        temp <- as.data.frame(envelopes[[i]]$rrl_F)
        observed <- as.data.frame(observed_values[[i]]$F)
        #temp$type <- as.factor(i)
        temp$type <- names(pattern.colors)[i]
        #temp$r = observed_values[[i]]$F$r
        #temp$mmean = observed_values[[i]]$F[,F_cor]
        temp$lo = (temp$lo) - (temp$mmean)
        temp$hi = (temp$hi) - (temp$mmean)
        temp$mmean = (observed[, F_cor]) - (temp$mmean)
        long <- rbind(long, temp)
        i <- i + 1
      }

      gplot <- long %>% ggplot(aes(
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
        xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(F)[g](r))) +
        ggtitle(title)
    }

    #######

    else {
      baseline = envelopes[[1]]$rrl_F$mmean
      long = data.frame("r" = c(),
                        "mmean" = c(),
                        "lo" = c(),
                        "hi" = c(),
                        "type" = c())

      i <- 1
      while (i <= length(envelopes)) {
        temp <- as.data.frame(envelopes[[i]]$rrl_F)
        temp$type <- names(pattern.colors)[i]
        temp$lo = (temp$lo) - (baseline)
        temp$hi = (temp$hi) - (baseline)
        temp$mmean = (temp$mmean) - (baseline)
        long <- rbind(long, temp)
        i <- i + 1

      }

      #
      if (is.data.frame(observed_values)) {
        observed <- as.data.frame(observed_values$F)
        observed <- data.frame(
          r = observed$r,
          mmean = (observed[, F_cor]) - (baseline),
          lo = (observed[, F_cor]) - (baseline),
          hi = (observed[, F_cor]) - (baseline),
          type = "Observed"
        )
      }
      else if (is.list(observed_values)) {
        observed = data.frame(row.names = c("r", "mmean", "lo", "hi", "type"))
        for (i in 1:length(observed_values)) { #
          obs = as.data.frame(observed_values[[i]]$F)
          obs_01 <- data.frame(
            r = obs$r,
            mmean = (obs[, F_cor]) - (baseline),
            lo = (obs[, F_cor]) - (baseline),
            hi = (obs[, F_cor]) - (baseline),
            type = names(pattern.colors)[length(envelopes) + i])
          observed = rbind(observed, obs_01)
        }
      }

      else {
        stop("observed_values must be either dataframe or list of dataframes")
      }


      long <- rbind(long, observed)

      gplot <- long %>% ggplot(aes(
        x = r, ymin = (lo) ,
        ymax = (hi) , color = type, fill = type
      )) +
        geom_ribbon(alpha = alpha, linewidth = linewidth, linetype = linetype) +
        geom_hline(yintercept = 0)


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
        xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(F)[G](r))) +
        ggtitle(title)
    }
  }

  ### if plotting GXGH
  else if (func == "GXGH") {


    # if there is an envelope for each observed value, then plot envelopes separately
    if ((length(envelopes) == length(observed_values)) && base_value != "first") {
      long = data.frame("r" = c(),
                        "mmean" = c(),
                        "lo" = c(),
                        "hi" = c(),
                        "type" = c())
      i <- 1
      while (i <= length(envelopes)) {
        temp <- as.data.frame(envelopes[[i]]$rrl_GXGH)
        observed <- as.data.frame(observed_values[[i]]$GXGH)
        #temp$type <- as.factor(i)
        temp$type <- names(pattern.colors)[i]
        #temp$r = observed_values[[i]]$GXGH$r
        #temp$mmean = observed_values[[i]]$GXGH[,GXGH_cor]
        temp$lo = (temp$lo) - (temp$mmean)
        temp$hi = (temp$hi) - (temp$mmean)
        temp$mmean = (observed[, GXGH_cor]) - (temp$mmean)
        long <- rbind(long, temp)
        i <- i + 1
      }

      gplot <- long %>% ggplot(aes(
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
        xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(G)[gh](r))) +
        ggtitle(title)
    }

    #######

    else {
      baseline = envelopes[[1]]$rrl_GXGH$mmean
      long = data.frame("r" = c(),
                        "mmean" = c(),
                        "lo" = c(),
                        "hi" = c(),
                        "type" = c())

      i <- 1
      while (i <= length(envelopes)) {
        temp <- as.data.frame(envelopes[[i]]$rrl_GXGH)
        temp$type <- names(pattern.colors)[i]
        temp$lo = (temp$lo) - (baseline)
        temp$hi = (temp$hi) - (baseline)
        temp$mmean = (temp$mmean) - (baseline)
        long <- rbind(long, temp)
        i <- i + 1

      }

      #
      if (is.data.frame(observed_values)) {
        observed <- as.data.frame(observed_values$GXGH)
        observed <- data.frame(
          r = observed$r,
          mmean = (observed[, GXGH_cor]) - (baseline),
          lo = (observed[, GXGH_cor]) - (baseline),
          hi = (observed[, GXGH_cor]) - (baseline),
          type = "Observed"
        )
      }
      else if (is.list(observed_values)) {
        observed = data.frame(row.names = c("r", "mmean", "lo", "hi", "type"))
        for (i in 1:length(observed_values)) { #
          obs = as.data.frame(observed_values[[i]]$GXGH)
          obs_01 <- data.frame(
            r = obs$r,
            mmean = (obs[, GXGH_cor]) - (baseline),
            lo = (obs[, GXGH_cor]) - (baseline),
            hi = (obs[, GXGH_cor]) - (baseline),
            type = names(pattern.colors)[length(envelopes) + i])
          observed = rbind(observed, obs_01)
        }
      }

      else {
        stop("observed_values must be either dataframe or list of dataframes")
      }


      long <- rbind(long, observed)

      gplot <- long %>% ggplot(aes(
        x = r, ymin = (lo) ,
        ymax = (hi) , color = type, fill = type
      )) +
        geom_ribbon(alpha = alpha, linewidth = env_linewidth, linetype = env_linetype) +
        geom_hline(yintercept = 0)



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
        xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(G)[gh](r))) +
        ggtitle(title)
    }
  }

  # Plot for GXHG
  else if (func == "GXHG")
  {
    # if there is an envelope for each observed value, then plot envelopes separately
    if ((length(envelopes) == length(observed_values)) && base_value != "first") {
      long = data.frame("r" = c(),
                        "mmean" = c(),
                        "lo" = c(),
                        "hi" = c(),
                        "type" = c())
      i <- 1
      while (i <= length(envelopes)) {
        temp <- as.data.frame(envelopes[[i]]$rrl_GXHG)
        observed <- as.data.frame(observed_values[[i]]$GXHG)
        #temp$type <- as.factor(i)
        temp$type <- names(pattern.colors)[i]
        #temp$r = observed_values[[i]]$GXHG$r
        #temp$mmean = observed_values[[i]]$GXHG[,GXHG_cor]
        temp$lo = (temp$lo) - (temp$mmean)
        temp$hi = (temp$hi) - (temp$mmean)
        temp$mmean = (observed[, GXHG_cor]) - (temp$mmean)
        long <- rbind(long, temp)
        i <- i + 1
      }

      gplot <- long %>% ggplot(aes(
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
        xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(G)[hg](r))) +
        ggtitle(title)
    }

    #######

    else {
      baseline = envelopes[[1]]$rrl_GXHG$mmean
      long = data.frame("r" = c(),
                        "mmean" = c(),
                        "lo" = c(),
                        "hi" = c(),
                        "type" = c())

      i <- 1
      while (i <= length(envelopes)) {
        temp <- as.data.frame(envelopes[[i]]$rrl_GXHG)
        temp$type <- names(pattern.colors)[i]
        temp$lo = (temp$lo) - (baseline)
        temp$hi = (temp$hi) - (baseline)
        temp$mmean = (temp$mmean) - (baseline)
        long <- rbind(long, temp)
        i <- i + 1

      }

      #
      if (is.data.frame(observed_values)) {
        observed <- as.data.frame(observed_values$GXHG)
        observed <- data.frame(
          r = observed$r,
          mmean = (observed[, GXHG_cor]) - (baseline),
          lo = (observed[, GXHG_cor]) - (baseline),
          hi = (observed[, GXHG_cor]) - (baseline),
          type = "Observed"
        )
      }
      else if (is.list(observed_values)) {
        observed = data.frame(row.names = c("r", "mmean", "lo", "hi", "type"))
        for (i in 1:length(observed_values)) { #
          obs = as.data.frame(observed_values[[i]]$GXHG)
          obs_01 <- data.frame(
            r = obs$r,
            mmean = (obs[, GXHG_cor]) - (baseline),
            lo = (obs[, GXHG_cor]) - (baseline),
            hi = (obs[, GXHG_cor]) - (baseline),
            type = names(pattern.colors)[length(envelopes) + i])
          observed = rbind(observed, obs_01)
        }
      }

      else {
        stop("observed_values must be either dataframe or list of dataframes")
      }


      long <- rbind(long, observed)

      gplot <- long %>% ggplot(aes(
        x = r, ymin = (lo) ,
        ymax = (hi) , color = type, fill = type
      )) +
        geom_ribbon(alpha = alpha, linewidth = env_linewidth, linetype = env_linetype) +
        geom_hline(yintercept = 0)

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
        legend.background = element_rect(fill = "white", color = "black"), ...
      ) +
        # guides(color = "none") +
        scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
        xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(tilde(G)[hg](r))) +
        ggtitle(title)
    }
  }


}


#' Plot summary functions raw V01 (Deprecated)
#' @param func Summary function to plot
#' @param observed_values observed summary function value
#' @param envelopes theoretical values found by function such as `calc_summary_funcs().` Should be list
#' @param pattern.colors colors and names for envelope lines
#' @param fill.colors colors and names for envelope fill.  Match names to `pattern.colors`
#' @param K_cor edge correction(s) to be used for K function
#' @param G_cor edge correction(s) to be used for G function
#' @param F_cor edge correction(s) to be used for F function
#' @param GXGH_cor edge correction(s) to be used for GXGH function
#' @param GXHG_cor edge correction(s) to be used for GXHG function
#'
#' @description Plot the observed value with envelopes for expected values for summary function
#' @export
plot_summary_raw_v01 <- function(func = "K",
                             observed_values, envelopes, ...,
                             pattern.colors = c("95.0% AI" = "pink", "Median" = "black", "Observed" = "blue"),
                             fill.colors = pattern.colors, unit = "nm",
                             K_cor = "trans", G_cor = "km", F_cor = "km",
                             GXGH_cor = "km", GXHG_cor = "km") {
  if (func == "K") {
    long <- as.data.frame(envelopes[[1]]$rrl_K)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_K)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }


    observed <- as.data.frame(observed_values$K)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, K_cor],
      lo = observed[, K_cor], hi = observed[, K_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = lo,
      ymax = hi, color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)
    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(K[g](r)))
  }

  ## plot for G
  else if (func == "G") {
    long <- as.data.frame(envelopes[[1]]$rrl_G)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    head(long)
    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_G)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }


    observed <- as.data.frame(observed_values$G)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, G_cor],
      lo = observed[, G_cor], hi = observed[, G_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo),
      ymax = (hi), color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(G[g](r)))
  }


  ### Plot for F
  else if (func == "F") {
    long <- as.data.frame(envelopes[[1]]$rrl_F)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_F)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }


    observed <- as.data.frame(observed_values$F)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, F_cor],
      lo = observed[, F_cor], hi = observed[, F_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo),
      ymax = (hi), color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(F[g](r)))
  }

  ### Plot for GXGH
  else if (func == "GXGH") {
    long <- as.data.frame(envelopes[[1]]$rrl_GXGH)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    head(long)
    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_GXGH)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }

    baseline <- (envelopes[[1]]$rrl_GXGH$mmean)

    observed <- as.data.frame(observed_values$GXGH)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, GXGH_cor],
      lo = observed[, GXGH_cor], hi = observed[, GXGH_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo),
      ymax = (hi), color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(G[gh](r)))
  }

  ### Plot for GXHG
  else if (func == "GXHG") {
    long <- as.data.frame(envelopes[[1]]$rrl_GXHG)
    long$type <- as.factor(1)
    long$type <- names(pattern.colors)[1]

    head(long)
    i <- 1
    while (i <= length(envelopes)) {
      temp <- as.data.frame(envelopes[[i]]$rrl_GXHG)
      temp$type <- as.factor(i)
      temp$type <- names(pattern.colors)[i]
      long <- rbind(long, temp)
      i <- i + 1
    }


    observed <- as.data.frame(observed_values$GXHG)
    observed <- data.frame(
      r = observed$r,
      mmean = observed[, GXHG_cor],
      lo = observed[, GXHG_cor], hi = observed[, GXHG_cor],
      type = "Observed"
    )
    long <- rbind(long, observed)

    gplot <- long %>% ggplot(aes(
      x = r, ymin = (lo),
      ymax = (hi), color = type, fill = type
    )) +
      geom_ribbon(alpha = 0.5) +
      geom_hline(yintercept = 0)

    gplot + theme(
      plot.title = element_text(hjust = 0.5, size = 60),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30),
      panel.border = element_rect(linetype = "solid", fill = NA),
      axis.title = element_text(size = 40),
      legend.key.size = unit(20, "pt"),
      legend.text = element_text(size = 20),
      # legend.title = element_text(size = 20, hjust = 0.5),
      legend.position = c(0.75, 0.80),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "white", color = "black"), ...
    ) +
      guides(color = "none") +
      scale_color_manual(values = pattern.colors) + scale_fill_manual(values = fill.colors) +
      xlab(paste("Radius (", unit, ")", sep = "")) + ylab(expression(G[hg](r)))
  }
}


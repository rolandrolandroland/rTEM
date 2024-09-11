#' Convert list of features to wide format
#' @param feats vector of feats
#' @description
#' Takes a list of lists and converts to single wide matrix
#' @export
feat_list_to_wide = function(feats, names, col_names, feat_names) {
  df = data.frame()
  for (i in 1:length(feats)) {
    for (j in 1:length(feats[[i]])) {
      for (k in 1:length(feats[[i]][[j]])) {
        temp = data.frame(feats[[i]][[j]][[k]])[,feat_names]
        temp = apply(temp, 2, max)
        temp[[names(col_names[1])]] = col_names[[1]][i]
        temp[[names(col_names[2])]] = col_names[[2]][j]
        temp[[names(col_names[3])]] = col_names[[3]][k]
        df = rbind(df, temp)
      }
    }
  }

  colnames(df) = c(feat_names, names(col_names))
  df
}

#' Combine ranks
#' @description
#' Takes a tall formatted matrix and combines the ranks of feature `feat_name` for functions `funcs`
#' @export
combine_rank = function(tall_data, feature_name, funcs = c("G", "K"), im = 1) {
  #dat = data.frame()
  dat = sapply(funcs, function(current_func) {
    print(current_func)
    # filter by feature name, the function, and the image number
    # the arrange by function value
    # add a "rank" column that is is the row number of the
    position = tall_data %>%
      filter(feat_name == feature_name & func == current_func & image_num == im)  %>%
      mutate(ind = row_number()) %>%
      arrange(as.numeric(value)) %>%
      mutate(rank = row_number())
    position$combo_num = as.numeric(position$combo_num)
    position = position %>% arrange((ind))
    position$rank
  })
  zweights = tall_data %>%
    filter(feat_name == feature_name & func == funcs[1] & image_num == im)  %>%
    mutate(ind = row_number()) %>%
    arrange(as.numeric(ind)) %>%
    select(zweight, combo_num)
  #score = dat[,1]
  score = rowSums(dat)
  #print(score)
  df = data.frame("score" = as.numeric(score),
                  combo_num = (zweights[,1]),
                  zweight = (zweights[,2]))
}

#' Get tall data from list of lists
#' @export
get_tall = function(iso_feats, vert_feats, feat_names,
                    tall_names = names, all_funcs = all_funcs
) {

  vert_c_names = list("image_num" = 1:length(vert_feats),
                      "combo_num" = 1:length(vert_feats[[1]]),
                      "func" = all_funcs)

  vert_wide = feat_list_to_wide(vert_feats, tall_names, vert_c_names, feat_names)

  iso_c_names = list("image_num" = 1:length(iso_feats),
                     "combo_num" = 1:length(iso_feats[[1]]),
                     "func" = all_funcs)
  iso_wide = feat_list_to_wide(iso_feats, tall_names, iso_c_names, feat_names)


  iso_tall = pivot_longer(iso_wide, feat_names,  names_to = "feat_name", values_to = "value")


  vert_tall =pivot_longer(vert_wide, feat_names,  names_to = "feat_name", values_to = "value")

  vert_tall$zweight = 0
  iso_tall$zweight = 1
  all_tall = rbind(vert_tall, iso_tall)
  all_tall$image_num = as.numeric(all_tall$image_num)
  all_tall$combo_num = as.numeric(all_tall$combo_num)
  all_tall$value = as.numeric(all_tall$value)
  all_tall
}

#' Get and save ordered tables for feature values of functions
#' @export
get_tables = function(all_tall, func_groups, feature_name, name_image, image_num, path) {
  ranks = lapply(1:length(func_groups), function(i) {
    ranks = combine_rank(all_tall, feature_name = feature_to_use, funcs = func_groups[[i]], im = image_num)
    ranks_sorted = ranks %>% arrange((score))
    cbind(size_fracs[as.numeric(ranks_sorted$combo_num),], ranks_sorted$zweight)[1:rows_to_print,]
  })
  ranks_mat = as.data.frame((ranks[[1]]))
  for (i in 2:length(ranks)) {
    ranks_mat = cbind(ranks_mat, ranks[[i]])
  }
  headers = sapply(func_groups, function(group) {
    if (length(group) == 1) {
      group
    }
    else {
      print(group)
      paste(group, collapse = " ")
    }
  })
  # Automatically create the header assignments, each with 4 columns
  header_assignments <- setNames(rep(5, length(headers)), headers)
  colnames(ranks_mat) = rep(c(1, 2, 3, 4, "z"), length(func_groups))
  caption_text <- paste0("<center>",
                         '<span style="font-size: 20px; color: black;">',
                         name_image, " Image ", image_num,
                         '</span>'
  )

  # Create the table with dynamic subheadings and visual gaps
  table_output =kable(ranks_mat, "html", col.names = colnames(ranks_mat), caption = caption_text) %>%
    add_header_above(header_assignments) %>%
    kable_styling(bootstrap_options = "striped", full_width = F) %>%
    column_spec(5, border_right = TRUE) %>%  # Adding a visual border for gaps between groups
    column_spec(10, border_right = TRUE) %>%
    column_spec(15, border_right = TRUE) %>%
    column_spec(20, border_right = TRUE) %>%
    column_spec(25, border_right = TRUE) %>%
    column_spec(30, border_right = TRUE) %>%
    kable_styling(latex_options="scale_down") %>%
    column_spec(1,width = "0.5in") %>%
    column_spec(4,width = "0.2in")
  type = ".html"
  file_name = paste(path, "combos_Image_", name_image, "_", image_num, type, sep = "")

  save_kable(table_output, file = file_name)


  type = ".png"
  png_name = paste(path, "combos_Image_", name_image, "_", image_num, type, sep = "")
  save_kable(table_output, png_name, vwidth = 2000, zoom = 2)
  table_output
}

#' Stack png files into a pdf
#' @export
stack_png_to_pdf = function(input_path, save_path, to_name, y_space = 50, resolution = 70) {
  png_files <- list.files(path = input_path, pattern = "\\.png$")
  images <- lapply(paste(input_path, png_files, sep = ""), readPNG)
  ln <- length(images)
  # Get total height and max width for the stacked image

  total_height <- sum(sapply(images, function(img) dim(img)[1])) + (y_space * (ln-1)) # Sum of heights
  max_width <- max(sapply(images, function(img) dim(img)[2]))     # Maximum width
  pdf(paste(save_path, to_name, ".pdf", sep = ""), width = max_width / resolution, height = total_height / resolution, paper = "special")
  plot(NA, xlim = c(0, max_width), ylim = c(0, total_height), type = "n", xaxt = "n", yaxt = "n", bty = "n", xaxs = "i", yaxs = "i")
  y_offset <- total_height
  for (i in 1:ln) {
    img_height <- dim(images[[i]])[1]
    img_width <- dim(images[[i]])[2]

    # Calculate the bottom left and top right corners for each image
    y_offset <- y_offset - img_height
    rasterImage(images[[i]], 0, y_offset, img_width, y_offset + img_height)

    y_offset = y_offset - y_space
  }

  # Close the PDF device
  dev.off()
}

#' Stack png files into a png file
#' @export
stack_png_to_png = function(input_path, save_path, to_name, y_space = 50, resolution = 70) {
  # Load PNG images
  png_files <- list.files(path = input_path, pattern = "\\.png$")
  images <- lapply(paste(input_path, png_files, sep = ""), readPNG)

  ln <- length(images)

  # Set space between images (in pixels)

  # Get total height and max width for the stacked image (including spaces between images)
  total_height <- sum(sapply(images, function(img) dim(img)[1])) + (y_space * (ln - 1))
  max_width <- max(sapply(images, function(img) dim(img)[2]))  # Maximum width

  # Open a PNG device with dimensions based on the total height and width
  # Adjust resolution with the `res` argument (e.g., 300 for high-resolution output)
  #png("stacked_tables.png", width = max_width, height = total_height, res = 72)
  png(paste(save_path, to_name, ".png", sep = ""), width = max_width , height = total_height, res = resolution)

  # Set up an empty plot with the right dimensions
  plot(NA, xlim = c(0, max_width), ylim = c(0, total_height), type = "n", xaxt = "n", yaxt = "n", bty = "n", xaxs = "i", yaxs = "i")

  # Initialize y_offset at the top (since plotting in R starts from bottom-left corner)
  y_offset <- total_height

  # Loop through the images and place them one after the other with space in between
  for (i in 1:ln) {
    img_height <- dim(images[[i]])[1]
    img_width <- dim(images[[i]])[2]

    # Calculate the bottom left and top right corners for each image, considering y_space
    y_offset <- y_offset - img_height
    rasterImage(images[[i]], 0, y_offset, img_width, y_offset + img_height)

    # Add space between images
    y_offset <- y_offset - y_space
  }

  # Close the PNG device
  dev.off()

}

#' assemble png files into a png square
#' @export
get_png_square = function(input_path, save_path, pat, to_name, y_space = 50, x_space = 50, resolution = 70) {
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




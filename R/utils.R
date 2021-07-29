####################################################
# utility functions for use in plotVisium() function
####################################################

.get_pal <- function(pal, val) {
    # (if missing) default to 'viridis' for continuous 
    # and 'libd_layer_colors' for discrete scales
    if (is.null(pal)) 
        pal <- ifelse(is.numeric(val), "viridis", "libd_layer_colors")
    # if length(palette) > 1, use palette as provided 
    # (either multiple colors for discrete labels, 
    # or length 2 for continuous gradient)
    if (length(pal) == 1) {
        pal <- switch(pal, 
            "libd_layer_colors" = c(
                "#F0027F", "#377EB8", "#4DAF4A", "#984EA3", 
                "#FFD700", "#FF7F00", "#1A1A1A", "#666666"), 
            "Okabe-Ito" = c(
                "#000000", "#E69F00", "#56B4E9", "#009E73", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), 
            # use 'scale_fill_viridis_c' instead of hex codes or color names
            "viridis" = pal, 
            # if providing a single color name (e.g. 'navy' or 'red'), 
            # combine with 'gray95' for continuous color scale
            c("gray95", pal)
        )
    }
    return(pal)
}

#' @importFrom SpatialExperiment imgData
.sub_imgData <- function(spe, sample_ids, image_ids) {
    .get_img_idx <- SpatialExperiment:::.get_img_idx
    if (is.null(image_ids)) {
        # default to first available image for each sample
        idx <- .get_img_idx(spe, TRUE, NULL)
    } else {
        if (length(image_ids) == 1) {
            idx <- .get_img_idx(spe, TRUE, image_ids)
        } else {
            stopifnot(length(image_ids) == length(sample_ids))
            idx <- mapply(s = sample_ids, i = image_ids,
                          function(s, i) .get_img_idx(spe, s, i))
        }
    }
    imgData(spe)[idx, ]
}

.y_reverse <- function(df, ix, y, img) {
    y_tmp <- df[ix, y]
    if (!is.null(img)) {
        y_tmp <- nrow(img) - y_tmp
    } else {
        y_tmp <- max(y_tmp) - y_tmp
    }
    df[ix, y] <- y_tmp
    return(df)
}

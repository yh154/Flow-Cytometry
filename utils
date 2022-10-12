fun_ExIPlot=function (object, features, type = "violin", idents = NULL, ncol = NULL, 
    sort = FALSE, assay = NULL, y.max = NULL, same.y.lims = FALSE, 
    adjust = 1, cols = NULL, pt.size = 0, group.by = NULL, split.by = NULL, 
    log = FALSE, slot = "data", stack = FALSE, combine = TRUE, 
    fill.by = NULL, flip = FALSE, raster = NULL, return_data=T) 
{
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    if (isTRUE(x = stack)) {
        if (!is.null(x = ncol)) {
            warning("'ncol' is ignored with 'stack' is TRUE", 
                call. = FALSE, immediate. = TRUE)
        }
        if (!is.null(x = y.max)) {
            warning("'y.max' is ignored when 'stack' is TRUE", 
                call. = FALSE, immediate. = TRUE)
        }
    }
    else {
        ncol <- ncol %||% ifelse(test = length(x = features) > 
            9, yes = 4, no = min(length(x = features), 3))
    }
    data <- FetchData(object = object, vars = features, slot = slot)
    pt.size <- pt.size %||% AutoPointSize(data = object)
    features <- colnames(x = data)
    if (is.null(x = idents)) {
        cells <- colnames(x = object)
    }
    else {
        cells <- names(x = Idents(object = object)[Idents(object = object) %in% 
            idents])
    }
    data <- data[cells, , drop = FALSE]
    idents <- if (is.null(x = group.by)) {
        Idents(object = object)[cells]
    }
    else {
        object[[group.by, drop = TRUE]][cells]
    }
    if (!is.factor(x = idents)) {
        idents <- factor(x = idents)
    }
    if (is.null(x = split.by)) {
        split <- NULL
    }
    else {
        split <- object[[split.by, drop = TRUE]][cells]
        if (!is.factor(x = split)) {
            split <- factor(x = split)
        }
        if (is.null(x = cols)) {
            cols <- hue_pal()(length(x = levels(x = idents)))
            cols <- Interleave(cols, InvertHex(hexadecimal = cols))
        }
        else if (length(x = cols) == 1 && cols == "interaction") {
            split <- interaction(idents, split)
            cols <- hue_pal()(length(x = levels(x = idents)))
        }
        else {
            cols <- Col2Hex(cols)
        }
        if (length(x = cols) < length(x = levels(x = split))) {
            cols <- Interleave(cols, InvertHex(hexadecimal = cols))
        }
        cols <- rep_len(x = cols, length.out = length(x = levels(x = split)))
        names(x = cols) <- levels(x = split)
        if ((length(x = cols) > 2) & (type == "splitViolin")) {
            warning("Split violin is only supported for <3 groups, using multi-violin.")
            type <- "violin"
        }
    }
    if (same.y.lims && is.null(x = y.max)) {
        y.max <- max(data)
    }
    if (isTRUE(x = stack)) {
        return(MultiExIPlot(type = type, data = data, idents = idents, 
            split = split, sort = sort, same.y.lims = same.y.lims, 
            adjust = adjust, cols = cols, pt.size = pt.size, 
            log = log, fill.by = fill.by, flip = flip))
    }
    plots <- lapply(X = features, FUN = function(x) {
        return(SingleExIPlot(type = type, data = data[, x, drop = FALSE], 
            idents = idents, split = split, sort = sort, y.max = y.max, 
            adjust = adjust, cols = cols, pt.size = pt.size, 
            log = log, raster = raster))
    })
    label.fxn <- switch(EXPR = type, violin = if (stack) {
        xlab
    } else {
        ylab
    }, splitViolin = if (stack) {
        xlab
    } else {
        ylab
    }, ridge = xlab, stop("Unknown ExIPlot type ", type, call. = FALSE))
    for (i in 1:length(x = plots)) {
        key <- paste0(unlist(x = strsplit(x = features[i], split = "_"))[1], 
            "_")
        obj <- names(x = which(x = Key(object = object) == key))
        if (length(x = obj) == 1) {
            if (inherits(x = object[[obj]], what = "DimReduc")) {
                plots[[i]] <- plots[[i]] + label.fxn(label = "Embeddings Value")
            }
            else if (inherits(x = object[[obj]], what = "Assay")) {
                next
            }
            else {
                warning("Unknown object type ", class(x = object), 
                  immediate. = TRUE, call. = FALSE)
                plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
            }
        }
        else if (!features[i] %in% rownames(x = object)) {
            plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
        }
    }
    if (combine) {
        plots <- patchwork::wrap_plots(plots, ncol = ncol)
        if (length(x = features) > 1) {
            plots <- plots & NoLegend()
        }
    }
    if(return_data){
    return(data)
    }else{return(plots)}
}

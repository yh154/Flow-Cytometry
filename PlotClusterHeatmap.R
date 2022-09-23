PlotClusterHeatmap = function (sce, features = rownames(sce), clusters = sce$cluster_id, 
    by_exprs_values = "exprs", fun = "median", scale = T, cluster_rows = T, 
    cluster_anno = F, draw_dend = T, draw_freqs = T, split_by = NULL, hm2 = NULL) 
{   require("ComplexHeatmap")
    u <- c("abundances", features)
    if (!is.null(hm2)) 
        stopifnot(hm2 %in% u)
    if (is.null(levels(clusters))) 
        clusters <- factor(clusters)
    sce$cluster_id <- clusters
    nk <- nlevels(sce$cluster_id)
    ms_by_k <- t(iMUBAC:::aggregateData(sce, by_exprs_values = by_exprs_values, 
        by = "cluster_id", fun = fun))
    d <- dist(ms_by_k[, features])
    if (cluster_rows) {
        row_clustering <- hclust(d, method = "complete")
    }
    else {
        row_clustering <- FALSE
    }
    if (cluster_anno) {
        anno <- levels(sce$cluster_id)
        if (nk > length(iMUBAC:::myCols)) {
            cols <- colorRampPalette(iMUBAC:::myCols)(nk)
        }
        else {
            cols <- iMUBAC:::myCols[seq_len(nk)]
        }
        cols <- setNames(cols, anno)
        cluster_anno <- ComplexHeatmap::Heatmap(matrix = anno, col = cols, name = "cluster_id", 
            rect_gp = grid::gpar(col = "white"), width = unit(0.4, 
                "cm"), cluster_rows = row_clustering, cluster_columns = T, 
            show_row_dend = draw_dend, row_dend_reorder = F)
    }
    many <- !is.null(split_by)
    cs <- seq_len(ncol(sce))
    if (many) 
    {groups <- split(cs, sce[[split_by]])}else{
          groups <- list(cs);print(length(groups))}
    if (scale) 
        assay(sce, by_exprs_values) <- iMUBAC:::scale_exprs(assay(sce, 
            by_exprs_values), 1)
    pals <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
    hm_cols <- colorRampPalette(pals)(100)
    hms <- lapply(seq_along(groups), function(i) {
        idx <- groups[[i]]
        cs_by_k <- split(idx, sce$cluster_id[idx])
        if (!many) {
            if (scale) {
                hm1_es <- t(iMUBAC:::aggregateData(sce, by_exprs_values = by_exprs_values, 
                  by = "cluster_id", fun = fun))
            }
            else {
                hm1_es <- ms_by_k
            }
        }
        else {
            hm1_es <- t(iMUBAC:::aggregateData(sce[, idx], by_exprs_values = by_exprs_values, 
                by = "cluster_id", fun = fun))
        }
        hm1 <- ComplexHeatmap::Heatmap(matrix = hm1_es[, features], col = hm_cols, 
            name = "expression", column_names_gp = grid::gpar(fontsize = 8), 
            rect_gp = grid::gpar(col = "white"), na_col = "lightgrey", 
            cluster_rows = row_clustering, cluster_columns = TRUE, 
            show_row_dend = draw_dend, column_title = names(groups)[i][many])
        freq_bars <- freq_anno <- NULL
        if (draw_freqs) {
            fq <- round(tabulate(sce$cluster_id[idx])/length(idx) * 
                100, 2)
            freq_bars <- ComplexHeatmap::rowAnnotation(`Freq [%]` = ComplexHeatmap::row_anno_barplot(fq, 
                axis = TRUE, border = FALSE, bar_with = 0.8, 
                gp = grid::gpar(fill = "grey50", col = "white")), 
                width = unit(2, "cm"))
            labs <- paste0(levels(sce$cluster_id), " (", fq, 
                "%)")
            freq_anno <- ComplexHeatmap::rowAnnotation(text = row_anno_text(labs), 
                width = max_text_width(labs))
        }
        p <- hm1 + freq_bars + freq_anno
        if (is(cluster_anno, "Heatmap")) 
            p <- cluster_anno + p
        if (!is.null(hm2)) {
            if (hm2 == "abundances") {
                cs <- table(sce$cluster_id[idx], sce$sample_id[idx])
                fq <- as.matrix(unclass(prop.table(cs, 2)))
                fq <- fq[, !is.na(colSums(fq)), drop = FALSE]
                p <- p + ComplexHeatmap::Heatmap(matrix = fq, name = "frequency", 
                  na_col = "lightgrey", rect_gp = grid::gpar(col = "white"), 
                  show_row_names = FALSE, column_names_gp = grid::gpar(fontsize = 8), 
                  cluster_rows = row_clustering, cluster_columns = FALSE)
            }
            else {
                for (ch in hm2) {
                  ms <- iMUBAC:::aggregateData(sce[ch, idx], by_exprs_values = by_exprs_values, 
                    by = c("cluster_id", "sample_id"), fun = fun)
                  ms <- do.call("rbind", ms)
                  rownames(ms) <- levels(sce$cluster_id)
                  p <- p + ComplexHeatmap::Heatmap(matrix = ms, col = hm_cols, 
                    na_col = "lightgrey", rect_gp = grid::gpar(col = "white"), 
                    show_heatmap_legend = FALSE, show_row_names = FALSE, 
                    cluster_rows = row_clustering, cluster_columns = FALSE, 
                    column_title = ch, column_names_gp = grid::gpar(fontsize = 8))
                }
            }
        }
        return(p)
    })
    hm_list <- NULL
    for (i in seq_along(hms)) hm_list <- hm_list + hms[[i]]
    draw(hm_list)
    invisible(hm_list)
}

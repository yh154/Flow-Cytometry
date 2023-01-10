BatchCorrection = function (sce, events = 1e+05, seed = 12345,batch="batch"){
  set.seed(seed)
  batches=table(sce[[batch]])
  slice_batch=events/length(batches)
  sce$event_id=1:ncol(sce)
  coldata=colData(sce) 
  coldata = split(coldata,coldata[[batch]])
  keep=list()
  for(i in seq_along(batches)){
     dt=coldata[[i]] %>% as.data.frame() %>% data.table::as.data.table() 
     fc=as.character(dt$sample_id) %>% unique() %>%length()
     slice_samp = floor(slice_batch/fc)
     keep[[i]]=dt[,sample(event_id,min(.N,slice_samp)), by=sample_id] %>%
       pull(V1)
  }
  sce <- sce[,sce$event_id%in%Reduce(c,keep)]
  res <- t(assay(sce, "exprs")) %>% scale() %>% 
    prcomp(retx = TRUE)
  data_mat = res$x
  colnames(data_mat)=rownames(sce)
  assay(sce, "normexprs") <- t(harmony::HarmonyMatrix(data_mat,         	      	
      meta_data=as.data.frame(colData(sce)),vars_use=batch, do_pca=F))
  return(sce)
}


RunUMAP = function (sce, 
                       by_exprs_values = "exprs", 
                       name = "UMAP", n_components = 2, 
                       n_neighbors = 100, min_dist = 0.5, 
                       scale = T, 
                       metric="euclidean",
                       n_threads = parallel::detectCores(logical = F),
                       seed = 12345)
{
  set.seed(seed)
  names <- paste0(name, 1:n_components)
  SingleCellExperiment::reducedDim(sce, type = name, withDimnames = T) <- 
    uwot::umap(t(assay(sce, by_exprs_values)), n_components = n_components, n_neighbors = n_neighbors,min_dist = min_dist, scale = scale, n_threads = n_threads,metric = metric,verbose = T) %>% magrittr::set_colnames(names)
  return(sce)
}

Clustering=function (sce, features = rownames(sce), by_exprs_values = "normexprs", 
    method = "SNNGraph", xdim = 20, ydim = 20, maxK = 40, save_graph=FALSE, output="",
    n_components = min(c(2,length(features))), n_neighbors = 10, min_dist = 0.1, resolution=1, seed = 12345) 
{
    set.seed(seed)
    if (identical(method, "FlowSOM")) {
        cat("1. FlowSOM clustering...\n")
        som <- FlowSOM::ReadInput(flowCore::flowFrame(t(assay(sce, 
            by_exprs_values))))
        som <- FlowSOM::BuildSOM(som, colsToUse = features, xdim = xdim, 
            ydim = ydim, silent = T)
        cat("2. ConsensusClusterPlus metaclustering...\n")
        pdf(NULL)
        mc <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(t(som$map$codes), 
            maxK = maxK, reps = 100, distance = "euclidean", 
            seed = seed, plot = NULL))
        dev.off()
        k <- xdim * ydim
        mcs <- seq_len(maxK)[-1]
        clstCodeDF <- data.frame(seq_len(k), purrr::map(mc[-1], 
            "consensusClass"))
        clstCodeDF <- dplyr::mutate_all(clstCodeDF, function(u) {
            factor(u, levels = sort(unique(u)))
        })
        colnames(clstCodeDF) <- c(sprintf("som%s", k), sprintf("meta%s", 
            mcs))
        clstIDs <- som$map$mapping[, 1]
        clstCodeDF <- data.frame(FlowSOM = clstCodeDF[[paste0("som", 
            k)]], Metaclustering = clstCodeDF[[paste0("meta", 
            maxK)]], stringsAsFactors = F)
        clstCodeDF <- clstCodeDF[which(clstCodeDF$FlowSOM %in% 
            unique(clstIDs)), ]
        sce$cluster_id <- factor(clstIDs, levels = clstCodeDF$FlowSOM, 
            labels = clstCodeDF$Metaclustering)
    }
    if (identical(method, "SNNGraph")) {
	      cat("SNN graph construction...\n")
        SingleCellExperiment::reducedDim(sce, type = "UMAP_SNN",withDimnames = T) <- t(assay(sce,"normexprs"))
        g <- scran::buildSNNGraph(sce, k=n_neighbors, use.dimred = "UMAP_SNN")
        if(save_graph){saveRDS(g,sprintf("%sSNNgraph.rds",output))}
        cat("Graph-based clustering...\n")
        for(i in seq_along(resolution)){
          g_comm <- igraph::cluster_louvain(g,resolution = resolution[i])
          sce[[paste0("cluster_id_res_",resolution[i])]]<- factor(g_comm$membership)
        }
    }
    return(sce)
}

DownSampleSCE = function(sce,
			 maxN,
			 group_by,
			 seed=800){
  set.seed(seed)
  sce$event_id=1:ncol(sce)
  dt <- colData(sce) %>% as.data.frame() %>% data.table::as.data.table()
  dt <- dt[, .SD[sample(.N, min(maxN,.N))], by = group_by]
  sce <- sce[, sce$event_id %in% dt$event_id]
  return(sce)
}

ClusterPropagation <- function(
    sce,
    sce_down,
    by_exprs_values="exprs",
    maxN=100,
    cluster_id="cluster_id",
    seed=12345
){
  require(caret)
  require(doParallel)
  batches=levels(sce$"batch")
  registerDoParallel()
  # The main function
  clusterPropagation_SingleBatch <- function(sce, sce_down, cluster_id, by_exprs_values,maxN,seed){
    
    set.seed(seed)
    sce_down <- DownSampleSCE(sce_down,maxN=maxN,group_by=cluster_id,seed=seed)
    
    ## Train a classifier
    cat("-Training a classifier...\n", sep="")
    trCtrls <- caret::trainControl(
      method="repeatedcv",
      number=10,
      repeats=5,
      sampling="up",
      verboseIter=F
    )
    mat <- t(assay(sce_down, by_exprs_values))
    tunegrid <- expand.grid(.mtry=5:15)
    #tunegrid <- expand.grid(.mtry=5:15, .numRandomCuts=1:2) ## for extraTrees
    caret.model <- caret::train(
      x=mat,
      y=droplevels(sce_down[[cluster_id]]),
      preProcess=c("center","scale"),
      method="rf",
      metric="Kappa",
      tuneGrid=tunegrid,
      trControl=trCtrls
    )
    ## Predict cell types
    cat("-Predicting clusters...\n", sep="")
    mat <- t(assay(sce, by_exprs_values))
    clusterIDs <- predict(caret.model, newdata=mat)
    clusterIDs <- as.numeric(as.character(clusterIDs)) ### remove factor levels
    return(clusterIDs)
  }
  
    clusterIDs=foreach(i=seq_along(batches)) %dopar% {
      clusterPropagation_SingleBatch(sce[,sce$batch==batches[i]],
                                     sce_down[,sce_down$batch==batches[i]],
                                     cluster_id = cluster_id,
                                     seed=seed,
                                     by_exprs_values = by_exprs_values,
                                     maxN = maxN)
      
    }

 sce$"predict_cluster_id" <- 0
 for(i in seq_along(batches)){
   sce[,sce$batch==batches[i]]$predict_cluster_id = clusterIDs[[i]]
 }
 sce$predict_cluster_id=factor(sce$predict_cluster_id)
 return(sce) 
}

PlotClusterHeatmap = function (sce, features = rownames(sce), clusters = sce$cluster_id, 
    by_exprs_values = "exprs", fun = "median", scale = T, cluster_rows = T, 
    cluster_anno = F, draw_dend = T, draw_freqs = T, split_by = NULL, hm2 = NULL,title=' ', main=' ') 
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
        cluster_anno <- ComplexHeatmap::Heatmap(matrix = anno, col = cols, name = title, 
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
            name = ifelse(fun=="mean", "avg.exp","median.exp"), column_names_gp = grid::gpar(fontsize = 8), 
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
    draw(hm_list, column_title = main)
    invisible(hm_list)
}

plotDensity = function(sce, exp="exprs",color="batch"){
  t(assay(sce, exp)) %>%
    as_tibble(rownames = "cell_id") %>%
    pivot_longer(!cell_id, names_to = "marker") %>%
    group_by(marker) %>%
    mutate(zscore = scale(value)) %>%
    ungroup() %>%
    # dplyr::filter(value > min_cutoff, value < max_cutoff) %>%
    dplyr::filter(zscore > -3, zscore < 3) %>%
    left_join(as_tibble(colData(sce), rownames = "cell_id"), by = "cell_id") %>%
    ggplot(aes_string(x = "value", color = color)) +
    geom_density() +
    facet_wrap(vars(marker), scales = "free") +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank()
    )
}


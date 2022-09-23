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

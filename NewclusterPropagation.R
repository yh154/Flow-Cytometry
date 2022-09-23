newClusterPropagation=function (sce, sce_down, by_exprs_values = "normexprs", maxN = 100, 
    numThreads = 4, seed = 12345) 
{
    clusterPropagation_SingleBatch <- function(sce, sce_down) {
        set.seed(seed)
        sce_down <- downsampleSCE(sce_down, maxN = maxN, group_by = c("batch", 
            "cluster_id"), indiv_by = NULL, seed = seed)
        cat("-Training a classifier...\n", sep = "")
        require(caret)
        require(doParallel)
        cl <- makePSOCKcluster(numThreads)
        registerDoParallel(cl)
        trCtrls <- caret::trainControl(method = "repeatedcv", 
            number = 10, repeats = 5, sampling = "up", verboseIter = F)
        mat <- t(assay(sce_down, by_exprs_values))
        tunegrid <- expand.grid(.mtry = 5:15)
        caret.model <- caret::train(x = mat, y = droplevels(sce_down$cluster_id), 
            preProcess = NULL, method = "rf", 
            metric = "Kappa", tuneGrid = tunegrid, trControl = trCtrls)
        stopCluster(cl)
        cat("-Predicting clusters...\n", sep = "")
        mat <- t(assay(sce,"exprs"))
        clusterIDs <- predict(caret.model, newdata = mat)
        clusterIDs <- as.numeric(as.character(clusterIDs))
        return(clusterIDs)
    }
    sce$cluster_id <- 0
    for (b in levels(sce$batch)) {
        cat("Batch: ", b, "\n", sep = "")
        clusterIDs <- clusterPropagation_SingleBatch(sce[, sce$batch == 
            b], sce_down[, sce_down$batch == b])
        sce[, sce$batch == b]$cluster_id <- clusterIDs
    }
    sce$cluster_id <- factor(sce$cluster_id, levels = unique(c(0, 
        as.numeric(as.character(levels(sce_down$cluster_id))))))
    return(sce)
}
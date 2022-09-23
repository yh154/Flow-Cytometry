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

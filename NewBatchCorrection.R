NewBatchCorrection = function (sce, events = 1e+05, seed = 12345,batch="batch"){
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

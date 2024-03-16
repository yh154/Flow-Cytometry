#!/usr/bin/env Rscript

'Flowcyto Batch Correction and Clustering.

Usage:
    flowcyto_bcc.R [options] (<csv>) [<output>]
Options:
    -h --help            Show this screen.
    --events=<bp>        Total events selected for anlaysis. 
                         (events/no-of-batch) of events are selected 
                         from each batch [default: 1e6].
    --batch=<bc>         Column name corresponds to "batch" [default: batch]. 
    -t --thread=<t>      Threads [default: 1].
    --group=<gp>         Column(s) to consider [default: batch,file_name]. Comma delimited. 
    --save_graph         Save SNN graph.
    --transfer           Perform arcsinh transfermation.
    --cofactor=<cf>      Transfermation cofactor [default: 5]
    --clustering         Perform clustering.
    --cluster_res=<res>  Clustering resolution [default: 1]. Comma delimited.    
    --propagate          Perform cluster propagation.
    --cluster_id         Which cluster id to propagate.

Arguments:
    csv Sample meta table with at least columns: full_name, batch, full_path, ...
    output Output directory. Default current working directory.

' -> doc

library(docopt)
opts <- docopt(doc)
print(opts)

if (!file.exists(opts$csv)) {
  stop(sprintf("%s does not exist.", opts$csv))
}

csv=opts$csv
events= as.integer(opts$events)
batch=opts$batch
output=opts$output
gp=as.character(unlist(strsplit(gsub("\\s+","",opts$group),",")))
if(is.null(output)){output=getwd()}
message(sprintf("Output file location:\n  %s\n",output))
thread=as.integer(opts$thread)
trans=opts$transfer
cofactor=as.integer(opts$cofactor)
res=as.numeric(as.character(unlist(strsplit(gsub("\\s+","",opts$cluster_res),","))))

if(trans){
    rds <- paste0(output, "/sce","_CF",cofactor,".rds")
    rds_down <- paste0(output, "/sce_down_",events,"_CF",cofactor,".rds")
 }else{
    rds <- paste0(output, "/sce",".rds")
    rds_down <- paste0(output, "/sce_down_",events,".rds")
}
if (file.exists(rds) | file.exists(rds_down)) {
  stop("sce object already exists.")
}

suppressPackageStartupMessages({
require(cowplot)
require(dplyr)
require(flowCore)
require(CATALYST)
require(sva)
require(ggplot2)
require(data.table)
require(caret)
require(doParallel)
})
source("utilities.R")

md=read.csv(csv) %>% data.table::setDT()
fs=read.flowSet(md[["full_path"]])
fsApply(fs,function(x){exprs(x) %>% range()})
CATALYST::guessPanel(fs[[1]])
message("\nCreate SCE ...")

## Create SingleCellExperiment object.
## md_cols and panel_cols may need edit accordingly
sce=CATALYST::prepData(
    fs,
    md=md,
    transform=trans,
    cofactor=cofactor,
    FACS = TRUE,
    md_cols = list(
        file = "file_name", id = "file_name",
        factors = gp), 
    panel_cols = list(channel = "fcs_colname",
        antigen = "antigen")
)

# logicle-transform
if(trans){
  message("\nPerform logicle-transform ...")
  ex <- fsApply(fsApply(fs, function(ff){
   lgcl <- estimateLogicle(ff, channels = colnames(ff))
   flowCore::transform(ff, lgcl)
  }),exprs)
  assay(sce, "counts", withDimnames=F) <- t(ex)
}

message(("\nDownSampling..."))
sce_down=DownSampleTotal(sce, total_events = events, seed=1234,batch = batch)
message("\nBatch correction...")
sce_down=BatchCorrection(sce_down,batch = batch)

##UMAP
message("\nRun UMAP ...")
sce_down <-RunUMAP(
  sce_down,
  by_exprs_values="counts",
  name="UMAPnorm",
  n_neighbors=25,
  min_dist=0.4,
  scale=F,
  n_threads=thread)
sce_down <-RunUMAP(
  sce_down,
  by_exprs_values="exprs",
  name="UMAP",
  n_neighbors=25,
  min_dist=0.4,
  scale=T,
  n_threads=thread)

message("\nPlotting ...")
p1=scater::plotReducedDim(sce_down, dimred = "UMAPnorm", colour_by = "batch",
                  point_size = 0.1, point_alpha = 1)+
  ggtitle("Before Batch Correction")+
  theme(plot.title = element_text(hjust = 0.5))
p2=scater::plotReducedDim(sce_down, dimred = "UMAP", colour_by = "batch",
                  point_size = 0.1, point_alpha = 1) +
  ggtitle("After Batch Correction")+
  theme(plot.title = element_text(hjust = 0.5))

message(paste0(output,"/UMAP_",events,".png"))
png(paste0(output,"/UMAP_",events,".png"),height = 400,width = 800)
cowplot::plot_grid(plotlist = list(p1,p2), ncol = 2)
dev.off()

##clustering
if(opts$clustering){
message("\nClustering ...")
sce_down<- Clustering(
    sce_down,
    method="SNNGraph",
    n_components=2,
    n_neighbors=30,
    save_graph=opts$save_graph,
    resolution=res
)
}

saveRDS(sce_down, rds_down)

## Cluster propagation
if(opts$propagate){
message("\nCluster propagation ...")
system.time(sce <- ClusterPropagation(
  sce,           
  sce_down,     
  by_exprs_values="exprs", 
  maxN=100,
  seed=12345,
  cluster_id = opts$cluster_id
))
}
saveRDS(sce, rds)

message("\nFinished Successfully!")

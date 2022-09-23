#!/usr/bin/env Rscript

'Flowcyto Batch Correction and Clustering.

Usage:
    flowcyto_bcc.R [--events output --thread] (--save_graph) <csv> <output>

Options:
    -h --help            Show this screen.
    --events=<bp>        Total events selected for batch correction.events/no-of-batch selected for each batch [default: 1e5].
    -t --thread=<t>      Threads [default: 1].
    --save_graph         Save SNN graph.

Arguments:
    csv Sample meta table with at least column: full_name, batch, full_path, ...
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
output=opts$output
if(is.null(output)){output=getwd()}
message(sprintf("Output file location:\n  %s\n",output))
thread=as.integer(opts$thread)

rds <- paste0(output, "/sce_", events, ".rds")
rds_down <- paste0(output, "/sce_down_",events,".rds")
if (file.exists(rds) | file.exists(rds_down)) {
  stop("sce object already exists.")
}

suppressPackageStartupMessages({
library(cowplot)
library(dplyr)
library(flowCore)
library(CATALYST)
library(ggplot2)
library(data.table)
})
source("BatchCorrection.R")
source("RunUMAP.R")
source("Clustering.R")

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
    transform=FALSE,
    FACS = TRUE,
    md_cols = list(
        file = "file_name", id = "file_name",
        factors = c("batch")), 
    panel_cols = list(channel = "fcs_colname",
        antigen = "antigen")
)

## logicle-transform
ex <- fsApply(fsApply(fs, function(ff){
    lgcl <- estimateLogicle(ff, channels = colnames(ff))
    flowCore::transform(ff, lgcl)
}),exprs)

assay(sce, "exprs", FALSE) <- t(ex)

sce_down=NewBatchCorrection(sce, events=100000)

##UMAP
message("\nRun UMAP ...")
sce_down <-NewRunUMAP(
  sce_down,
  by_exprs_values="normexprs",
  name="UMAPnorm",
  n_neighbors=30,
  min_dist=0.3,
  scale=F,
  metric = 'cosine',
  n_threads=thread)
sce_down <-NewRunUMAP(
  sce_down,
  by_exprs_values="exprs",
  name="UMAP",
  n_neighbors=30,
  min_dist=0.3,
  scale=T,
  n_threads=thread)

message("\nPlotting ...")
p1=scater::plotReducedDim(sce_down, dimred = "UMAP", colour_by = "batch",
                  point_size = 0.1, point_alpha = 0.2)+
  ggtitle("Before Batch Correction")+
  theme(plot.title = element_text(hjust = 0.5))
p2=scater::plotReducedDim(sce_down, dimred = "UMAPnorm", colour_by = "batch",
                  point_size = 0.1, point_alpha = 0.2) +
  ggtitle("After Batch Correction")+
  theme(plot.title = element_text(hjust = 0.5))

message(paste0(output,"/UMAP_",events,".png"))
png(paste0(output,"/UMAP_",events,".png"),height = 400,width = 800)
cowplot::plot_grid(plotlist = list(p1,p2), ncol = 2)
dev.off()

##clustering
message("\nClustering ...")
sce_down<- NewClustering(
    sce_down,
    method="SNNGraph",
    n_components=2,
    n_neighbors=30,
    save_graph=opts$save_graph,
    resolution=c(1,1.6)
)
saveRDS(sce_down, paste0(output,"/sce_down_",events,".rds"))
saveRDS(sce, paste0(output,"/sce_",events,".rds"))

#!/usr/bin/env Rscript

'Flowcyto Batch Correction and Clustering.

Usage:
  flowcyto_bcc.R csv [events output thread].
  flowcyto_bcc.R (-h | --help)

Options:
  -h --help     Show this screen.
  --csv  sample.csv contains columns: full_path, batch, full_path and et al.
  --events      Total events selected for batch correction. events/no-of-batch selected for each batch.
  --output    Output file. [default: ""].
  -t --thread   threads. [default: 1].

' -> doc

library(docopt)
opts <- docopt(doc)

csv=read.csv(opts$csv)
events= as.integer(opts$events)
output=opts$output
thread=opts$thread
if (!file.exists(csv)) {
  stop("csv does not exist: ", csv)
}

rds <- paste0(output, "/sce_",maxN,"".rds")
rds_down <- paste0(output, "/sce_down_",maxN,".rds")
if (file.exists(rds) | file.exists(rds_down)) {
  stop("sce object already exists.")
}

suppressPackageStartupMessages({
library(cowplot)
library(tidyr)
library(flowCore)
library(CATALYST)
library(iMUBAC)
library(data.table)
})

md=read.csv(csv) %>% data.table::setDT()
fs=read.flowSet(md[["full_path"]])
fsApply(fs,function(x){exprs(x) %>% range()})
CATALYST::guessPanel(fs[[1]]) 
message("\nCreate SCE ...")
sce=CATALYST::prepData(fs,md=md,
    transform=FALSE,
    FACS = TRUE,
    md_cols = list(
        file = "file_name", id = "file_name",
        factors = c("batch", "treatment","outcome","disease")),
    panel_cols = list(channel = "fcs_colname",
        antigen = "fcs_colname")
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
p1=iMUBAC::plotDR(sce_down, dimred = "UMAP", colour_by = "batch",
                  point_size = 0.1, point_alpha = 0.2)+
  ggtitle("Before Batch Correction")+
  theme(plot.title = element_text(hjust = 0.5))
p2=iMUBAC::plotDR(sce_down, dimred = "UMAPnorm", colour_by = "batch",
                  point_size = 0.1, point_alpha = 0.2) +
  ggtitle("After Batch Correction")+
  theme(plot.title = element_text(hjust = 0.5))
png(paste0(output,"UMAP_",events,".png",maxN),height = 400,width = 800)
cowplot::plot_grid(plotlist = list(p1,p2), ncol = 2)
dev.off()

##clustering
message("\nClustering ...")
sce_down<- NewClustering(
    sce_down,
    method="SNNGraph",
    n_components=2,
    n_neighbors=30,
    save_graph=TRUE,
    resolution=c(0.6,1,1.4,1.6)
)
saveRDS(sce_down, paste0(output,"sce_down_",events,".rds"))
saveRDS(sce, paste0(output,"sce_",events,".rds"))


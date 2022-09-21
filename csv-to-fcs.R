# Required packages
require('flowCore')
require('Biobase')
require('data.table')
CSV_Directory="/path/to/csv/"
FileNames <- list.files(path=CSV_Directory, pattern = ".csv")     
as.matrix(FileNames) 

out_dir="/output/dir/"

markers=read.table("markers_order.txt",sep="\t",header=F)$V1

DataList=list() 
for (File in FileNames) { 
  tempdata <- data.table::fread(paste0(PrimaryDirectory,File), check.names = FALSE, skip = 341)[,8:41]
  cns = sapply(colnames(tempdata),function(x){
    as.character(unlist(strsplit(x,"::")))[2]
  }) 
  cns = gsub("\\s+","",cns)
  colnames(tempdata)=cns 
  if(!all(sort(cns)==sort(markers))){
     cat(File)
     stop("markers different!")
  }
  tempdata=data.table::setcolorder(tempdata, markers)
  File <- gsub(".csv", "", File)
  DataList[[File]] <- tempdata
}

rm(tempdata)
AllSampleNames <- names(DataList)

for(i in c(1:length(AllSampleNames))){
  data_subset <- DataList[i]
  data_subset <- data.table::rbindlist(as.list(data_subset))
  dim(data_subset)

  metadata <- data.frame(name=dimnames(data_subset)[[2]],desc=paste('column',dimnames(data_subset)[[2]],'from dataset'))
  
  ## Create FCS file metadata - ranges, min, and max settings
  #metadata$range <- apply(apply(data_subset,2,range),2,diff)
  metadata$minRange <- apply(data_subset,2,min)
  metadata$maxRange <- apply(data_subset,2,max)
  
  data_subset.ff <- new("flowFrame",exprs=as.matrix(data_subset), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
  #head(data_subset.ff)
  write.FCS(data_subset.ff, paste0(out_dir,AllSampleNames[i],".fcs"))
}



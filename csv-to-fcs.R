# Required packages
require('flowCore')
require('Biobase')
require('data.table')
require('tibble')
require('dplyr')
csv_dir="/path/to/csv/"
FileNames <- list.files(path=csv_dir, pattern = ".csv")     
as.matrix(FileNames) 

out_dir="/output/dir/"

# markers of interests. 
markers=read.table("markers_order.txt",sep="\t",header=F)$V1 %>% toupper()

DataList=list() 
for (File in FileNames) { 
  message(File)
  tempdata <- data.table::fread(File, check.names = FALSE, skip = 391) ### check which rows to skip
  cns = sapply(colnames(tempdata),function(x){
    x=gsub("\\s+","",x)
    as.character(unlist(strsplit(x,"::")))[2]
  }) %>% toupper()
  
  colnames(tempdata)=cns 
  tempdata=tempdata[, ..markers]
  if(!all(sort(colnames(tempdata))==sort(markers))){
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
  if(nrow(DataList[[i]])>0){ 
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
}



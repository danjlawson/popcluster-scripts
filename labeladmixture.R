#!/usr/bin/env Rscript
## labeladmixture.R
## Author: Daniel Lawson dan.lawson@bristol.ac.uk (2023)
## Licence: GPLv3
## Purpose: Assign population labels to an admixture estimate
## Usage: call 'Rscript labeladmixture.R' or see help below.

args=commandArgs(trailingOnly=TRUE)
mycat=function(line,file=stdout()){
    cat(paste0(line,"\n"),file=file)
}

mincount=10
threshold=0.7
digits=2

if(length(args)>=2) {
    labelsf=args[1]
    admixturef=args[2]
    if(length(args)>2){
        threshold=as.numeric(args[3])
    }
}else{
    mycat("Showing help:",stderr())
    mycat("Usage: Rscript labeladmixture.R labels.csv admixture.csv (threshold) > output.csv",stderr())
    mycat("Purpose: annotate inferred admixture populations by the proportion of labels found in that population.")
    mycat("labels.csv : a file with individual labels as the first column and population labels as the rest of the line, comma separated, with a header.",stderr())
    mycat("admixture.csv : a csv file with a header, the first column is the individual label ids matching labels.txt and the rest are population admixture estimates",stderr())
    mycat("threshold : optional threshold for printing proportions, default 0.7")
    stop()
}

labels=read.csv(labelsf,header=TRUE,row.names=1)
admixture=as.matrix(read.csv(admixturef,header=TRUE,row.names=1))
labelmat=model.matrix(~0+LABEL,data=labels)
colnames(labelmat)<-gsub("LABEL","",colnames(labelmat))
labelmat = labelmat[,colSums(labelmat)>=mincount,drop=FALSE]


rawres=t(labelmat) %*% admixture
normres=rawres/rowSums(rawres)

## 
res=apply(normres,2,function(x){
    if(max(x)<threshold){
        tthreshold=max(x)
    }else  tthreshold=threshold
    tnames=rownames(normres)[x>=tthreshold]
    x=x[x>=tthreshold]
    x[order(x,decreasing=TRUE)]
},simplify=FALSE)

for(r in 1:length(res)){
    cat(paste0(names(res)[r],":\n"))
    for(i in 1:length(res[[r]])){
        if(i>1)cat(" + ")
        cat(paste(format(res[[r]][i],digits=digits),"*", names(res[[r]])[i]))
    }
    cat("\n")
}


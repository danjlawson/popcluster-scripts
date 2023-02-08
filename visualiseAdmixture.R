##

## I run with:
## Rscript  ~/code/popcluster-scripts/visualiseAdmixture.R MH7_comparison -n labels.txt MH_K_7_R_1_admixture.csv SNPs+AIMs_K7_fulldata_K_7_R_1_admixture.csv MH_top1000_pop7_K_7_R_1_admixture.csv SNPs_top1000_pop7_K_7_R_1_admixture.csv

mycat=function(line,file=stdout()){
    cat(paste0(line,"\n"),file=file)
}
args=commandArgs(trailingOnly=TRUE)
labelf=""
outputf="tmp"
mycat("Running...")
if(length(args)>=2) {
    outputf=args[1]
    for (i in 1:length(args)){
        if((args[i])=="-n"){
            labelf=args[i+1]
            args=args[-(i:(i+1))]
            break;
        }
    }
    files=args[-1]
}else{
    mycat("Showing help:",stderr())
    mycat("Usage: Rscript visualiseAdmixture.R output_file_root (-n labels.csv) (list of input files)",stderr())
    mycat("Purpose: visualise admixture file sets.")
    mycat("output_file_root: location for output files (.png used for graphs, reorded population matrices also created)")
    mycat("labels.csv : (optional) a file with population labels as the first column and population labels as the rest of the line, comma separated.",stderr())
    mycat("list of input files: _admixture.csv files as creted by popcluster2csv.R")
    stop()
}

mycat(paste0("Labels file: ",labelf))
mycat(paste0("Output root: ",outputf))
mycat(paste0("Input files: ",paste(files,collapse=",")))
mycat("Reading data...")

## files=c("MH-3-10_K_7_R_1_admixture.csv",
##         "SNPs+AIMs_K7_fulldata_K_7_R_1_admixture.csv",
##         "MH_top1000_pop7_K_7_R_1_admixture.csv",
##         "SNPs_top1000_pop7_K_7_R_1_admixture.csv")
qlist=lapply(files,function(x) as.data.frame(read.csv(x,row.names=1)))

popnames=colnames(qlist[[1]])
if(labelf!=""){
    mycat("Reading population labels...")
    popnames=read.csv(labelf,as.is=TRUE,header=FALSE)[,1]
}
mycat(paste0("labels: ",paste(popnames,collapse=",")))

##qlistA <-alignK(qlist,type="across")

align=function(q1,q2,fn=cor){
    if(dim(q1)[2]>dim(q2)[2])stop("Give the first admixture matrix as the smallest")
    q1=q1[,order(colSums(q1),decreasing=TRUE)]
    tcor=fn(q1,q2)
    indices=1:dim(q2)[2]
    res=numeric()
    for(i in 1:dim(tcor)[1]){
        tmp=order(tcor[i,],decreasing=TRUE)
        found=tmp[tmp%in%indices][1]
        res=c(res,found)
        indices=indices[!indices%in%res]
    }
    if(length(indices)>0) {
        cs=colSums(q2[,indices,drop=FALSE])
        indices=indices[order(cs,decreasing=TRUE)]
        res=c(res,indices)
    }
    ret=list(q1,q2[,res])
    ret
}
alignlist=function(qlist,fn=cor,ref="last"){
    res=align(qlist[[1]],qlist[[2]])
    if(length(qlist)>2) {
        for(i in 2:(length(qlist)-1)){
            iref=i
            if(ref=="first")iref=1
            res[[i+1]]=align(qlist[[iref]],qlist[[i+1]],fn)[[2]]
        }
    }
    names(res)=names(qlist)
    res
}
scorefn=function(qlist){
    weights=colSums(qlistA[[1]])
    rawscores=sapply(qlistA,function(x){
        (diag(cor(qlistA[[1]],x)))
    })
    weightedscores=apply(rawscores,2,function(x)x*weights/sum(weights))
    colSums(weightedscores)
}

mycat("Aligning...")
qlistA=alignlist(qlist,ref="first")
names(qlistA)<-gsub("_admixture.csv","",files)

porder=sapply(colnames(qlistA[[1]]),function(x){which(colnames(qlist[[1]])==x)})
popnames=popnames[porder]

scores=scorefn(qlistA)
scores2dp=round(scores,2)

tcor=cor(qlistA[[1]],qlistA[[2]])

mycat("Creating figure...")
cols=c("red","blue","green","orange","yellow","purple","cyan")
## cols <- c("#F2F3F4","#222222","#F3C300","#875692","#F38400","#A1CAF1","#BE0032","#C2B280","#848482","#008856","#E68FAC","#0067A5","#F99379","#604E97","#F6A600","#B3446C","#DCD300","#882D17","#8DB600", "#654522","#E25822","#2B3D26")
png(paste0(outputf,".png"),height=1200,width=1200)
k=length(qlistA)
par(mfrow=c((k+1),1),cex=1.5,mar=c(1,4,4,1))
plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab="",ylab="")
legend("center",fill=cols,legend=popnames,horiz=TRUE)
for(i in 1:k){
    barplot(t(as.matrix(qlistA[[i]])),horiz=FALSE,space=0,border=NA,
            col=cols,main=paste0(names(qlistA)[i],", score=",scores2dp[i]),
            xlab="",axisnames=FALSE,ylab="")
    axis(2)
}
dev.off()

mycat("Writing matrices...")
for(i in 1:k){
    write.csv(qlistA[[i]],file=paste0(outputf,"_",files[i]))
}


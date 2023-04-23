##
## I run with:
## Rscript  ~/code/popcluster-scripts/visualiseAdmixture.R MH7_comparison -n labels.txt -l individual_labels_asclusters.csv -c MH,SNPs+AIMs,MH\ top1000,SNPS\ top1000 MH_K_7_R_1_admixture.csv SNPs+AIMs_K7_fulldata_K_7_R_1_admixture.csv MH_top1000_pop7_K_7_R_1_admixture.csv SNPs_top1000_pop7_K_7_R_1_admixture.csv
## Rscript  ~/code/popcluster-scripts/visualiseAdmixture.R MH7_comparison -n labels.txt -c MH,SNPs+AIMs,MH\ top1000,SNPS\ top1000 MH_K_7_R_1_admixture.csv SNPs+AIMs_K7_fulldata_K_7_R_1_admixture.csv MH_top1000_pop7_K_7_R_1_admixture.csv SNPs_top1000_pop7_K_7_R_1_admixture.csv

## Editable quantities
marbarplot=c(0,4,6,1)
marlegend=c(8,4,0,1)
cex=0.8
legendheight=1.5
cexlegend=2
## Labels
cexlabels=1.75
srtlabels=60 # rotation
## title position
cextitle=2
linetitle=2
## y-axis
cexyaxis=2
lineyaxis=-1
cols=c("red","blue","green","orange","yellow","purple","cyan")
## cols <- c("#F2F3F4","#222222","#F3C300","#875692","#F38400","#A1CAF1","#BE0032","#C2B280","#848482","#008856","#E68FAC","#0067A5","#F99379","#604E97","#F6A600","#B3446C","#DCD300","#882D17","#8DB600", "#654522","#E25822","#2B3D26")

####################
mycat=function(line,file=stdout()){
    cat(paste0(line,"\n"),file=file)
}
args=commandArgs(trailingOnly=TRUE)
labelf=""
indlabelf=""
outputf="tmp"
sort=FALSE
captions=""
mycat("Running...")
if(length(args)>=2) {
    outputf=args[1]
    i=1
    while (i<length(args)){
        if((args[i])=="-n"){
            labelf=args[i+1]
            args=args[-(i:(i+1))]
        }else if((args[i])=="-l"){
            indlabelf=args[i+1]
            sort=TRUE
            args=args[-(i:(i+1))]
        }else if((args[i])=="-s"){
            sort=FALSE
            args=args[-i]
        }else if((args[i])=="-c"){
            captions=strsplit(args[i+1],",")[[1]]
            args=args[-(i:(i+1))]
        }else{
            i=i+1
        }
    }
    files=args[-1]
}else{
    mycat("Showing help:",stderr())
    mycat("Usage: Rscript visualiseAdmixture.R output_file_root (-n labels.csv) (-l indlabels.csv) (-s) (-c list of captions) (list of input files)",stderr())
    mycat("Purpose: visualise admixture file sets.")
    mycat("output_file_root: location for output files (.png used for graphs, reorded population matrices also created)")
    mycat("labels.csv : (optional) a file with population labels as the first column and latent population labels as the rest of the line, comma separated.",stderr())
    mycat("indlabels.csv : (optional) a file with individual labels as the first column and sample population labels as the second, comma separated.",stderr())
    mycat("-s: if providing indlabels, the output wil be sorted into clusters by default. Use this to disable sorting.")
    mycat("list of input files: _admixture.csv files as creted by popcluster2csv.R")
    stop()
}

####################
## Functions
##qlistA <-alignK(qlist,type="across")

align=function(q1,q2,fn=cor){
    ## Takes two admixture matrices q1,q2
    ## And a score function
    ## And aligns the columns of q2 to maximally match the columns of q1
    ## By sequentially aligning the columns in decreasing sum order
    ## Assuming dim(q2)[2]>=dim(q1)[2]
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
    ## Aligns a list of matrices qlist sequentially
    ## To either the first (ref=="first")
    ## Or the previous one (ref=="last"))
    ## Using the score function fn
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
scorefn=function(qlist,fn=cor){
    ## Scores a set of aligned admixture matrices by their alignment to the FIRST
    ## according to the score function fn
    ## By weighting the scores of each column by the total sum in that column, according to the first matrix
    weights=colSums(qlist[[1]])
    rawscores=sapply(qlist,function(x){
        (diag(fn(qlist[[1]],x)))
    })
    weightedscores=apply(rawscores,2,function(x)x*weights/sum(weights))
    colSums(weightedscores)
}
reorderlabel=function(A){
    ## Takes an admixture (sub)matrix
    ## And reorders the rows (individuals) according to the most frequent ancestry in that matrix
    ## Ties if they exist are broken by the second most frequent ancestry
    cm=colMeans(A)
    tA=A[,order(cm,decreasing=T)]
    row.names(A)[order(tA[,1],tA[,2],decreasing=TRUE)]
}
matrixorder=function(A,c){
    ## Takes an admixture matrix and a clustering of the rows in that matrix, expressed as a character string or other set of levels
    ## Returns the matrix ordered by those levels, and by admixture within those levels
    ## The levels appear in the same order of the FIRST occurance of the level
    lev=unique(c[,1])
    r=do.call("rbind",sapply(lev,function(x){
        cbind(reorderlabel(A[c[,1]==x,]),x)
    }))
    r
}
####################

## Start of processing
mycat(paste0("Population Labels file: ",labelf))
mycat(paste0("Individual Labels file: ",indlabelf))
mycat(paste0("Output root: ",outputf))
mycat(paste0("Input files: ",paste(files,collapse=",")))

mycat("Parsing captions...")
if(captions[1]==""){
    captions<-gsub("_admixture.csv","",files)
}
mycat(paste0("Using captions:", paste(captions,collapse=",")))

mycat("Reading data...")
qlist=lapply(files,function(x) as.data.frame(read.csv(x,row.names=1)))

popnames=colnames(qlist[[1]])
if(labelf!=""){
    mycat("Reading population labels...")
    popnames=read.csv(labelf,as.is=TRUE,header=FALSE)[,1]
}
mycat(paste0("labels: ",paste(popnames,collapse=",")))

if(indlabelf!=""){
    mycat("Reading individual labels...")
    ilabels=read.csv(indlabelf,row.names=1,sep=",")
    ilabels=ilabels[,dim(ilabels)[2],drop=FALSE]
    ulabels=unique(ilabels[,1])
    if(sort){
        mycat("Aligning using individual labels...")
        indorder=matrixorder(qlist[[1]],ilabels)
        changes=which(diff(as.numeric(as.factor(indorder[,2])))!=0)
        poplabels=ulabels
        mids=c(0,changes)+diff(c(0,changes,dim(indorder)[1]))/2
        qlist=lapply(qlist,function(x){x[indorder[,1],,drop=FALSE]})
    }else{
        mycat("NOT aligning using individual labels...")
        indorder=rownames(ilabels)
    }

}

mycat("Aligning...")
qlistA=alignlist(qlist,ref="first")

porder=sapply(colnames(qlistA[[1]]),function(x){which(colnames(qlist[[1]])==x)})
popnames=popnames[porder]

scores=scorefn(qlistA)
scores2dp=round(scores,2)

tcor=cor(qlistA[[1]],qlistA[[2]])

mycat("Creating figure...")
png(paste0(outputf,".png"),height=1500,width=1500)
k=length(qlistA)
layout(matrix(1:(k+1),ncol=1),
       heights=c(rep(1,k),legendheight))
##par(mar=c(0,0,0,0),cex=1.5)
##par(xpd=NA)
##plot(c(0,1),c(0,1),type="n",axes=FALSE,xlab="",ylab="")
##par(xpd=FALSE)
par(mar=marbarplot,cex=cex)
for(i in 1:k){
    barplot(t(as.matrix(qlistA[[i]])),horiz=FALSE,space=0,border=NA,axes=FALSE,
            col=cols,main="",
            xlab="",axisnames=FALSE,ylab="")
    axis(2,las=2,cex.axis=cexyaxis,line=lineyaxis)
    mtext(paste0(letters[1:k][i],") ",
                 captions[i],
                 ", score=",scores2dp[i]),3,
          adj=0,cex=cextitle,line=linetitle)
    if(sort){
        axis(1,at=mids,labels=rep("",length(mids)))
        abline(v=changes)
    }
}
par(mar=marlegend,cex=cex)
plot(c(1,dim(qlistA[[1]])[1]),c(0,1),type="n",axes=FALSE,xlab="",ylab="")
    if(sort){
        par(xpd=NA)
        text(mids,1,labels=poplabels,srt=srtlabels,adj=1,cex=cexlabels)
        par(xpd=FALSE)
    }
legend("bottomleft",fill=cols,legend=popnames,horiz=TRUE,cex=cexlegend)
dev.off()

mycat("Writing matrices...")
for(i in 1:k){
    write.csv(qlistA[[i]],file=paste0(outputf,"_",files[i]))
}


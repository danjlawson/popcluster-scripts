#!/usr/bin/env Rscript
## popcluster2csv.R
## 
## Author: Daniel Lawson dan.lawson@bristol.ac.uk (2023)
## Licence: GPLv3
## Purpose: Assign population labels to an admixture estimate
## Usage: call 'Rscript popcluster2csv.R' or see help below.
##
## Compute C_l from https://www.sciencedirect.com/science/article/pii/S1872497318303661
## (in Appendix A. Computational algorithms)
## Which is the total sum of KL divergence between all populations at that locus
## Compute this for every locus

## Usage cat popcluster_results_K_X_R_X | Rscript popcluster2csv.R <usage> > output.csv
## Where <usage> = Loci or Admixture (Admixture by default)
## e.g. 
## cat popcluster_results_K_X_R_X | Rscript popcluster2csv.R Loci > popcluster_results_K_X_R_X_loci.csv
## cat popcluster_results_K_X_R_X | Rscript popcluster2csv.R Admixture > popcluster_results_K_X_R_X_admixture.csv

args=commandArgs(trailingOnly=TRUE)
mycat=function(line,file=stdout()){
    cat(paste0(line,"\n"),file=file)
}
getscores<-function(mymatrix,eps){
  ## Return a matrix of information score comparisons, from a matrix of loci frequencies
  resforallele=sapply(1:npops,function(i){
    ii=i
    pi=mymatrix[,ii]
    pi[pi<eps]=eps
    logpi=log10(pi)
    piterm=pi*logpi
    sapply(1:npops,function(j){
      jj=j
      pj=mymatrix[,jj]
      pj[pj<eps]=eps
      logpj=log10(pj)
      pjterm=pi*logpj
      return(sum(abs(piterm-pjterm)))
    })
  })
  resforallele
}
eps=1e-4 # Minimum frequency for alleles to contribute to calculation 
ninds=-1
nloci=-1
npops=-1
locuson=0
alleleon=0
nalleles=0
screenprint=1000
extract="Admixture" ## options: "Admixture" or "Loci"
if(length(args)>0){
  if(args[1]=="Loci") {
    extract="Loci"
  }else if(args[1]=="Admixture"){
    extract="Admixture"
  }else{
    mycat("Showing help:",stderr())
    mycat("Usage: cat popclusterlogfile > Rscript popcluster2csv.R (mode) > output.csv",stderr())
    mycat("Where (mode) = Admixture : extract a csv of the estimated admixture",stderr())
    mycat("Where (mode) = Loci : extract a csv of the estimated Loci informativeness",stderr())
    stop()
  }
}
allelemat=matrix()
mode="NONE"
f <- file("stdin")
open(f, blocking=TRUE)

while(length(line <- readLines(f,n=1)) > 0) {
    if (startsWith( line,"--") ) {
        line <- readLines(f,n=1)
        if(startsWith(line,"PopCluster")){ ## Read the header
            mode="HEADER"
        }else if(startsWith(line,"Running parameters and options")){
            mycat(line,stderr())
            mode="PARAM"
        }else if(startsWith(line,"Summary of analyses")){
            mode="SUMMARY"
        }else if(startsWith(line,"Clustering analysis results")){
            mode="CLUSTERDETAILS"
        }else if (startsWith(line,"Inferred clusters")){
            mode="CLUSTERSUMMARY"
        }else if (startsWith(line,"Inferred Fst")){
            mode="FST"
        }else if (startsWith(line,"Admixture analysis")){
            mode="ADMIXTURE"
        }else if (startsWith(line,"Pairwise co-assignment probability")){
            mode="PAIRWISEPROB"
        }else if (startsWith(line,"Locus specific F-Statistics analysis")){
            mode="LOCUSFST"
        }else if (startsWith(line,"Allele frequency analysis")){
           ## if((ninds<0)|(nloci<0)|(npops<0))stop("Failed to parse header before reaching Allele Frequency Analysis!")
          if(extract=="Loci"){
            mycat(paste("!!! using parameters: ninds =",ninds,", nloci =",nloci,", npops =",npops),stderr())
            header=c("locus",paste0("pops:",apply(expand.grid(1:npops,1:npops),1,paste,collapse="_")),"mean")
            mycat(paste(header,collapse=","))
          }
          mode="ALLELEFREQ"
        }

        while(TRUE){ ## Skip until the content
            line <- readLines(f,n=1)
            if(length(line)==0)stop("Invalid File!")
            if (startsWith( line,"--") ) break;
        }
    }else{
        if(mode=="PARAM"){
            if(startsWith(line,"  Number of individuals")){
                ninds=as.integer(strsplit(line,"=")[[1]][2])
                mycat(paste("!!! found ninds=",ninds),stderr())
            }else if(startsWith(line,"  Number of loci")){
                nloci=as.integer(strsplit(line,"=")[[1]][2])
                mycat(paste("!!! found nloci=",nloci),stderr())
            }else if(startsWith(line,"  Number of assumed populations")){
                npops=as.integer(strsplit(line,"=")[[1]][2])
                mycat(paste("!!! found npops=",npops),stderr())
            }
        }else if((mode=="ADMIXTURE")&&(extract=="Admixture")) {
          if(startsWith(line,"Inferred")){
            mycat(paste("!!! processing admixture..."),stderr())
          }else if(startsWith(line,"Index")){
            mycat(paste0(c("Label",paste0("Pop",1:npops)),collapse = ","))
          }else{
            splitline=strsplit(line,"\\s+")[[1]]
            if(length(splitline)==7+npops){
              vals=splitline[c(4,7+(1:npops))]
              mycat(paste0(vals,collapse=","))
            }
          }
        }else if((mode=="ALLELEFREQ")&&(extract=="Loci")){
            if(startsWith(line,"Locus")){
                splitline=strsplit(line," ")[[1]]
                locuson=as.integer(gsub(":","",splitline[2]))
                nalleles=as.integer(splitline[3])
                alleleon=0
                if(( (locuson-1) %% screenprint)==0) {
                  mycat(paste("!!! processing locus",locuson,"with",nalleles,"alleles ..."),stderr())
                }
                mymatrix=read.table(f,nrows=nalleles,row.names=1)
                resforallele=getscores(mymatrix,eps)
#                print(mymatrix)
#                print(resforallele)
                res=as.numeric(resforallele)
                diag(resforallele)=NA
                res=c(locuson,res,mean(as.numeric(resforallele),na.rm = TRUE))
                mycat(paste(res,collapse=","))
            }
        }
    }
  # process line
}

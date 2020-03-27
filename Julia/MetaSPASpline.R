#!/usr/bin/env Rscript

# @Wallace Wang, wavefancy@gmail.com
#------------------------------

"
=============================================================
* Do meta-analysis based on SPA spline based method.

Usage:
MetaSPASpline.R
MetaSPASpline.R -h [--help]

Options:


Notes:
1. Read data from stdin and output results to stdout.
   Read in format as (at least 2 columns):
      ID spldata1 spldata2 ...
2. ref: 1. Dey R, et al. Robust meta‐analysis of biobank‐based genome‐wide 
   association studies with unbalanced binary phenotypes. Genet Epidemiol. 2019;43:462–476. 
3. SPAteset version >= 3.0.2.
=============================================================
" -> doc

suppressMessages(library(docopt))
suppressMessages(library(SPAtest))
suppressMessages(library(stringr))

opts <- docopt(doc)

##### Taking parameters ########


##### Start the main function ########
#data = read.table("./test.metaspline.txt", header = T,stringsAsFactors=F)
data = read.table(file("stdin"),header = T,stringsAsFactors=F,check.names=F)

fmt  = function(x){formatC(x, digits = 2, format = "f")}
fmtP = function(x){formatC(x, digits = 2, format = "e")}

n_meta = length(data[1,]) -1
out = c(colnames(data)[1],'PValue','NStudy','Is.Converge')
for(i in 1:nrow(data)){
    d = data[i,]
    n = 0
    meta_data = c()
    
    for(j in 1:n_meta){
        m_data = suppressWarnings(as.numeric(unlist(str_split(as.character(d[[1+j]]),','))))
        #print(m_data)
        if(sum(is.na(m_data))==0){
            meta_data = c(meta_data,m_data)
            n = n + 1
        }
    }
    
    # some gene may don't have a valid data.
    if(n ==0 ){
      out = rbind(out,c(d[[1]],'NA','NA','NA'))
    }else{
      # capture the output from a function to check converge.
      capture <- capture.output(p <- SPAmeta(spldata = matrix(meta_data,nrow=n, byrow=T)))
      #p = SPAmeta(spldata = matrix(meta_data,nrow=n, byrow=T))
      #print(length(capture))
      #print(c(d[[1]],fmtP(p),n))
      if (length(capture) == 0){
          out = rbind(out,c(d[[1]],fmtP(p),n,'YES'))
      }else{
          out = rbind(out,c(d[[1]],fmtP(p),n,'NO'))
      }  
    }
}
write.table(out,'',row.names = F, col.names = F, quote = F,sep = '\t')

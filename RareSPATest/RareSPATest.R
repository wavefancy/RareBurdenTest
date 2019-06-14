'
Rare variant burden/collapse test based on SPATest, which is a score test method corrected for
(qusi)-separation and inbalance case/control ratio. The overalll performance is very close to
Firth biased-corrected method, but with about 100x faster than Firth method.

Ref: A Fast and Accurate Algorithm to Test for Binary Phenotypes and Its Application to PheWAS

Notes:
  Read gene level score from stdin and output results to stdout.

Usage:
  RareSPATest.R -f file -p pheno -i id [-c covariates]

Options:
  -f file       Covariates file (including ID, Phenotype and covariates).
                  TSV with header, missing coded as NA.
  -p pheno      Phenotype name.
  -i id         Column name for individual id.
  -c covariates Specifiy covariates, eg.: cov1,cov2|cov1.
  ' -> doc

# load the docopt library and parse options.
suppressPackageStartupMessages(library(docopt))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(SPAtest))
opts <- docopt(doc)

# str(opts)

pedfile = ifelse(is.null(opts$f)==F, opts$f, NULL)
pheno   = ifelse(is.null(opts$p)==F, opts$p, NULL)
idname  = ifelse(is.null(opts$i)==F, opts$i, NULL)
covnames= ifelse(is.null(opts$c)==F, opts$c, NULL)

# https://stackoverflow.com/questions/15784373/process-substitution
# A method to help read process substitution.
OpenRead <- function(arg) {
   if (arg %in% c("-", "/dev/stdin")) {
      file("stdin", open = "r")
   } else if (grepl("^/dev/fd/", arg)) {
      fifo(arg, open = "r")
   } else {
      file(arg, open = "r")
   }
}
tfile = OpenRead(pedfile)
# **** IMPORTANT: read.table will auto ignore comment lines.
ped = read.table(tfile, header=T) # sep as white spaces: more spaces, tabs, newlines or carriage returns.
close(tfile)
# load data from stdin line by line,
# do the test gene by gene.
input <- file('stdin', 'r')
myhead <- stri_split_regex( readLines(input, n=1), "\\s+" )[[1]]

# subset to the individuals we have.
COVS = ifelse( is.null(covnames) == F, stri_replace_all_fixed(covnames,",","+"), "1")

# merge two data.frames.
ids = data.frame(ID=myhead)
# print(ids)
# print(colnames(ped))
datanull = merge(ids, ped, by.x="ID", by.y=idname, sort=T)
# print(datanull)
# print(pheno)
# print(COVS)
form = paste(pheno,COVS,sep="~")
# print(as.formula(form))
# SPANULL =  ScoreTest_wSaddleApprox_NULL_Model(as.formula(form), data=datanull)
# print(SPANULL)

datacol = 7
outtitle = T
COVS = stri_split_fixed(COVS,'+')[[1]]
for (line in readLines(input)){
    ss = stri_split_regex( line, "\\s+" )[[1]]
    if(outtitle==T){ #output title line
        out = c(myhead[1:datacol], c("NS", "FRAC_WITH_RARE", "PVALUE", "BETA", "SEBETA"))
        cat(out,sep="\t")
        cat("\n")
        outtitle=F
    }

    gene = data.frame(ID=myhead, SCORE = ss)
    data = merge(datanull, gene, by="ID", sort=T)

    genos = as.numeric(as.character(data[["SCORE"]]))
    # print(genos)

    out = c(ss[1:datacol], length(genos), length(which(genos>0))*1.0/length(genos))
    if(is.null(covnames))
    {
        fit = ScoreTest_SPA_wMeta(genos,data[pheno],minmac=1,Cutoff=0.1,output="metaZ")
    }else{
        # print(data)
        # print(as.character(data[["SCORE"]]))
        # print(as.numeric(as.character(data[["SCORE"]])))
        # print((data[[COVS]]))

        # fit = ScoreTest_SPA_wMeta(as.numeric(data[["SCORE"]]),data[pheno],data[COVS],SPANULL,minmac=1,Cutoff=0.1,output="metaZ")

        #first convert to character, and then convert to number, as 'score' is level in the data frames.
        # print(as.numeric(as.character(data[["SCORE"]])))
        fit = ScoreTest_SPA_wMeta(genos,data[pheno],data[COVS],
            minmac=1,Cutoff=0.1,output="metaZ",
            beta.out=T, beta.Cutoff = 1e-3
        ) # using fastSPA=0.1
    }

    out = c(out, fit$p.value, fit$beta, fit$SEbeta)
    cat(out,sep="\t")
    cat("\n")
}

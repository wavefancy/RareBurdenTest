# Restore the envirobment:
# conda env create -n burden -f conda_burden.yml
#
# OR by install by conda:
# conda create -n burden r-essentials r-base scipy pandas rpy2=2.9.4 docopt
# conda activate burden
# conda install -c conda-forge r-lmtest tzlocal

cat test.gene.txt | python3 RareNetworkBurdenRPython.py -c COV1,COV2 -f test.ped.txt -p PHENO -i ID -n test.network.txt

#GENES   BETAS   BETA_SE DF      CHISQ   PV
#Y,T     245.661,-343.925        614495,888558   2       6.46234 0.0395113
#M,T     -122.83,24.5661 307247,207146   2       6.46234 0.0395113

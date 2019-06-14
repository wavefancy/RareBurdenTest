wecho "
    cat test.gene.txt | Rscript RareSPATest.R -f test.ped.txt -p PHENO -i ID -c COV1
    # load the full gene data matrix into memory.
    # 3x faster than the read one by one version.
    !!cat test.gene.txt | Rscript RareSPATest.loadall.R -f test.ped.txt -p PHENO -i ID -c COV1
"


parallel -j1 -q wecho "
    cov=\`cat ../Julia/testgene/{2}.ped | head -n 1 | wcut -f3- | sed 's|{2}\t||' | sed 's|\t|,|g'\`
    &&
    gzcat ../Julia/testgene/{1}-{2}.score.gz
    | Rscript ./RareSPATest.R -f <(cat ../Julia/testgene/{2}.ped | sed 's|#||') -p {2} -i IND_ID
    #     -c 'pc1,pc2'
        -c \$cov
    > spa.{1}-{2}.txt
" ::: pcsk9 ldlr ::: MI EOMI

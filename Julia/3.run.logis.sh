
parallel -j1 -q wecho "
    cov=\`cat {2}.ped | head -n 1 | wcut -f3- | sed 's|{2}\t||' | sed 's|\t|,|g'\`
    &&
    gzcat testgene/{1}-{2}.score.gz
    #| grep -i pc
    | julia LogisticRegression.jl -f {2}.ped -p {2} -i IND_ID
    #     -c 'pc1,pc2'
        -c \$cov
    > testgene/logis.{1}-{2}.txt
" ::: pcsk9 ldlr ::: MI EOMI


parallel -j1 -q wecho "
    cov=\`cat {2}.ped | head -n 1 | wcut -f3- | sed 's|{2}\t||' | sed 's|\t|,|g'\`
    &&
    gzcat testgene/{1}-{2}.score.gz
    #| grep -i pc
    # There is one individual in LDLR gene has alt-homo which make the results different with epacts.
    | ppawk.py -u 'x=[min(1,float(f[i])) for i in range(7,len(f))];f[:7]+x'
    | julia LogisticRegression.R.jl -f {2}.ped -p {2} -i IND_ID
    #     -c 'pc1,pc2'
        -c \$cov
    > testgene/logis.{1}-{2}.txt
" ::: pcsk9 ldlr ::: MI EOMI

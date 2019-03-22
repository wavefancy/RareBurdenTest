
parallel -j1 -q wecho "
    java -jar ../target/RareBurden-1.0-SNAPSHOT-jar-with-dependencies.jar
        -g testgene/java.group.test.gene.txt -v testgene/{1}.gene.vcf.gz --test b.collapse
        -s <(cat testgene/{2}.ped | wcut -t 'IND_ID,{2}' | grep -vi NA | wcut -f1| tail -n +2)
        --max-maf 0.01
        | gzip > testgene/{1}-{2}.score.gz
" ::: pcsk9 ldlr ::: MI EOMI

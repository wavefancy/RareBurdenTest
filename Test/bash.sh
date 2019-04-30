wecho "
    java -jar ../target/RareBurden-1.0-SNAPSHOT-jar-with-dependencies.jar -g test.group.txt -v test.indels.missing.vcf.gz --test b.burden -s samples.txt > out.test.txt
"
wecho "
    java -jar ../target/RareBurden-1.0-SNAPSHOT-jar-with-dependencies.jar
        -g ./j.group.LOFHC.gene.txt
        -v ./test.sites.vcf.gz
        --test b.burden.weight
        -s <(cat ./ped.CAD.ped | wcut -t IND_ID | tail -n +2)
    | bgzip > out.ANGPTL8.txt.gz
"


dir="/medpop/esp2/wallace/projects/Amit/MIGEN_V11/LOFTEE_Burden"
wecho "
    scp $xfer:$dir/testgene/pcsk9.gene.vcf.gz testgene
    !scp $xfer:$dir/testgene/{EOMI-chr1.maf.0.01.epacts,MI-chr1.maf.0.01.epacts,MI-chr19.maf.0.01.epacts,EOMI-chr19.maf.0.01.epacts} testgene
    !scp $xfer:$dir/data/*.ped testgene
    !scp $xfer:$dir/testgene/ldlr.gene.vcf.gz testgene
    !scp $xfer:$dir/data/group.test.gene.txt testgene
    !cat testgene/group.test.gene.txt | sed -e 's|_|:|g' -e 's|/|:|g' > testgene/java.group.test.gene.txt

"

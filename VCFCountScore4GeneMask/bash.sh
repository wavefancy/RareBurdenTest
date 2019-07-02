# Test the capture of indels.
wecho "
    python3 VCFCountScore4GeneMask.py
        -g ./test.group.txt
        -v ./test.indels.missing.vcf.gz
        --weight maf
        --max-maf 0.01,1
        -s ./samples.txt
    | bgzip > out.test.txt.gz
"
# Test the -k keep option, supass the maf filterring.
wecho "
    python3 VCFCountScore4GeneMask.py
        -g ./test.group.txt
        -v ./test.indels.missing.vcf.gz
        --weight maf
        --max-maf 0.1
        -k <(echo -e '20:1234590:G:GTC\n20:14370:G:A')
        -s ./samples.txt
    | bgzip > out.keep.txt.gz
"
wecho "
    python3 VCFCountScore4GeneMask.py
        -g ./j.group.LOFHC.gene.txt
        -v ./hg19.gene.ANGPTL8.vcf.gz
        --weight file
        --max-maf 0.01,1
        -s <(bgzip -dc ./ped.CAD.ped.gz | wcut -t IND_ID | tail -n +2)
    | bgzip > out.ANGPTL8.txt.gz
"

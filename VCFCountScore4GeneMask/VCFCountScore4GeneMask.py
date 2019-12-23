#!/usr/bin/env python3

"""

    Aggregate/weight the number of rare variant for each individual in each gene.
        - A single pass can classify variant into multiple maf or mac bins.
        - Facilitate downstream burden based test for checking multiple conditions.
        *** Do not support multi-allelic sites, please norm the VCF.

    @Author: wavefancy@gmail.com

    Usage:
        VCFCountScore4GeneMask.py -g file -v bgzfile [--weight text] [-s file] [--max-maf float] [-n int] [-k file] [--max-mac int] [--maf-bin floats] [--mac-bin ints]
        VCFCountScore4GeneMask.py -h | --help | -v | --version | -f | --format

    Notes:
        1. Output to stdout, each line for each gene mask.
        2. *** Do not support multi-allelic sites in the current version.

    Options:
        -g file          Group file, each line is a group of variants.
        -v bgzfile       Bgziped and tabix indexed vcf file.
        --weight txt     The way for aggregate variant. Default: file.
                             file: sum(alt_count*variant_weight), variant_weight from group file.
                                   Set as 1 if group file don't have weight.
                             MAF:  Set the weight as 1/(MAF*(1-MAF))^0.5. Madsen and Browning (2009).
        -s file          Sample file, only count the score for those individuals decleared in this file.
        -n int           Set the number of threads, Default 2, no impove as more threads.
        --max-maf float  MAF filtering for alt allele (alt allele frequency), eg. 0.01.
        --max-mac int    MAC filtering for alt allele (alt allele count), eg. 5.
        --maf-bin floats MAF cut-off for alt allele (alt allele frequency), eg. 0.01|0.01,0.05.
                            For a single pass, the variants will be checked if meet any 'maf-bin' cut-off,
                            IF so, will contribute to the score for those bins it met the condition.
        --mac-bin ints   MAC cut-off for alt allele (alt allele count), eg. 5|1,3.
                            For a single pass, the variants will be checked if meet any 'mac-bin' cut-off,
                            IF so, will contribute to the score for those bins it met the condition.
        -k file          Always keep the variants listed in this file, surpass MAF and MAC filtering.
                            Put id line by line, CHR:POS:REF:ALT format.
        -h --help        Show this screen.
        --version        Show version.
        -f --format      Show format example.
"""
import sys
from docopt import docopt
from signal import signal, SIGPIPE, SIG_DFL
from cyvcf2 import VCF
from collections import OrderedDict
import numpy as np
import numexpr as ne
signal(SIGPIPE, SIG_DFL)

def ShowFormat():
    '''Input File format example:'''
    print('''
    ''');

if __name__ == '__main__':
    args = docopt(__doc__, version='1.1')
    # version 1.1
    # Check format at each line. FORMAT may different line by line.
    #print(args)

    if(args['--format']):
        ShowFormat()
        sys.exit(-1)

    group_file = args['-g']
    vcf_file = args['-v']
    weight_way = args['--weight'] if args['--weight'] else 'file'
    weight_way = weight_way.lower()
    if weight_way not in set(['file','maf']):
        sys.stderr.write('Please set a proper value for --weight')
    # set the number of threads
    N_threads = 2
    if args['-n']:
        N_threads = int(args['-n'])
        ne.set_num_threads(N_threads)

    samples = []
    if args['-s']:
        with open(args['-s'], 'r') as sfile:
            for line in sfile:
                line = line.strip()
                if line:
                    ss = line.split()
                    [samples.append(x) for x in ss]
    # convert list to set for fast checking
    # samples = set(samples)

    # The variant list we always want to keep, no matter what's the MAF and MAC filter is.
    KeepVIDs = set()
    if args['-k']:
        with open(args['-k'], 'r') as content_file:
            content = content_file.read().strip().split()
            KeepVIDs = set(content)
    #print(KeepVIDs)

    # Will check if a varinat met any maf or mac bin.
    mafs = [float(x) for x in args['--maf-bin'].split(',')] if args['--maf-bin'] else []
    macs = [int(x)   for x in args['--mac-bin'].split(',')] if args['--mac-bin'] else []
    MAX_MAF = float(args['--max-maf']) if args['--max-maf'] else 1
    MAX_MAC = int(args['--max-mac'])   if args['--max-mac'] else 9000000000 # I don't think we may have sample size more than this.
    MAX_MAF = min(MAX_MAF, max(mafs))  if mafs else MAX_MAF
    MAX_MAC = min(MAX_MAC, max(macs))  if macs else MAX_MAC
    
    # gts012: bool
    #    if True, then gt_types will be 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN. If False, 3, 2 are flipped.
    # strict_gt: True: half missing set as UNKNOWN. Otherwise half missing set at HET.
    # API: https://brentp.github.io/cyvcf2/docstrings.html#api
    # lazy or not has no big change on performance.
    if samples: #automatic substract the subset of samples we want.
        vcf = VCF(vcf_file, gts012=True, strict_gt=True, samples=samples, lazy=True, threads = N_threads)
    else:
        vcf = VCF(vcf_file, gts012=True, strict_gt=True, lazy=True, threads = N_threads)

    # print(vcf.samples)
    out = "#CHROM  BEGIN   END     MARKER_ID       NUM_ALL_VARS    NUM_PASS_VARS   NUM_SING_VARS MAF/MAC_CUT".split()
    sys.stdout.write('%s\n'%('\t'.join(out + vcf.samples)))
    # go through the group file and count the scores.
    with open(group_file, 'r') as gfile:

        for line in gfile: # Process by group file, line by line, each line is a gene.
            line = line.strip()
            if line: # for each gene
                ss = line.split()
                # record weight map
                NUM_ALL_VARS = 0  # the number of sites defined in group file, but also found in vcf file.
                record_weight = {}
                positions = []
                for snp in ss[1:]:
                    t = snp.split(':')
                    key = ':'.join(t[:4])
                    w = float(t[4]) if len(t) == 5 else 1.0
                    record_weight[key] = w
                    positions.append(int(t[1]))

                # results for different MAFs and MACs.
                MAF_RESULTS = OrderedDict()        # Individual scores.
                MAF_RESULTS_PASS = OrderedDict()   # Number of passed variants.
                MAF_RESULTS_SING = OrderedDict()   # Number of singleton.
                MAC_RESULTS = OrderedDict()
                MAC_RESULTS_PASS = OrderedDict()
                MAC_RESULTS_SING = OrderedDict()
                if mafs:
                    for x in mafs:
                        MAF_RESULTS[str(x)] = np.zeros(len(vcf.samples))
                        MAF_RESULTS_PASS[str(x)] = 0
                        MAF_RESULTS_SING[str(x)] = 0
                if macs:
                    for x in macs:
                        MAC_RESULTS[str(x)] = np.zeros(len(vcf.samples))
                        MAC_RESULTS_PASS[str(x)] = 0
                        MAC_RESULTS_SING[str(x)] = 0

                #META INFO for a chunk of VCF.
                chr  = ss[1].split(':')[0]
                minp, maxp = (min(positions), max(positions))
                out = [chr, str(minp), str(maxp), ss[0]]
                # print(record_weight)
                # iterate through variant
                query = '%s:%d-%d'%(chr, minp, maxp)
                for variant in vcf(query): # query in a gene region.
                    # e.g. REF='A', ALT=['C', 'T']
                    # convert position from 0 based to 1 based.
                    # print(variant.start+1)
                    # print(variant.ALT)
                    ALT_ALLELE = variant.ALT[0] if variant.ALT else '.'
                    id = '%s:%s:%s:%s'%(variant.CHROM, variant.start+1, variant.REF, ALT_ALLELE)
                    # print(id)
                    if id in record_weight:
                        # print(id)
                        NUM_ALL_VARS += 1
                        aaf = variant.aaf # alt allele frequency across samples in this VCF.
                        # v_maf = min(aaf, 1-aaf)
                        v_maf = aaf
                        alt_count = variant.num_het + variant.num_hom_alt*2
                        v_mac = alt_count

                        if v_maf == 0: continue # I don't think this kind of sites make sense for testing.
                        weight = np.power(1/(aaf * (1-aaf)),0.5) if weight_way == 'maf' else record_weight[id]

                        # Fake the maf(mac) as 0, to pass the maf filtering, in order to always keep this variant.
                        if KeepVIDs and (id in KeepVIDs):
                            v_maf = 0.0
                            v_mac = 0
                        if v_maf > MAX_MAF or v_mac > MAX_MAC: #max-maf or max-mac filtering.
                            continue

                        # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
                        genos = variant.gt_types
                        # alt_count = variant.num_het + variant.num_hom_alt*2
                        # impute missing as HOM_REF
                        genos[genos == 3] = 0
                        # convert geno count to weight
                        genos = ne.evaluate('genos * weight')
                        #genos[genos >0 ] *= weight

                        # filter by maf and do computation:
                        if mafs:
                            # aggregate results
                            for maf in mafs:
                                if v_maf <= maf:
                                    aid = str(maf)
                                    MAF_RESULTS_PASS[aid] = MAF_RESULTS_PASS[aid] + 1
                                    t = MAF_RESULTS[aid]
                                    MAF_RESULTS[aid] = ne.evaluate('t + genos')

                                    if alt_count == 1:
                                        MAF_RESULTS_SING[aid] = MAF_RESULTS_SING[aid] + 1

                        # filter by mac and do computation:
                        if macs:
                            # aggregate results
                            for mac in macs:
                                if v_mac <= mac:
                                    aid = str(mac)
                                    MAC_RESULTS_PASS[aid] = MAC_RESULTS_PASS[aid] + 1
                                    t = MAC_RESULTS[aid]
                                    MAC_RESULTS[aid] = ne.evaluate('t + genos')

                                    if alt_count == 1:
                                        MAC_RESULTS_SING[aid] = MAC_RESULTS_SING[aid] + 1

                # print(MAF_RESULTS)
                #output results
                for maf in mafs:
                    aid = str(maf)
                    myout = out + [str(NUM_ALL_VARS), str(MAF_RESULTS_PASS[aid]), str(MAF_RESULTS_SING[aid]), aid] + ['%g'%(x) for x in MAF_RESULTS[aid]]
                    sys.stdout.write('%s\n'%('\t'.join(myout)))
                for mac in macs:
                    aid = str(mac)
                    myout = out + [str(NUM_ALL_VARS), str(MAC_RESULTS_PASS[aid]), str(MAC_RESULTS_SING[aid]), aid] + ['%g'%(x) for x in MAC_RESULTS[aid]]
                    sys.stdout.write('%s\n'%('\t'.join(myout)))

sys.stdout.flush()
sys.stdout.close()
sys.stderr.flush()
sys.stderr.close()

#!/usr/bin/env python3

"""

    Network guided burden test, using loglikelihood ratio test (LRT) for deriving P value.
    Chi-square of LRT were also output for meta-analysis, sum of Chisq values.

    Model1: phenotype ~ covariates + 1
    Model2: phenotype ~ g1_score + g2_score ... + covariates +1
    LRT(Model1, Model2)

    @Author: wavefancy@gmail.com

    Usage:
        RareNetworkBurdenRPython.py [-c txt] -f file -p pheno -i txt -n file -g int
        RareNetworkBurdenRPython.py -h | --help | -v | --version | -f | --format

    Notes:
        1. Read gene scores from stdin, and output results to stdout,
        2. See example by -f.

    Options:
        -c txt         Covariate names, eg. cov1|cov1,cov2. No covariate if this is off.
        -f file        Covariate and phenotype file.
        -p pheno       Column name for phenotype name.
        -i txt         Column name for individual names.
        -g int         Column index for gene name, starts from 1.
        -n file        Network defining file, group genes, each group per line in file, white spaced.
        -h --help      Show this screen.
        -v --version   Show version.
        --format       Show input/output file format example.
    
    Dependencies:
        R: lmtest
        Python: rpy2, scipy, tzlocal.
        
"""
import sys
from docopt import docopt
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) #prevent IOError: [Errno 32] Broken pipe. If pipe closed by 'head'.

def ShowFormat():
    '''Input File format example:'''
    print('''
    ''')

if __name__ == '__main__':
    args = docopt(__doc__, version='1.0')
    #print(args)

    if(args['--format']):
        ShowFormat()
        sys.exit(-1)

    import pandas as pd
    COVARS   = args['-c'].split(',') if args['-c'] else []
    PHENO_N  = args['-p']
    ID_N     = args['-i']
    # pandas will auto impute the data type.
    df = pd.read_csv(args['-f'],sep=r'\s+')
    # check all covariates are in loaded phenotype and covariates files.
    m_cov = [c for c in COVARS + [PHENO_N, ID_N] if c not in df.columns]
    if m_cov:
        sys.stderr.write('ERROR! not found column name in cov file: %s\n'%('\t'.join(m_cov)))
        sys.exit(-1)
    COVARS += ['1']
    
    # check phenotype to determine use normal regression or logistic regress.
    LOGIS_REGRESS = True if len(set(df[PHENO_N])) == 2 else False
    link = 'binomial' if LOGIS_REGRESS else 'gaussian'
    if link == 'binomial':
        sys.stderr.write('The phenotype only has 2 distinct values, using logistic regression\n')
    else:
        sys.stderr.write('The phenotype has more than 2 distinct values, using linear regression\n')
    
    # load gene score file, and map file from gene_name -> line values.
    title = ''
    score_map = {} # genename -> line values.
    col_genename = 3
    col_genename = int(args['-g']) -1
    id_index = [] # value index for each individual in cov file.
    for line in sys.stdin:
        line = line.strip()
        if line:
            ss = line.split()
            if not title:
                title = ss
                # get index map for individuals.
                title_map = {}
                for x,y in enumerate(title):
                    if y in title_map:
                        sys.stderr.write('ERROR: duplicate title entries in score stream: %s\n'%(y))
                        sys.exit(-1)
                    title_map[y] = x
                
                # subset the df to the common ids.
                id_subset = [(str(x) in title_map) for x in df[ID_N]] # true/false array.
                # print(id_subset)
                df = df.loc[id_subset,] # update the dataframe to the subset.
                # sample index in the score file.
                id_index = [title_map[str(x)] for x in df[ID_N]]
            else:
                if ss[col_genename] in score_map:
                    sys.stderr.write('WARN: duplicated gene name in score file: %s\n'%(ss[3]))
                else:
                    scores = [float(ss[x]) for x in id_index]
                    if sum(scores) > 0: # do not need to test if no rare variants.
                        score_map[ss[col_genename]] = scores # score for each individual, and align to the order of cov.
  
    # load gene network structure, and do association test.

    # https://stackoverflow.com/questions/15777951/how-to-suppress-pandas-future-warning
    import warnings
    warnings.simplefilter(action='ignore', category=FutureWarning)

    from rpy2.robjects import r, pandas2ri
    # Need explicitly activate this, otherwise will fail.
    # https://pandas.pydata.org/pandas-docs/version/0.22/r_interface.html
    pandas2ri.activate()
    from rpy2.robjects import FloatVector,IntVector
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri

    from scipy.stats import chi2
    import numpy as np

    #import rpy2
    lrt = importr("lmtest")
    sts = importr("stats")
    rbase = importr("base")
    # print(df)
    r_null_df = pandas2ri.py2ri(df)
    # print(PHENO_N + '~' + '+'.join(COVARS))
    null_model = sts.glm(PHENO_N + '~' + '+'.join(COVARS), family = link, data=r_null_df)

    # plugin gene score, and do network guided burden test.
    sys.stdout.write('GENES\tBETAS\tBETA_SE\tDF\tCHISQ\tPV\n')
    with open(args['-n'], 'r') as nfile:
        for line in nfile:
            line = line.strip()
            if line:
                gnames = line.split()

                new_df = df.copy(deep=True)
                NUM_GENE = 0
                for g in gnames:
                    if g in score_map: # some genes may don't have scores.
                        new_df[g] = score_map[g]
                        NUM_GENE += 1
                if NUM_GENE < 2: # No enough data test for this network.
                    sys.stdout.write('%s\t%s\t%s\t%s\n'%(','.join(gnames), 'NA', 
                    'NA',
                    '\t'.join(['NA','NA','NA'])))
                    continue

                # print(new_df)
                # testing the alternative model and compare the LRT.
                alt_model = sts.glm(PHENO_N + '~' + '+'.join(gnames + COVARS), family = link, data=new_df)
                # print(type(alt_model)) , rpy2.robjects.vectors.ListVector
                # coefficient for the gene score.
                # print(type(sts.coefficients(alt_model))), rpy2.robjects.vectors.FloatVector
                # beta = sts.coefficients(alt_model)[1:1+len(gnames)]
                sout = rbase.summary(alt_model) 
                b = np.array(sout[sout.names.index('coefficients')]) # rpy2.robjects.vectors.Matrix, matrix to array.
                # print(type(b))
                beta   = b[1:1+len(gnames),0]
                betase = b[1:1+len(gnames),1]

                lrt_res = lrt.lrtest(null_model, alt_model)
                # print(type(lrt_res)) # rpy2.robjects.vectors.DataFrame
                # print(pandas2ri.ri2py(lrt_res))
                degree = lrt_res[lrt_res.names.index('Df')][1]
                cs = lrt_res[lrt_res.names.index('Chisq')][1]
                mychi2 = chi2(degree)
                pv = mychi2.sf(cs)

                sys.stdout.write('%s\t%s\t%s\t%s\n'%(','.join(gnames), ','.join(['%g'%(x) for x in beta]), 
                    ','.join(['%g'%(x) for x in betase]),
                    '\t'.join(['%g'%(x) for x in (degree, cs, pv)])))
                
                new_df = None
                

sys.stdout.flush()
sys.stdout.close()
sys.stderr.flush()
sys.stderr.close()

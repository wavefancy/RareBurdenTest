# Do logistic regression
# binary_phenotype ~ burden_score + covariates.
# Using R do the logistic regression, in order to match the results of epacts.
# SPA Ref: A Fast and Accurate Algorithm to Test for Binary Phenotypes and Its Application to PheWAS

doc = """Regression for rare variant burden test.

phenotype(binary/quantitative) ~ burden_score + covariates.

*** Call R glm to do [logistic]regression, to match the results of epacts,
*** Slow compared to Julia version of GLM
*** Auto detect binary or quantitative regression based the the number of uniq number in phenotype.
*** With function to use SPAtest for logistic regression,
    Firth logistic regression for estimating effect size.

Usage:
    RareRegressionAuto.R.jl -f file -p pheno -i id [-c covariates] [--om file] [--spa] [--firthp float] [--dc int]
    RareRegressionAuto.R.jl -h | --help | --version

Notes:
    1. Read gene burden scores from stdin, and output results to stdout.
    2. Load and use all of the covariates from -c file.

Options:
    -f file         Covariates file (including ID, Phenotype and covariates).
                        TSV with header, missing coded as NA.
    -p pheno        Phenotype name.
    -i id           Column name for individual id.
    --spa           Using SPAtest for binary trait, instead of GLM(default).
                        SPAtest optimized for unbalance data, can handle (qusi)-separation.
    --firthp float  P value cut-off for estimating effect size by Firth logistic regression
                        in SPAtest, default 0.001.
    -c covariates   Specifiy covariates, eg.: cov1,cov2|cov1.
    --om file       Output data matrix for regression analysis, default not.
    --dc int        Specify the start column index for data, default 7.
                        The first '--dc' columns will be copied to stdout.
    -h --help       Show this screen.
    --version       Show version.
"""

using DocOpt
# using GLM
using RCall
using DataFrames
using CSV
using Printf
using Distributions

args = docopt(doc, version=v"1.0")
# print(args)

path=args["-f"]
pheno   = args["-p"]
idname  = args["-i"]
covars  = args["-c"] == nothing ? []  : split(args["-c"],",")
SPAtest = args["--spa"]
firthP  = args["--firthp"] == nothing ? 0.001 : parse(Float64,args["--firthp"])
if SPAtest
    R"library(SPAtest)"
end
# println(covars)

df = CSV.File(path,header=1,delim='\t',missingstring="NA") |> DataFrame
# convert to string to do match later based on idname.
df[Symbol(idname)] = map(string,df[Symbol(idname)])

# check pheno and covariates whether in ped file.
check_names = vcat(pheno, idname, covars)
errors = [x for x in check_names if (Symbol(x) in Set(names(df))) == false]
# println(repr(errors))
# println(check_names)

# we want to load the rare variant score name as "RARESCORE" in DataFrame.
x = join(vcat("RARESCORE",covars),"+")
# y = Symbol(pheno)
# fm = Formula(y,Meta.parse(x))
fm = join([pheno, x], "~ 1 + ")
# as the glm model will auto add the interception.
# in order to make the log clear, we output "1 +" here.
println(stderr,"Regression model: "*fm)

# detect use binary or quantitative regression.
FAMILY = "binomial"
PHENO_LENGTH = length( Set([x for x in df[Symbol(pheno)] if ismissing(x)==false]) )
# https://www.statmethods.net/advstats/glm.html
if PHENO_LENGTH == 2
    println(stderr,"Using logistic regression as only detected 2 possible values for phenotype.")
else
    FAMILY = "gaussian"
    println(stderr,"Using gaussian(quantitative) regression as more than 2 possible values for phenotype.")
end

# println(covars)
#load data and to testing. first line as title for define ID, 2-n as meta and scores.
mytitle=[]
# Column start index for data.
datacol  = args["--dc"] == nothing ? 7 : parse(Int,args["--dc"])
# cache previous line result, if the same no need to compute again
# this is very useful for computing different maf together. As the result may exact the same.
cache_data = []
cache_out  = []

for line in eachline(stdin)
    # declear as global if we want to update a global values,
    # global values if view only for "for-loop".
    global mytitle
    global df
    global cache_data
    global cache_out

    ss = split(chomp(line))

    if length(mytitle) == 0
        # title = DataFrame(i=ss[datacol:end])
        # names!(title, Symbol(idname))
        mytitle = ss[datacol:end]
        if SPAtest == true
            out_title = vcat(ss[1:datacol], ["NS", "FRAC_WITH_RARE", "PVALUE", "BETA", "SEBETA", "WALD_PVALUE"])
        else
            out_title = vcat(ss[1:datacol], ["NS", "FRAC_WITH_RARE", "PVALUE", "BETA", "SEBETA", "ZSTAT"])
        end
        println(join(out_title,"\t"))

    else
        d = DataFrame()
        d[Symbol(idname)] = mytitle
        data_values = ss[datacol:end]
        # the same data, just copy result
        if data_values == cache_data
            println(join(vcat(ss[1:datacol],cache_out),"\t"))
            continue
        end
        cache_data = data_values

        d[Symbol("RARESCORE")] = [parse(Float64,x) for x in data_values]

        data = join(d, df, on=Symbol(idname))
        # println(first(data, 3))
        # println(data)
        if data == nothing
            println(stderr, "No common ids between RARESCORE from stdin and covariates file -f.")
            exit(-1)
        else
            score = data[Symbol("RARESCORE")]
            # the proportion of people with rare.
            f_rare = length([x for x in score if x != 0])*1.0/length(score)
            n_sample = length(score)
            # no information at grouping.
            out = []
            if f_rare == 0
                out = vcat(repr(n_sample), ["0", "NA", "NA", "NA", "NA"])
                cache_out = out
            else
                # print(SPAtest)
                if SPAtest == false
                    # fitting the GLM by R.
                    @rput data
                    # R"out = glm(as.formula($fm),data = data, family = binomial(link='logit'))"
                    #match epacts
                    #R"out = glm(as.formula($fm),data = data, family = binomial)"
                    R"out = glm(as.formula($fm),data = data, family = $FAMILY)"
                    # copy back the results
                    #              Estimate Std. Error     z value  Pr(>|z|)
                    #    (Intercept)  39.13214   21508.03  0.00181942 0.9985483
                    #    X           -19.56607   10754.01 -0.00181942 0.9985483
                    r = rcopy(R"matrix(coef(summary(out)),ncol=4)")
                    # 2×4 Array{Float64,2}:
                    # 39.1321  21508.0   0.00181942  0.998548
                    # -19.5661  10754.0  -0.00181942  0.998548

                    beta = r[2,1]
                    beta_se = r[2,2]
                    z = r[2,3]
                    p = r[2,4]

                    # fitting the GLM by Julia GLM.
                    # l = glm(fm, data, Binomial(), LogitLink())
                    # r = coeftable(l)
                    #
                    # beta = coeftable(l).cols[1][2]
                    # beta_se = coeftable(l).cols[2][2]
                    # z = coeftable(l).cols[3][2]
                    # # typeof(p) is StatsBase.PValue, using .v to access it value.
                    # # check an example of : x = StatsBase.PValue.(rand(1)).
                    # p = coeftable(l).cols[4][2].v

                    #out = vcat(ss[1:datacol], [n_sample, f_rare, p, beta, beta_se, z])
                    out = vcat(repr(n_sample), [@sprintf("%.4E", x) for x in [ f_rare, p, beta, beta_se, z]])
                else

                    if length(covars) == 0
                        tcov = ones(length(score))
                    else
                        tcov = data[map(Symbol,covars)]
                    end
                    tcov = data[map(Symbol,covars)]
                    # print(tcov)
                    # print(pheno)
                    # println(data[Symbol(pheno)])

                    tpheno = data[Symbol(pheno)]
                    # println(tpheno)
                    @rput(tpheno, score, tcov)
                    # # CALL fastSPA-0.1
                    r = rcopy(R"ScoreTest_SPA_wMeta(score,tpheno,tcov,minmac=1,Cutoff=0.1,output='metaZ',beta.out=T, beta.Cutoff = $firthP)")
                    # output is a dict
                    # :p_value     => 0.032222 , p value from SPAtest
                    # :p_value_NA  => 0.035015 , p value from traditional score test.
                    # :Is_converge => true
                    # :beta        => missing
                    # :SEbeta      => missing
                    # println(r[:beta])
                    beta =   ismissing(r[:beta])   ? "NA" : @sprintf("%.4E", r[:beta])
                    sebeta = ismissing(r[:SEbeta]) ? "NA" : @sprintf("%.4E", r[:SEbeta])
                    p = @sprintf("%.4E", r[:p_value])
                    wald_p = beta == "NA" ? "NA" : @sprintf("%.4E", ccdf(Normal(),abs(r[:beta]/r[:SEbeta]))*2) # Z to P, two sides.
                    out = vcat(repr(n_sample), [@sprintf("%.4E", f_rare), p, beta, sebeta, wald_p])
                end
                #output DataFrame for testing the association.
                if args["--om"] != nothing
                    #println(stderr, data)
                    CSV.write(args["--om"],data,delim="\t")
                end

            end
            cache_out = out
            println(join(vcat(ss[1:datacol],out),"\t"))

        end
    end
end

# Do logistic regression
# binary_phenotype ~ burden_score + covariates.
# Using R do the logistic regression, in order to match the results of epacts.

doc = """Regression for rare variant burden test.

phenotype(binary/quantitative) ~ burden_score + covariates.

*** Call R glm to do [logistic]regression, to match the results of epacts,
*** Slow compared to Julia version of GLM
*** Auto detect binary or quantitative regression based the the number of uniq number in phenotype.

Usage:
    RareRegressionAuto.R.jl -f file -p pheno -i id [-c covariates] [--om file]
    RareRegressionAuto.R.jl -h | --help | --version

Notes:
    1. Read gene burden scores from stdin, and output results to stdout.
    2. Load and use all of the covariates from -c file.

Options:
    -f file       Covariates file (including ID, Phenotype and covariates).
                    TSV with header, missing coded as NA.
    -p pheno      Phenotype name.
    -i id         Column name for individual id.
    -c covariates Specifiy covariates, eg.: cov1,cov2|cov1.
    --om file     Output data matrix for regression analysis, default not.
    -h --help     Show this screen.
    --version     Show version.
"""

using DocOpt
# using GLM
using RCall
using DataFrames
using CSV
using Printf

args = docopt(doc, version=v"1.0")
# print(args)

path=args["-f"]
pheno = args["-p"]
idname = args["-i"]
covars = args["-c"] == nothing ? []  : split(args["-c"],",")
# println(covars)

df = CSV.File(path,header=1,delim='\t',missingstring="NA") |> DataFrame
# convert to string to do match later based on idname.
df[Symbol(idname)] = map(string,df[Symbol(idname)])

# df = load(Stream(format"CSV", path),delim='\t',header_exists=true) |> DataFrame
# load(Stream(format"CSV", io)
# println(first(df, 3))

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
println(stderr,"GLM model: "*fm)

# detect use binary or quantitative regression.
FAMILY = "binomial"
PHENO_LENGTH = length( Set([x for x in df[Symbol(pheno)] if ismissing(x)==false]) )
# https://www.statmethods.net/advstats/glm.html
if PHENO_LENGTH == 2
    println(stderr,"Using logistic regression as only detect 2 possible values for phenotype.")
else
    FAMILY = "gaussian"
    println(stderr,"Using gaussian(quantitative) regression as more than 2 possible values for phenotype.")
end

#load data and to testing. first line as title for define ID, 2-n as meta and scores.
mytitle=[]
# println(length(mytitle))
datacol=7
for line in eachline(stdin)
    # declear as global if we want to update a global values,
    # global values if view only for "for-loop".
    global mytitle
    global df
    ss = split(chomp(line))

    if length(mytitle) == 0
        # title = DataFrame(i=ss[datacol:end])
        # names!(title, Symbol(idname))
        mytitle = ss[datacol:end]
        out = vcat(ss[1:datacol], ["NS", "FRAC_WITH_RARE", "PVALUE", "BETA", "SEBETA", "ZSTAT"])
        println(join(out,"\t"))

    else
        d = DataFrame()
        d[Symbol(idname)] = mytitle
        d[Symbol("RARESCORE")] = [parse(Float64,x) for x in ss[datacol:end]]

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
            if f_rare == 0
                out = vcat(ss[1:datacol], repr(n_sample), ["0", "NA", "NA", "NA", "NA"])
            else
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
                out = vcat(ss[1:datacol], repr(n_sample), [@sprintf("%.4E", x) for x in [ f_rare, p, beta, beta_se, z]])
                #output DataFrame for testing the association.
                if args["--om"] != nothing
                    #println(stderr, data)
                    CSV.write(args["--om"],data,delim="\t")
                end

            end
            println(join(out,"\t"))

        end
    end
end

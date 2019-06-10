##################################################################
## GENE-LEVEL BURDEN TEST
## INPUT VARIABLES: 
##   n        : total # of individuals
##   genos    : genotype matrix for each gene
##   NS       : number of called samples for each marker
##   AC       : allele count for each marker
##   MAC      : minor allele count for each marker
##   MAF      : minor allele frequency
##   vids     : indices from 1:n after AF/AC threshold
## EXPECTED OUTPUT : list(p, addcols, addnames) for each genos row
##   p        : p-value
##   add      : additional column values
##   cname    : additional column names
##################################################################

## gene.b.collapse() : Collapsing test using logistic regression
## KEY FEATURES : 0/1 collapsing variables to dichotomous traits
## TRAITS : BINARY
## RETURNS : PVALUE, BETA, SEBETA, CHISQ
## MISSING VALUE : IMPUTED AS MAJOR ALLELES

##################################################################
## HELPER FUNCTION:  fast.logistf.fit
## This function calculates the maximum penalized likelihood estimates
##################################################################
 
fast.logistf.fit <- function (x, y, weight = NULL, offset = NULL, firth = TRUE, col.fit = NULL, 
    init = NULL, control) {
    n <- nrow(x)
    k <- ncol(x)
    if (is.null(init)) 
        init = rep(0, k)
    if (is.null(col.fit)) 
        col.fit = 1:k
    if (is.null(offset)) 
        offset = rep(0, n)
    if (is.null(weight)) 
        weight = rep(1, n)
    if (col.fit[1] == 0) 
        maxit <- 0
    if (missing(control)) 
        control <- fast.logistf.control()
    maxit <- control$maxit
    maxstep <- control$maxstep
    maxhs <- control$maxhs
    lconv <- control$lconv
    gconv <- control$gconv
    xconv <- control$xconv
    beta <- init
    l.change <- 5
    iter <- 0
    pi <- as.vector(1/(1 + exp(-x %*% beta - offset)))
    loglik <- sum(weight[y == 1] * log(pi[y == 1])) + sum(weight[y == 
        0] * log(1 - pi[y == 0]))
    if (firth) {
        #XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)
        XW2 <- t(x * (weight * pi * (1-pi))^0.5)
        Fisher <- crossprod(t(XW2))
        loglik <- loglik + 0.5 * determinant(Fisher)$modulus[1]
    }
    evals <- 1
    repeat {
        loglik.old <- loglik
        beta.old <- beta
        #XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)
        XW2 <- t(x * (weight * pi * (1-pi))^0.5)
        #Fisher <- crossprod(t(XW2))
        #covs <- fast.invFisher(Fisher)
        #H <- crossprod(XW2, covs) %*% XW2
        myQR <- qr(t(XW2))
        Q <- qr.Q(myQR)
        h <- (Q*Q) %*% rep(1, ncol(Q))        
        if (firth) 
            #U.star <- crossprod(x, weight * (y - pi) + diag(H) * (0.5 - pi))
            U.star <- crossprod(x, weight * (y - pi) + h * (0.5 - pi))
        else U.star <- crossprod(x, weight * (y - pi))
        XX.covs <- matrix(0, k, k)
        if (col.fit[1] != 0) {
            #XX.XW2 <- crossprod(x[, col.fit, drop = FALSE], diag(weight * pi * (1 - pi))^0.5)
            XX.XW2 <- t(x[, col.fit, drop=FALSE] * (weight * pi * (1-pi))^0.5)
            XX.Fisher <- crossprod(t(XX.XW2))
            XX.covs[col.fit, col.fit] <- fast.invFisher(XX.Fisher)   ###### HERE IS THE PROBLEM!!!
        }
        if(all(is.na(XX.covs)) == T) {
            break
        }  
        delta <- as.vector(XX.covs %*% U.star)
        delta[is.na(delta)] <- 0
        mx <- max(abs(delta))/maxstep
        if (mx > 1) 
            delta <- delta/mx
        evals <- evals + 1
        if (maxit > 0) {
            iter <- iter + 1
            beta <- beta + delta
            for (halfs in 1:maxhs) {
                pi <- as.vector(1/(1 + exp(-x %*% beta - offset)))
                loglik <- sum(weight[y == 1] * log(pi[y == 1])) + 
                  sum(weight[y == 0] * log(1 - pi[y == 0]))
                if (firth) {
                  #XW2 <- crossprod(x, diag(weight * pi * (1 - pi))^0.5)
                  XW2 <- t(x * (weight * pi * (1-pi))^0.5)  
                  Fisher <- crossprod(t(XW2))
                  loglik <- loglik + 0.5 * fast.logDet(Fisher)
                }
                if(is.na(loglik)==T) {
                    break
                }
                evals <- evals + 1
                l.change <- loglik - loglik.old
                if (loglik > loglik.old) 
                  break
                beta <- beta - delta * 2^(-halfs)
            }
            if(is.na(loglik)==T) {
                break
            }
        }
        if(is.na(loglik)==T) {
            break
        }
        if (iter == maxit | ((max(abs(delta)) <= xconv) & (all(abs(U.star[col.fit]) < 
            gconv)) & (all(l.change < lconv)))) 
            break
    }
    # Error catching (if chol(x) not positive definite)
    if(all(is.na(XX.covs))==T | is.na(loglik)==T) {
        var <- XX.covs
        list(beta = NA, var = var, pi = NA, hat.diag = NA, 
        loglik = NA, iter = NA, evals = NA, conv = c(NA, 
            NA, NA))
    } else {
        var <- XX.covs
        list(beta = beta, var = var, pi = pi, hat.diag = h, 
            loglik = loglik, iter = iter, evals = evals, conv = c(l.change, 
                max(abs(U.star)), max(abs(delta))))
    }
}


##################################################################
## HELPER FUNCTION:  fast.logistf.control
## This function provides the convergence control parameters
##################################################################
fast.logistf.control <- function (maxit = 50, maxhs = 15, maxstep = 15, lconv = 1e-05, 
    gconv = 1e-05, xconv = 1e-05) 
{
    list(maxit = maxit, maxhs = maxhs, maxstep = maxstep, lconv = lconv, 
        gconv = gconv, xconv = xconv)
}

##################################################################
## HELPER FUNCTION:  fast.logDet
##################################################################
fast.logDet <- function (x) {
    my.chol <- tryCatch(chol(x),error=function(e) {NA})
    if (all(is.na(my.chol))==T) {
        return(NA)
    } else {
        return (2 * sum(log(diag(my.chol))))
    }
}   

##################################################################
## HELPER FUNCTION:  fast.invFisher
##################################################################
fast.invFisher <- function(x) {
  my.chol <- tryCatch(chol(x),error=function(e) {NA})
    if (all(is.na(my.chol))==T) {
        return(NA)
    } else {
        return (chol2inv(my.chol))
    }
  #ifelse(is.na(my.chol), NA, chol2inv(my.chol))
}

group.b.burdenFirth <- function() {

  m <- nrow(genos)
  cname <- c("BETA","SEBETA","CHISQ") #,"N.CASE","N.CTRL","AF.CASE","AF.CTRL")

  if ( m > 0 ) {

    p <- rep(NA,1) # store p-values for m markers
    k <- ncol(cov)
    add <- matrix(NA,1,3) # extra columns:  Beta, SE, Chisq

    #BEGIN JASON 
    #genos <- as.double(colSums(genos,na.rm=T) > 0)
    if (!is.null(weights)) {
      genos <- apply(genos * weights, 2, sum, na.rm=T)
    } else {
      genos <- apply(genos, 2, sum, na.rm=T)
    }
    #sg <- sum(genos)
    sg <- sum(genos, na.rm=T)
    #END JASON 

    if ( ( sg > 0 ) && ( sg < n ) ) {

        #print (rownames(genos)[i])
        #print(vMAF[i]);
        # Full model pheno ~ intercept + geno + cov
        fit.full <- fast.logistf.fit(x = cbind(t(t(genos)), cov), y=pheno)
        # Reduced model pheno ~ intercept + geno + cov (omit genotype using "col.fit")
        # dimensions of x:  intercept(1) + geno(1) + cov(k) = k+2 columns
        if(all(is.na(fit.full$beta))==F) {
            fit <- fast.logistf.fit(x = cbind(t(t(genos)), cov), col.fit = (1:(k+1))[-1], y=pheno)
            p <- 1 - pchisq(2 * (fit.full$loglik - fit$loglik), 1)
            add[1:3] <- c(fit.full$beta[1], sqrt(fit.full$var[1,1]), 2*(fit.full$loglik - fit$loglik))
        } else {
            print("Design matrix X is not positive definite!")
            print(paste("Skipping SNP"))
            p <- NA
            add[1,1:3] <- rep(NA,3)
        }
      }


      return(list(p=p,
                  add=add,
                  cname=cname))
  }
  return(list(p=NA,add=rep(NA,length(cname)),cname=cname))
}

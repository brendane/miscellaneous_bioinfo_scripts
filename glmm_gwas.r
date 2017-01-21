### The function in this file is based on GMMAT::glmm.wald with some
### modifications to allow interactions and stripping away of a lot
### of messy bits of code that I don't think I need. I also took some
### ideas from the lrgpr package.
###
### It will work well with the data reading functions from lrgpr that
### store genotype information in a big matrix format on disk to save
### memory.
###
### It always omits individuals with missing data, and it only handles
### one random effects matrix.
###
### The various fitting settings available are the same as glmm.wald, but
### I removed the code that tells you if you've specified something wrong.
###
### Args:
###   fixed       A standard formula giving the response and fixed effects,
###               where the symbol SNP is used as a placeholder for the SNP
###               being tested.  e.g.: y ~ age + sex + SNP
###   X           A matrix of genotypes, where SNPs are columns and
###               individuals are rows.
###   snpidx      The index (in X) for the SNP to be tested.
###   kmat        The kinship matrix - same as "GRM" argument for
###               glmm.wald. This has to be in the same order as the vectors
###               in fixed.
###   Terms       The terms on which to perform the Wald test - should be
###               in the same order as in fixed. Interaction terms end up
###               at the end, and the intercept is the first term. Note
###               that the Wald test is a test of all the terms - so if
###               you specify multiple terms, you still get one p-value.
###
###   All other options are from glmm.wald and the same defaults and
###   limits.
###
### Written by Brendan Epstein 2016-08-03.
###

glmm_gwas = function(fixed, X, snpidx, kmat, Terms, family, method='REML',
                     method.optim='AI', maxiter=500, tol=1e-5, taumin=1e-5,
                     taumax=1e5, tauregion=10, verbose=FALSE, ...) {

    require(aod)
    require(GMMAT)

    ## SNP is a column of the matrix X, which can be an attached
    ## big.matrix. We make sure it is numeric here.
    SNP = as.numeric(X[, snpidx])

    ## Convert fixed to formula.
    fixed = as.formula(deparse(fixed))

    ## Then subset the kinship matrix as needed to deal with missing
    ## data and count the number of individuals left (N).
    idx = match(rownames(model.frame(formula=fixed, na.action=na.omit)),
                rownames(model.frame(formula=fixed, na.action=na.pass)))
    tmpkmat = kmat[idx, idx]
    N = nrow(tmpkmat)
    if(method.optim != 'Brent')
        tmpkmat = list(kins1=tmpkmat)

    ## First do a GLM without the kinship matrix.
    ## It would be nice to be able to use lrgpr::glmApply here to speed
    ## up the calculations, but it only returns the p-values.
    fit0 = glm(formula=fixed, family=family, ...)

    ## And then use glmmkin.fit to include the K matrix.
    fit <- try(GMMAT:::glmmkin.fit(fit0, tmpkmat, method=method,
                                   method.optim=method.optim,
                                   maxiter=maxiter,
                                   tol=tol,
                                   taumin=taumin,
                                   taumax=taumax,
                                   tauregion=tauregion,
                                   verbose=verbose))

    ## Do a Wald test on the specified terms. The terms are tested
    ## simultaneously, so you don't get multiple p-values for multiple
    ## terms. The hypothesis being tested is that the coefficients are
    ## not equal to zero.
    term_names = names(coef(fit0))[Terms]
    Beta = numeric(length(Terms)) * NaN
    se   = numeric(length(Terms)) * NaN
    if(class(fit) != 'try-error') {
        convgd = fit$converged
        Beta   = fit$coefficients[Terms]
        se     = sqrt(diag(fit$cov)[Terms])
        ## This line (commented out) is the code from GMMAT::glmm.wald.
        ## I replaced it with aod::wald.test to deal with multiple terms.
        ## aod::wald.test seems to give the same results for one term.
        ## pval   = pchisq((Beta / se)^2, 1, lower.tail=FALSE)
        pval   = wald.test(fit$cov, fit$coefficients, Terms=Terms)$result$chi2['P']
    } else {
        pval   = NaN
        convgd = NA
    }

    names(Beta) = term_names
    names(se) = term_names
    ## Return a list of useful stuff.
    cn = colnames(X)
    if(!is.null(cn)) {
        snp = cn[snpidx]
    } else {
        snp = NA
    }
    list(beta=Beta, se=se, pval=pval, converged=convgd, N=N,
         terms=term_names, snp=snp)
}

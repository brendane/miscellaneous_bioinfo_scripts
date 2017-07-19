#!/usr/bin/env Rscript
###
### Use MASS::stepAIC to perform backward selection on variants.
###
### Imputes missing genotype data by using the mean; ignores strains with
### missing phenotype data.
###

library(data.table)
library(MASS)

argv = commandArgs(trailingOnly=TRUE)
genotype_file = argv[1]
phenotype_file = argv[2]
variant_file = argv[3]
n_vars = as.numeric(unlist(strsplit(argv[4], ',')))
min_maf = as.numeric(argv[5])
output_prefix = argv[6]

oh = file(paste0(output_prefix, '.stats.tsv'), 'wb')
ohf = file(paste0(output_prefix, '.fitted.tsv'), 'wb')

n_random_vars = 10000

genotypes = fread(genotype_file, header=TRUE)
setkey(genotypes, 'rs')
phenotypes = fread(phenotype_file)
variants = scan(variant_file, what='character')

maf = apply(genotypes[, -(1:3)], 1,
            function(x) min(sum(x == 0, na.rm=TRUE), sum(x == 1, na.rm=TRUE))/sum(!is.na(x)))
genotypes = genotypes[maf >= min_maf, ]
cat(nrow(genotypes), 'variants left after MAF filtering\n', file=stderr())

if(length(variants) < max(n_vars)) {
    stop("Can't choose that many variants")
}

X = apply(t(as.matrix(genotypes[variants, -(1:3)])), 2,
          function(x) { x[is.na(x)] = mean(x, na.rm=TRUE); x})
p = phenotypes[match(rownames(X), as.character(phenotypes[[ 2]])), ]
y = p[[3]]

X_rand = t(as.matrix(genotypes[sample(1:nrow(genotypes), n_random_vars, FALSE), -(1:3)]))

cat('n_vars_start\tr2\tr2_unadj\tgenomic_r2\tgenomic_r2_unadj\tn_vars_final\t',
    'ranks_final\tvars_final\tcoefficients\n',
    file=oh, sep='')
cat('n_vars_start\tphenotype\tfitted\n', file=ohf)
for(nv in n_vars) {
    rhs = paste('~', paste(paste0('X[,', 1:nv, ']'), collapse=' + '))
    form = paste('y', rhs)
    m_full = lm(form)
    m_final = stepAIC(m_full, scope=list(lower=~1), direction='backward',
                      trace=FALSE)
    s = summary(m_final)
    r2 = s$r.squared
    r2a = s$adj.r.squared
    rhs_final = as.character(m_final$call$formula)[3]
    ranks_final = as.numeric(gsub('X\\[,', '',
                                  gsub('\\]', '', names(m_final$coefficients)[-1])))
    coefs_final = m_final$coefficients
    variants_final = variants[ranks_final]
    ft = m_final$fitted.values
    gr2 = 0
    gr2a = 0
    for(i in 1:n_random_vars) {
        ss = summary(lm(paste('X_rand[, i] ~', rhs_final)))
        gr2 = gr2 + ss$r.squared
        gr2a = gr2a + ss$adj.r.squared
    }
    cat(nv, '\t', round(r2a, 3), '\t', round(r2, 3),
        '\t', round(gr2a/n_random_vars, 3),
        '\t', round(gr2/n_random_vars, 3),
        '\t', length(variants_final),
        '\t', paste(ranks_final, collapse=','),
        '\t', paste(variants_final, collapse=','),
        '\t', paste(coefs_final, collapse=','),
        '\n', sep='', file=oh)
    write.table(cbind(rep(nv, sum(!is.na(y))), y[!is.na(y)], ft),
                col.names=FALSE, row.names=FALSE, quote=FALSE,
                sep='\t', file=ohf)
}

close(oh)
close(ohf)

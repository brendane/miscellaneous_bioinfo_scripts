#!/usr/bin/env Rscript
###
### R script to do run GWAS on pairwise interactions. Choice
### of logistic or gaussian regression, and optional covariates.
###
### DETAILS
###
### For the initial pass, it uses the default settings of lrgpr::glmApply,
### which set missing genotypes to the mean and uses a Wald Chi-square
### test. It also tests all the SNPs - the filters are applied before
### testing interactions.
###
### The second pass also uses lrgpr::glmApply for speed.
###
#==============================================================================#

library(lrgpr)
library(optparse)

#==============================================================================#

opts = list(make_option('--output'),
            make_option('--test-percent', type='double', default=1),
            make_option('--min-maf', type='double', default=0))
parser = OptionParser(option_list=opts)
argv = parse_args(parser, positional_arguments=TRUE)

outprefix = argv$options[['output']]
test_percent = argv$options[['test-percent']]
min_maf = argv$options[['min-maf']]
tfam_file = argv$args[1]
tped_file = argv$args[2]
model_type = argv$args[3]
if(length(argv$args) > 3) {
    covariates_file = argv$args[4]
    covariates = read.table(covariates_file)
    use_covariates = TRUE
} else {
    use_covariates = FALSE
}

## Figure out the type of regression to run
if(model_type == 'gaussian') {
    fam = gaussian(link='identity')
} else if(model_type == 'logistic') {
    fam = binomial(link='logit')
} else {
    stop('Model type not recognized')
}

## Make the binary file
convertToBinary(tped_file, './binary_genotypes', 'TPED')

## Read in data
X = attach.big.matrix('./binary_genotypes_descr', readonly=TRUE)
y = read.table(tfam_file)[, 6]

## Get minor allele frequency and other information
maf = getAlleleFreq(X)
maf[maf > 0.5] = 1 - maf[maf > 0.5]
pos = as.numeric(gsub('.+-', '', colnames(X)))
chr = gsub('-.+', '', colnames(X))

## Set up formula for initial pass
if(use_covariates) {
    f = paste('y ~ ', paste(paste0('covariates[, ', 1:ncol(covariates), ']'),
                            collapse=' + '), '+ SNP')
    f = as.formula(f)
    term = ncol(covariates) + 2
} else {
    f = as.formula('y ~ SNP')
    term = 2
}

## Run initial pass - one SNP at a time
initial_pvalues = glmApply(f, features=X,  terms=term, family=fam)$pValues[, 1]

## Apply filters
snps_to_test = initial_pvalues <= quantile(initial_pvalues,
                                           test_percent, na.rm=TRUE)
snps_to_test[is.na(snps_to_test)] = FALSE
snps_to_test = snps_to_test &
                maf >= min_maf
snps_to_test[is.na(snps_to_test)] = FALSE
cat('SNPs to test <=', sum(snps_to_test), '\n',
    file=stderr())

write.table(data.frame(chr=chr, pos=pos, snp=colnames(X), freq=maf,
                       interact=as.numeric(snps_to_test),
                       p=initial_pvalues),
            file=paste0(outprefix, '.single.tsv'),
            sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)

## Run interaction tests
if(use_covariates) {
    f_cov = paste(paste(paste0('covariates[, ', 1:ncol(covariates), ']'),
                        collapse=' + '), '+')
} else {
    f_cov = ''
}
f = paste('y ~', f_cov, 'Xi * SNP')

outhandle = file(paste0(outprefix, '.interactions.tsv'), 'wb')
cat('chr1\tpos1\tsnp1\tfreq1\t',
    'chr2\tpos2\tsnp2\tfreq2\t',
    'pi\tp1\tp2\n', sep='', file=outhandle)

for(i in 1:(ncol(X) - 1)) {

    if(!snps_to_test[i])
        next
    stt = snps_to_test
    stt[i] = FALSE

    chr1 = chr[i]
    pos1 = pos[i]
    maf1 = maf[i]
    snp1 = colnames(X)[i]

    Xi = set_missing_to_mean(X[, i])

    p = glmApply(f, features=X, terms=term + 2, family=fam,
                 cincl=which(stt))$pValues[, 1]

    output = data.frame(rep(chr1, sum(stt & !is.na(p))),
                        rep(pos1, sum(stt & !is.na(p))),
                        rep(snp1, sum(stt & !is.na(p))),
                        rep(maf1, sum(stt & !is.na(p))),
                        chr[stt & !is.na(p)],
                        pos[stt & !is.na(p)],
                        colnames(X)[stt & !is.na(p)],
                        maf[stt & !is.na(p)],
                        p[!is.na(p)],
                        rep(initial_pvalues[i], sum(stt & !is.na(p))),
                        initial_pvalues[stt & !is.na(p)])
    write.table(output, file=outhandle, col.names=FALSE,
                row.names=FALSE, quote=FALSE, sep='\t')
}
close(outhandle)

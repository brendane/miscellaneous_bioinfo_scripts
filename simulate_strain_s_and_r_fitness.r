#!/usr/bin/env Rscript
#
# Randomly draw distributions of fitness effects for variants and use those
# to simulate S&R evolution.
#
# Can either specify "all" to draw fitness effects from an exponential
# distribution, or a number to specify the number of SNPs that will have
# an effect (all others set to zero).
#

library(data.table)
library(dplyr)
library(dtplyr)

## Program options:
cargs = commandArgs(trailingOnly=TRUE)
genotype_file = cargs[1]
n_snps = cargs[2]
n_generations = as.numeric(cargs[3])
min_maf = as.numeric(cargs[4])
max_miss = as.numeric(cargs[5])
output_prefix = cargs[6]

## Probability that the reference allele is more fit than the alternate
ref_allele_more_fit = 0.8
if(length(cargs) > 6) {
    ref_allele_more_fit = as.numeric(cargs[7])
}


## Read in data
genotypes = fread(genotype_file, header=TRUE) %>%
    select(-pos, -contig, -rs) %>%
    as.matrix()

genotype_names = fread(genotype_file, header=TRUE) %>%
    select(pos, contig, rs)

initial_allele_freqs = apply(genotypes, 1, function(x) {
                             ref = sum(!(is.na(x)) & x == 1)
                             alt = sum(!(is.na(x)) & x == 0)
                             ref / (ref + alt)
})


## Filter data
maf = ifelse(initial_allele_freqs < 0.5, initial_allele_freqs,
             1 - initial_allele_freqs)
missingness = apply(genotypes, 1, function(x) sum(is.na(x)) / length(x))

keep = maf > min_maf & missingness < max_miss
genotypes = genotypes[keep, ]
genotype_names = genotype_names[keep, ]
initial_allele_freqs = initial_allele_freqs[keep]
maf = maf[keep]
missingness = missingness[keep]


## Draw fitness effects
if(!is.na(as.numeric(n_snps))) {
    ## Set a few SNPs to fitness effect of 1, all others = 0
    ns = as.numeric(n_snps)
    x = rep(0, nrow(genotypes))
    x[sample(1:nrow(genotypes), ns, FALSE)] = sample(c(1, -1), ns, TRUE,
                                                     prob=c(1-ref_allele_more_fit, ref_allele_more_fit))
    snp_fitness_effects = x
} else {
    ## Draw from exponential distribution for all SNPs
    snp_fitness_effects = rexp(nrow(genotypes)) * ifelse(sample(c(TRUE, FALSE), nrow(genotypes), TRUE,
                                                                prob=c(1-ref_allele_more_fit, ref_allele_more_fit)),
                                                         -1, 1)
}


## Simulate fitness of strains and calculate allele frequency and
## allele frequency change. Takes fitness effects as Wrightian (average
## number of offspring).

## 1) Column-by-column, multiply the fitness_effect of the alternate
##    allele by whether the alternate allele is present. Do this in
##    a way that doesn't give any negative numbers (if alternate allele
##    effect is negative, make it a positive effect of the ref. allele).
## 2) Row-by-row, set missing genotypes to have the average effect
## 3) Add up all the fitness effects in a strain and divide by the
##    mean across strains.
## Note that in the genotype input file 1 = reference, 0 = alternate
snp_effects =
    apply(genotypes, 2, function(x) {
          snp_fitness_effects * ifelse(snp_fitness_effects < 0, x, 1-x)
         }) %>%
    apply(., 1, function(x) {
          mw = mean(x, na.rm=TRUE)
          ifelse(is.na(x), mw, x)
         }) %>%
    t()
fitnesses = colSums(snp_effects) / mean(colSums(snp_effects))

## Calculate final strain frequencies assuming equal starting frequencies.
## Uses Wrightian fitness (fitness = average number of offspring).
initial_strain_freqs = rep(1/length(fitnesses[[1]]), length(fitnesses[[1]]))
n_inds = (fitnesses^n_generations)
strain_freqs = n_inds / sum(n_inds)
relative_fitnesses = log2(strain_freqs / initial_strain_freqs)
## If fitness = 0, need to do a rough correction:
relative_fitnesses[is.na(relative_fitnesses) | is.infinite(relative_fitnesses)] =
    min(relative_fitnesses[!is.na(relative_fitnesses) & !(is.infinite(relative_fitnesses))]) - 2


## Calculate final allele frequencies
stopifnot(all(names(strain_freqs) == colnames(genotypes)))
allele_freqs = apply(genotypes, 1, function(x) {
                     ref = sum(strain_freqs[!(is.na(x)) & x == 1])
                     alt = sum(strain_freqs[!(is.na(x)) & x == 0])
                     ref / (ref + alt)
         })


## Write output:

##  1. Phenotype file for Plink with fitness of strains
##     (for external association testing)
write.table(cbind(names(fitnesses), names(fitnesses), fitnesses),
            file=paste0(output_prefix, '.fitness.txt'),
            sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)

##  2. Phenotype file with relative fitness
write.table(cbind(names(relative_fitnesses), names(relative_fitnesses), relative_fitnesses),
            file=paste0(output_prefix, '.relative_fitness.txt'),
            sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)

##  3. Table of fitness effects and allele frequency start and end for
##     each SNP
write.table(data.frame(genotype_names, 'effect'=snp_fitness_effects,
                       'initial'=initial_allele_freqs,
                       'final'=allele_freqs),
            file=paste0(output_prefix, '.snp_effects.txt'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

##  4. Table of strain frequencies and fitnesses
write.table(data.frame(strain=names(fitnesses), 'fitness'=fitnesses,
                       'relative_fitness'=relative_fitnesses,
                       'initial'=initial_strain_freqs,
                       'final'=strain_freqs),
            file=paste0(output_prefix, '.strain_frequencies.txt'),
            sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

## 5. Phenotype file with change in strain frequency
write.table(cbind(names(strain_freqs), names(strain_freqs),
                  strain_freqs - initial_strain_freqs),
            file=paste0(output_prefix, '.strain_freqs.txt'),
            sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)

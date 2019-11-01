#!/usr/bin/env Rscript
#
# This script compares various pairwise gene divergence statistics for
# a set of genes to the other genes in the genome.
#

library(ape)
library(magrittr)
library(optparse)

rcv = function(f, ...) {
    read.csv(f, sep='\t', header=TRUE, comment.char='#', as.is=TRUE, check.names=FALSE, ...)
}

optlist = list(
               make_option('--output'),
               make_option('--n-random-draws', type='double',
                           help='number of random draws when comparing to background'),
               make_option('--species',
                           help='tsv file with species names; columns named "species" and "strain"'),
               make_option('--strains',
                           help='text file with one strain per line'),
               make_option('--targets',
                           help='text file with list of genes to target; if not given, use symbiotic column'),
               make_option('--random-targets',
                           help='tsv file with lists of background genes, one list per column; if not given, background genes are matched on other characteristics'),
               make_option('--included-targets',
                           help='text file with list of genes to limit the analysis to'),
               make_option('--kaks',
                           help='pairwise Ka/Ks values'),
               make_option('--gene-id-column', default='subset',
                           help='column of input file to use as gene IDs; defaults to "subset"; note that most analyses assume that each row is a gene'),

               make_option('--count-genes', type='logical', action='store_true', default=FALSE,
                           help='Make a table of gene counts'),
               make_option('--compare-strain-count', type='logical', action='store_true', default=FALSE,
                           help='Make plots comparing number of strains in targets and background'),
               make_option('--compare-copy-number', type='logical', action='store_true', default=FALSE,
                           help='Make plots comparing number of copies in targets and background'),
               make_option('--compare-protein-divergence', type='logical', action='store_true', default=FALSE,
                           help='Make plots comparing pairwise divergence in targets and background'),
               make_option('--compare-relative-protein-divergence', type='logical', action='store_true', default=FALSE,
                           help='Make plots comparing relative pairwise divergence in targets and background'),
               make_option('--compare-paralog-protein-divergence', type='logical', action='store_true', default=FALSE,
                           help='Make plots comparing divergence of paralogs in targets and background'),
               make_option('--compare-kaks', type='logical', action='store_true', default=FALSE,
                           help='Make plots comparing pairwise Ka/Ks in targets and background'),
               make_option('--compare-paralog-kaks', type='logical', action='store_true', default=FALSE,
                           help='Make plots comparing pairwise Ka/Ks of paralogs in targets and background'),
               make_option('--compare-expected-divergence', type='logical', action='store_true', default=FALSE,
                           help='Make plots comparing phylogenetic distribution of targets and background'),

               make_option('--protein-divergence-sampling-tolerance', type='double',
                           help='Tolerance for matching random gene sets in number of strains and exp. divergence (non-paralog prot. div, Ka/Ks)'),
               make_option('--paralog-divergence-sampling-tolerance', type='double',
                           help='Tolerance for matching random gene sets in number of sequences and copy-number for paralog comparison')
               )


## Parse options
parser = OptionParser(option_list=optlist)
opts = parse_args(parser, positional_arguments=TRUE)

prot_div_tol = opts[['options']][['protein-divergence-sampling-tolerance']]
para_cn_ng_tol = opts[['options']][['paralog-divergence-sampling-tolerance']]

orthofile = opts[['args']][1]   # tsv file with one line per gene set and summaries of annotation and divergence stats 'table/orthosets_data.73strains_alpha_complete.nopseudo.2019-07-31.tsv'
outpre = opts[['options']][['output']]
if(is.null(opts[['options']][['random-targets']])) {
    background = NULL
    N = opts[['options']][['n-random-draws']]
} else {
    background = read.csv(opts[['options']][['random-targets']],
                          sep='\t', header=FALSE, as.is=TRUE)
    N = ncol(background)
}

strain_file = opts[['options']][['strains']]
spp_file = opts[['options']][['species']] # ('table/73strains_alpha_complete.species.tsv')
ka_ks_file = opts[['options']][['kaks']] # ('table/orthoset_pairwise_distances_2019-07-11.kaks.nopseudo.tsv')
targets_file = opts[['options']][['targets']]
gene_column = opts[['options']][['gene-id-column']]

if(opts[['options']][['compare-paralog-kaks']]) opts[['options']][['compare-kaks']] = TRUE


## Read in data
osets = rcv(orthofile)
if(!is.null(opts[['options']][['included-targets']])) {
    inc = scan(opts[['options']][['included-targets']], what='character', sep='\n')
    osets = osets[osets[, gene_column] %in% inc, ]
}
strains = scan(strain_file, what='character', sep='\n')
s = rcv(spp_file)
species = structure(s[, 'species'], names=s[, 'strain'])
pairwise_kaks = rcv(ka_ks_file)
if(is.null(targets_file)) {
    targeted = osets[, 'symbiotic'] == 1
} else {
    targeted = osets[, gene_column] %in% scan(targets_file, what='character', sep='\n')
}


if(opts[['options']][['count-genes']]) {
    ## Count of genes
    counts = data.frame(
                        'Measure'=c('Total genes',
                                    'Target genes',
                                    'Single copy genes',
                                    'Core genes',
                                    'Single-copy core genes',
                                    'Mean no. strains',
                                    'Mean no. genes',
                                    'Unique to one strain',
                                    'Single copy target genes',
                                    'Core target genes',
                                    'Single copy core target genes',
                                    'Unique target genes'),

                        'Counts'=c(nrow(osets),
                                   sum(targeted),
                                   sum(osets[, 'single_copy']),
                                   sum(osets[, 'core']),
                                   sum(osets[, 'core'] == 1 & osets[, 'single_copy'] == 1),
                                   round(mean(osets[, 'n_strains']), 2),
                                   round(mean(osets[, 'n_genes']), 2),
                                   sum(osets[, 'n_strains'] == 1),
                                   sum(targeted & osets[, 'single_copy'] == 1),
                                   sum(targeted & osets[, 'core'] == 1),
                                   sum(targeted & osets[, 'single_copy'] == 1 & osets[, 'core'] == 1),
                                   sum(targeted & osets[, 'n_strains'] == 1)),
                        check.names=FALSE)
    write.table(counts, file=paste0(outpre, '.counts.tsv'),
                sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)
}


if(opts[['options']][['compare-strain-count']]) {
    ## Compare number of strains in targeted genes to the number of
    ## strains in other genes.
    pdf(paste0(outpre, '.strain_count_comparison.pdf'), width=8, height=8, useDingbats=FALSE)
    op = par(no.readonly=TRUE)

    ## Distribution of number of strains per orthoset
    par(mfrow=c(2, 1))
    hist(osets[, 'n_strains'],
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='All genes',
         ylab='Number of genes',
         xlab='', breaks=seq(0, length(strains), 1),
         mar=c(2.1, 4.1, 3.1, 2.1))
    s = summary(osets[, 'n_strains'])
    text(length(strains)/2, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    hist(osets[targeted, 'n_strains'],
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='Targets',
         ylab='Number of genes',
         xlab='Number of strains', breaks=seq(0, length(strains), 1),
         mar=c(4.1, 4.1, 3.1, 2.1))
    s = summary(osets[targeted, 'n_strains'])
    text(length(strains)/2, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))

    par(op)


    ## Compare target genes to random draws
    ## Draw ngenes from background without replacement. No restrictions on
    ## which background genes to draw from.
    par(mfrow=c(4, 2))
    ngenes = sum(targeted)
    rands = matrix(nrow=N, ncol=4)
    colnames(rands) = c('p25', 'median', 'mean', 'p75')
    deciles = matrix(nrow=N, ncol=10)
    for(i in 1:N) {
        if(is.null(background)) {
            r = sample(osets[!targeted, 'n_strains'], ngenes, FALSE)
        } else {
            r = osets[osets[, gene_column] %in% background[, i], 'n_strains']
        }
        rands[i, 'p25'] = quantile(r, 0.25)
        rands[i, 'p75'] = quantile(r, 0.75)
        rands[i, 'mean'] = mean(r)
        rands[i, 'median'] = median(r)
        deciles[i, ] = quantile(r, seq(0.1, 1, 0.1))
    }

    par(mfcol=c(4, 2))
    hist(rands[, 'p25'], col='gray', border='white', xlab='# Strains',
         ylab='# genes', main='25th percentile', breaks=seq(0, length(strains), 1),
         xaxs='i', yaxs='i')
    abline(v=quantile(osets[targeted, 'n_strains'], 0.25))
    hist(rands[, 'median'], col='gray', border='white', xlab='# Strains',
         ylab='# genes', main='median', breaks=seq(0, length(strains), 1),
         xaxs='i', yaxs='i')
    abline(v=median(osets[targeted, 'n_strains']))
    hist(rands[, 'mean'], col='gray', border='white', xlab='# Strains',
         ylab='# genes', main='mean', breaks=seq(0, length(strains), 1),
         xaxs='i', yaxs='i')
    abline(v=mean(osets[targeted, 'n_strains']))
    hist(rands[, 'p75'], col='gray', border='white', xlab='# Strains',
         ylab='# genes', main='75th percentile', breaks=seq(0, length(strains), 1),
         xaxs='i', yaxs='i')
    abline(v=quantile(osets[targeted, 'n_strains'], 0.75))
    par(op)


    plot(quantile(osets[targeted, 'n_strains'], seq(0.1, 0.9, 0.1)),
         pch=19, cex=1.5, col='black', ylim=c(0, length(strains)+1), bty='n',
         xaxs='i', yaxs='i', xaxt='n', type='b', xlim=c(0, 10),
         ylab='# Strains', xlab='decile',
         main='# strains / orthoset')
    boxplot(deciles[, 1:9], add=TRUE, xaxt='n', frame=FALSE, yaxt='n', col='gray',
            border='gray60')
    legend('topleft', legend=c('Targets', 'Random'), fill=c('black', 'gray'))
    axis(side=1, at=1:9)


    par(op)
    dev.off()

    write.table(rbind(quantile(osets[targeted, 'n_strains'], seq(0.1, 1, 0.1)),
                      deciles),
                file=paste0(outpre, '.strain_count_comparison.deciles.tsv'),
                sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
}


if(opts[['options']][['compare-copy-number']]) {
    ## Compare copy number in targeted genes to copy number in other genes.
    pdf(paste0(outpre, '.copy_number_comparison.pdf'), width=8, height=8, useDingbats=FALSE)
    op = par(no.readonly=TRUE)

    mean_copy_number = osets[, 'n_genes'] / osets[, 'n_strains']
    M = ceiling(max(mean_copy_number))

    ## Distribution of number of copies per strain by orthoset
    par(mfrow=c(2, 1))
    hist(mean_copy_number,
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='All genes',
         ylab='# genes',
         xlab='', breaks=seq(1, M, 0.5),
         mar=c(2.1, 4.1, 3.1, 2.1))
    s = summary(mean_copy_number)
    text(max(mean_copy_number) / 2, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    hist(mean_copy_number[targeted],
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='Targets',
         ylab='# genes',
         xlab='genes / strain',
         breaks=seq(1, M, 0.5),
         mar=c(4.1, 4.1, 3.1, 2.1))
    s = summary(mean_copy_number[targeted])
    text(max(mean_copy_number)/2, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))

    par(op)


    ## Compare target genes to random draws
    ## Match the number of strains in the randomly chosen orthosets to the number
    ## of strains in the target orthosets. Do this in a way that avoids sampling
    ## the same orthoset more than once per random draw.
    par(mfrow=c(4, 2))
    ngenes = sum(targeted)
    n_strains_draw = table(osets[targeted, 'n_strains'])
    rands = matrix(nrow=N, ncol=4)
    colnames(rands) = c('p25', 'median', 'mean', 'p75')
    deciles = matrix(nrow=N, ncol=10)
    for(i in 1:N) {
        if(is.null(background)) {
            r = numeric(ngenes)
            k = 1
            for(j in seq_along(n_strains_draw)) {
                s = sample(mean_copy_number[!targeted & osets[, 'n_strains'] == as.numeric(names(n_strains_draw)[j])],
                           n_strains_draw[j],
                           FALSE)
                r[k:(k+length(s)-1)] = s
                k = k + length(s)
            }
        } else {
            r = mean_copy_number[osets[, gene_column] %in% background[, i]]
        }
        rands[i, 'p25'] = quantile(r, 0.25)
        rands[i, 'p75'] = quantile(r, 0.75)
        rands[i, 'mean'] = mean(r)
        rands[i, 'median'] = median(r)
        deciles[i, ] = quantile(r, seq(0.1, 1, 0.1))
    }

    par(mfcol=c(4, 2))
    hist(rands[, 'p25'], col='gray', border='white', xlab='',
         ylab='# genes', main='25th percentile', breaks=seq(0, M, 0.5),
         xaxs='i', yaxs='i')
    abline(v=quantile(mean_copy_number[targeted], 0.25))
    hist(rands[, 'median'], col='gray', border='white', xlab='',
         ylab='# genes', main='median', breaks=seq(0, M, 0.5),
         xaxs='i', yaxs='i')
    abline(v=median(mean_copy_number[targeted]))
    hist(rands[, 'mean'], col='gray', border='white', xlab='',
         ylab='# genes', main='mean', breaks=seq(0, M, 0.5),
         xaxs='i', yaxs='i')
    abline(v=mean(mean_copy_number[targeted]))
    hist(rands[, 'p75'], col='gray', border='white', xlab='mean copy number',
         ylab='# genes', main='75th percentile', breaks=seq(0, M, 0.5),
         xaxs='i', yaxs='i')
    abline(v=quantile(mean_copy_number[targeted], 0.75))
    par(op)


    plot(quantile(mean_copy_number[targeted], seq(0.1, 0.9, 0.1)),
         pch=19, cex=1.5, col='black', ylim=c(0.9, max(mean_copy_number[targeted], max(deciles))), bty='n',
         xaxs='i', yaxs='i', xaxt='n', type='b', xlim=c(0, 10),
         ylab='mean copy number', xlab='decile',
         main='copy number', xpd=NA)
    boxplot(deciles[, 1:9], add=TRUE, xaxt='n', frame=FALSE, yaxt='n', col='gray',
            border='gray60', xpd=NA)
    legend('topleft', legend=c('Targets', 'Random'), fill=c('black', 'gray'))
    axis(side=1, at=1:9)
    par(op)


    par(mfrow=c(2, 2))
    plot(mean_copy_number[!targeted],
         osets[!targeted, 'n_strains'],
         xlab='mean copy number',
         ylab='number of strains',
         col='gray',
         yaxs='i', xaxs='i', bty='n', xpd=NA,
         main='Background')
    plot(mean_copy_number[targeted],
         osets[targeted, 'n_strains'],
         xlab='mean copy number',
         ylab='number of strains',
         yaxs='i', xaxs='i', bty='n', xpd=NA,
         main='Targets')
    par(op)

    dev.off()

    write.table(rbind(quantile(mean_copy_number[targeted], seq(0.1, 1, 0.1)),
                      deciles),
                file=paste0(outpre, '.copy_number_comparison.deciles.tsv'),
                sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
}



if(opts[['options']][['compare-protein-divergence']]) {
    ## Compare summary of pairwise protein divergence in targets to
    ## background.
    ## Remove NAs before calculating statistics (these are genes with
    ## only one sequence). Also take into account number of copies from
    ## each strain when calculating expected protein divergence.
    pdf(paste0(outpre, '.protein_divergence_comparison.pdf'), width=8, height=8, useDingbats=FALSE)
    op = par(no.readonly=TRUE)

    median_prot_dist = osets[, 'median']
    M = ceiling(max(median_prot_dist, na.rm=TRUE))

    ## Distribution of median pairwise protein distance
    ## Note that NA's are orthosets with only one gene
    par(mfrow=c(2, 1))
    hist(median_prot_dist,
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='All genes',
         ylab='# genes',
         xlab='', breaks=seq(0, M, 0.5),
         mar=c(2.1, 4.1, 3.1, 2.1))
    s = summary(median_prot_dist)
    text(M/2, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    hist(median_prot_dist[targeted],
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='Target genes',
         ylab='# genes',
         xlab='median pairwise protein divergence (subs/site)',
         breaks=seq(0, M, 0.5),
         mar=c(4.1, 4.1, 3.1, 2.1))
    s = summary(median_prot_dist[targeted])
    text(M/2, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    par(op)


    ## Compare target genes to random draws.
    ## Match within a certain amount on number of sequences and
    ## summaries of expected divergence, which are related to
    ## the phylogenetic distribution of sequences
    par(mfrow=c(4, 2))
    ngenes = sum(targeted)
    n_genes_draw = osets[targeted, 'n_genes']
    med_exp_draw = osets[targeted, 'expected_median_copynumber']
    min_exp_draw = osets[targeted, 'expected_min_copynumber']
    max_exp_draw = osets[targeted, 'expected_max_copynumber']
    rands = matrix(nrow=N, ncol=4)
    colnames(rands) = c('p25', 'median', 'mean', 'p75')
    deciles = matrix(nrow=N, ncol=10)
    for(i in 1:N) {
        if(is.null(background)) {
            r = numeric(ngenes)
            for(j in seq_along(n_genes_draw)) {
                ng = n_genes_draw[j]
                de = med_exp_draw[j]
                ne = min_exp_draw[j]
                xe = max_exp_draw[j]
                choices =
                !targeted &
                osets[, 'expected_median_copynumber'] <= de * (1+prot_div_tol) &
                osets[, 'expected_median_copynumber'] >= de * (1-prot_div_tol) &
                osets[, 'expected_min_copynumber'] <= ne * (1+prot_div_tol) &
                osets[, 'expected_min_copynumber'] >= ne * (1-prot_div_tol) &
                osets[, 'expected_max_copynumber'] <= xe * (1+prot_div_tol) &
                osets[, 'expected_max_copynumber'] >= xe * (1-prot_div_tol) &
                osets[, 'n_genes'] <= ng * (1+prot_div_tol) &
                osets[, 'n_genes'] >= ng * (1-prot_div_tol)
                r[j] = sample(median_prot_dist[choices], 1)
            }
        } else {
            r = median_prot_dist[osets[, gene_column] %in% background[, i]]
        }
        rands[i, 'p25'] = quantile(r, 0.25, na.rm=TRUE)
        rands[i, 'p75'] = quantile(r, 0.75, na.rm=TRUE)
        rands[i, 'mean'] = mean(r, na.rm=TRUE)
        rands[i, 'median'] = median(r, na.rm=TRUE)
        deciles[i, ] = quantile(r, seq(0.1, 1, 0.1), na.rm=TRUE)
    }

    par(mfcol=c(4, 2))
    hist(rands[, 'p25'], col='gray', border='white', xlab='',
         ylab='# genes', main='25th percentile', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=quantile(median_prot_dist[targeted], 0.25, na.rm=TRUE))
    hist(rands[, 'median'], col='gray', border='white', xlab='',
         ylab='# genes', main='median', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=median(median_prot_dist[targeted], na.rm=TRUE))
    hist(rands[, 'mean'], col='gray', border='white', xlab='',
         ylab='# genes', main='mean', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=mean(median_prot_dist[targeted], na.rm=TRUE))
    hist(rands[, 'p75'], col='gray', border='white', xlab='median pairwise protein distance',
         ylab='# genes', main='75th percentile', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=quantile(median_prot_dist[targeted], 0.75, na.rm=TRUE))
    par(op)


    plot(quantile(median_prot_dist[targeted], seq(0.1, 0.9, 0.1), na.rm=TRUE),
         pch=19, cex=1.5, col='black', ylim=c(0, max(deciles[, 9])*1.1), bty='n',
         xaxs='i', yaxs='i', xaxt='n', type='b', xlim=c(0, 10),
         ylab='median pairwise protein distance', xlab='decile',
         main='divergence', xpd=NA)
    boxplot(deciles[, 1:9], add=TRUE, xaxt='n', frame=FALSE, yaxt='n', col='gray',
            border='gray60', xpd=NA)
    legend('topleft', legend=c('Targets', 'Random'), fill=c('black', 'gray'))
    axis(side=1, at=1:9)
    par(op)

    dev.off()

    write.table(rbind(quantile(median_prot_dist[targeted], seq(0.1, 1, 0.1), na.rm=TRUE),
                      deciles),
                file=paste0(outpre, '.protein_divergence_comparison.deciles.tsv'),
                sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
}


if(opts[['options']][['compare-relative-protein-divergence']]) {
    ## Compare summary of relative pairwise protein divergence in
    ## targets to background. Note that the relative divergence
    ## is pre-calculated.
    ## Remove NAs before calculating statistics (these are genes with
    ## only one sequence or one strain, paralogs were not included in
    ## the original calculation). Also take into account number of
    ## copies from each strain when calculating expected protein
    ## divergence for random sampling.
    pdf(paste0(outpre, '.relative_protein_divergence_comparison.pdf'), width=8, height=8, useDingbats=FALSE)
    op = par(no.readonly=TRUE)

    median_rel_prot_dist = osets[, 'median_scc']
    lmedian_rel_prot_dist = log2(median_rel_prot_dist)
    lmedian_rel_prot_dist[is.infinite(lmedian_rel_prot_dist)] =
    min(lmedian_rel_prot_dist[!is.infinite(lmedian_rel_prot_dist)], na.rm=TRUE) * 1.1
    M = ceiling(max(median_rel_prot_dist, na.rm=TRUE))
    n = floor(min(median_rel_prot_dist, na.rm=TRUE))
    lM = ceiling(max(lmedian_rel_prot_dist, na.rm=TRUE))
    ln = floor(min(lmedian_rel_prot_dist, na.rm=TRUE))

    ## Distribution of median pairwise protein distance
    ## Note that NA's are orthosets with only one gene or one strain
    par(mfrow=c(2, 1))
    hist(median_rel_prot_dist,
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='All genes',
         ylab='# genes',
         xlab='', breaks=seq(0, M, 5),
         mar=c(2.1, 4.1, 3.1, 2.1))
    s = summary(median_rel_prot_dist)
    text(M/2, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    hist(median_rel_prot_dist[targeted],
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='Target genes',
         ylab='# genes',
         xlab='median relative pairwise protein divergence',
         breaks=seq(0, M, 5),
         mar=c(4.1, 4.1, 3.1, 2.1))
    s = summary(median_rel_prot_dist[targeted])
    text(M/2, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    par(op)

    par(mfrow=c(2, 1))
    hist(lmedian_rel_prot_dist,
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='All genes',
         ylab='# genes',
         xlab='', breaks=seq(ln, lM, 0.5),
         mar=c(2.1, 4.1, 3.1, 2.1))
    s = summary(lmedian_rel_prot_dist)
    text(-5, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    hist(lmedian_rel_prot_dist[targeted],
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='Target genes',
         ylab='# genes',
         xlab='log2(median relative pairwise protein divergence)',
         breaks=seq(ln, lM, 0.5),
         mar=c(4.1, 4.1, 3.1, 2.1))
    s = summary(lmedian_rel_prot_dist[targeted])
    text(-5, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    par(op)

    ## Compare target genes to random draws.
    ## The matching method here does not try to avoid sampling the same background
    ## gene more than once per draw.
    par(mfrow=c(4, 2))
    ngenes = sum(targeted)
    n_genes_draw = osets[targeted, 'n_genes']
    med_exp_draw = osets[targeted, 'expected_median_copynumber']
    min_exp_draw = osets[targeted, 'expected_min_copynumber']
    max_exp_draw = osets[targeted, 'expected_max_copynumber']
    rands = matrix(nrow=N, ncol=4)
    colnames(rands) = c('p25', 'median', 'mean', 'p75')
    lrands = rands
    deciles = matrix(nrow=N, ncol=10)
    ldeciles = deciles
    for(i in 1:N) {
        if(is.null(background)) {
            r = numeric(ngenes)
            lr = numeric(ngenes)
            for(j in seq_along(n_genes_draw)) {
                ng = n_genes_draw[j]
                de = med_exp_draw[j]
                ne = min_exp_draw[j]
                xe = max_exp_draw[j]
                choices =
                !targeted &
                osets[, 'expected_median_copynumber'] <= de * (1+prot_div_tol) &
                osets[, 'expected_median_copynumber'] >= de * (1-prot_div_tol) &
                osets[, 'expected_min_copynumber'] <= ne * (1+prot_div_tol) &
                osets[, 'expected_min_copynumber'] >= ne * (1-prot_div_tol) &
                osets[, 'expected_max_copynumber'] <= xe * (1+prot_div_tol) &
                osets[, 'expected_max_copynumber'] >= xe * (1-prot_div_tol) &
                osets[, 'n_genes'] <= ng * (1+prot_div_tol) &
                osets[, 'n_genes'] >= ng * (1-prot_div_tol)
                r[j] = sample(median_rel_prot_dist[choices], 1)
                lr[j] = sample(lmedian_rel_prot_dist[choices], 1)
            }
        } else {
            r = median_rel_prot_dist[osets[, gene_column] %in% background[, i]]
            lr = lmedian_rel_prot_dist[osets[, gene_column] %in% background[, i]]
        }
        rands[i, 'p25'] = quantile(r, 0.25, na.rm=TRUE)
        rands[i, 'p75'] = quantile(r, 0.75, na.rm=TRUE)
        rands[i, 'mean'] = mean(r)
        rands[i, 'median'] = median(r, na.rm=TRUE)
        deciles[i, ] = quantile(r, seq(0.1, 1, 0.1), na.rm=TRUE)
        lrands[i, 'p25'] = quantile(lr, 0.25, na.rm=TRUE)
        lrands[i, 'p75'] = quantile(lr, 0.75, na.rm=TRUE)
        lrands[i, 'mean'] = mean(lr, na.rm=TRUE)
        lrands[i, 'median'] = median(lr, na.rm=TRUE)
        ldeciles[i, ] = quantile(lr, seq(0.1, 1, 0.1), na.rm=TRUE)
    }

    par(mfcol=c(4, 2))
    hist(rands[, 'p25'], col='gray', border='white', xlab='',
         ylab='# genes', main='25th percentile', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=quantile(median_rel_prot_dist[targeted], 0.25, na.rm=TRUE))
    hist(rands[, 'median'], col='gray', border='white', xlab='',
         ylab='# genes', main='median', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=median(median_rel_prot_dist[targeted], na.rm=TRUE))
    hist(rands[, 'mean'], col='gray', border='white', xlab='',
         ylab='# genes', main='mean', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=mean(median_rel_prot_dist[targeted], na.rm=TRUE))
    hist(rands[, 'p75'], col='gray', border='white', xlab='median pairwise protein distance',
         ylab='# genes', main='75th percentile', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=quantile(median_rel_prot_dist[targeted], 0.75, na.rm=TRUE))
    hist(lrands[, 'p25'], col='gray', border='white', xlab='',
         ylab='# genes', main='25th percentile', breaks=seq(ln, lM, 0.1),
         xaxs='i', yaxs='i')
    abline(v=quantile(lmedian_rel_prot_dist[targeted], 0.25, na.rm=TRUE))
    hist(lrands[, 'median'], col='gray', border='white', xlab='',
         ylab='# genes', main='median', breaks=seq(ln, lM, 0.1),
         xaxs='i', yaxs='i')
    abline(v=median(lmedian_rel_prot_dist[targeted], na.rm=TRUE))
    hist(lrands[, 'mean'], col='gray', border='white', xlab='',
         ylab='# genes', main='mean', breaks=seq(ln, lM, 0.1),
         xaxs='i', yaxs='i')
    abline(v=mean(lmedian_rel_prot_dist[targeted], na.rm=TRUE))
    hist(lrands[, 'p75'], col='gray', border='white', xlab='median pairwise protein distance',
         ylab='# genes', main='75th percentile', breaks=seq(ln, lM, 0.1),
         xaxs='i', yaxs='i')
    abline(v=quantile(lmedian_rel_prot_dist[targeted], 0.75, na.rm=TRUE))
    par(op)


    plot(quantile(median_rel_prot_dist[targeted], seq(0.1, 0.9, 0.1), na.rm=TRUE),
         pch=19, cex=1.5, col='black', ylim=c(0, max(deciles[, 9])*1.1), bty='n',
         xaxs='i', yaxs='i', xaxt='n', type='b', xlim=c(0, 10),
         ylab='median relative pairwise protein distance', xlab='decile',
         main='relative divergence', xpd=NA)
    boxplot(deciles[, 1:9], add=TRUE, xaxt='n', frame=FALSE, yaxt='n', col='gray',
            border='gray60', xpd=NA)
    legend('topleft', legend=c('Targets', 'Random'), fill=c('black', 'gray'))
    axis(side=1, at=1:9)
    par(op)


    plot(quantile(lmedian_rel_prot_dist[targeted], seq(0.1, 0.9, 0.1), na.rm=TRUE),
         pch=19, cex=1.5, col='black', ylim=c(ln, lM), bty='n',
         xaxs='i', yaxs='i', xaxt='n', type='b', xlim=c(0, 10),
         ylab='log2(median relative pairwise protein distance)', xlab='decile',
         main='relative divergence', xpd=NA)
    boxplot(ldeciles[, 1:9], add=TRUE, xaxt='n', frame=FALSE, yaxt='n', col='gray',
            border='gray60', xpd=NA)
    legend('topleft', legend=c('Targets', 'Random'), fill=c('black', 'gray'))
    axis(side=1, at=1:9)
    par(op)

    dev.off()

    write.table(rbind(quantile(median_rel_prot_dist[targeted], seq(0.1, 1, 0.1), na.rm=TRUE),
                      deciles),
                file=paste0(outpre, '.relative_protein_divergence_comparison.deciles.tsv'),
                sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
}


if(opts[['options']][['compare-paralog-protein-divergence']]) {
    ## Compare pairwise divergence among paralogs within target genes
    ## to the background.
    pdf(paste0(outpre, '.paralog_protein_divergence_comparison.pdf'), width=8, height=8, useDingbats=FALSE)
    op = par(no.readonly=TRUE)

    op = par(no.readonly=TRUE)

    median_para_prot_dist = osets[, 'median_dups']
    copy_number = osets[, 'n_genes'] / osets[, 'n_strains']
    M = ceiling(max(median_para_prot_dist, na.rm=TRUE))

    ## Distribution of median pairwise paralog protein distance
    par(mfrow=c(2, 1))
    hist(median_para_prot_dist,
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='All genes',
         ylab='# genes',
         xlab='', breaks=seq(0, M, 0.5),
         mar=c(2.1, 4.1, 3.1, 2.1))
    s = summary(median_para_prot_dist)
    text(M/2, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    hist(median_para_prot_dist[targeted],
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='Target genes',
         ylab='# genes',
         xlab='median pairwise prot. div. between paralogs',
         breaks=seq(0, M, 0.5),
         mar=c(4.1, 4.1, 3.1, 2.1))
    s = summary(median_para_prot_dist[targeted])
    text(M/2, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    par(op)


    ## Compare target genes to random draws.
    ## The matching method here does not try to avoid sampling the same background
    ## gene more than once per draw.
    ## Match on average copy number and total number of sequences
    ngenes = sum(targeted)
    n_genes_draw = osets[targeted, 'n_genes']
    cn_draw = copy_number[targeted]
    rands = matrix(nrow=N, ncol=4)
    colnames(rands) = c('p25', 'median', 'mean', 'p75')
    deciles = matrix(nrow=N, ncol=10)
    for(i in 1:N) {
        if(is.null(background)) {
            r = numeric(ngenes)
            for(j in seq_along(n_genes_draw)) {
                ng = n_genes_draw[j]
                cn = cn_draw[j]
                choices =
                !targeted &
                copy_number <= cn * (1+para_cn_ng_tol) &
                copy_number >= cn * (1-para_cn_ng_tol) &
                osets[, 'n_genes'] <= ng * (1+para_cn_ng_tol) &
                osets[, 'n_genes'] >= ng * (1-para_cn_ng_tol)
                r[j] = sample(median_para_prot_dist[choices], 1)
            }
        } else {
            r = median_para_prot_dist[osets[, gene_column] %in% background[, i]]
        }
        rands[i, 'p25'] = quantile(r, 0.25, na.rm=TRUE)
        rands[i, 'p75'] = quantile(r, 0.75, na.rm=TRUE)
        rands[i, 'mean'] = mean(r, na.rm=TRUE)
        rands[i, 'median'] = median(r, na.rm=TRUE)
        deciles[i, ] = quantile(r, seq(0.1, 1, 0.1), na.rm=TRUE)
    }

    par(mfcol=c(4, 2))
    hist(rands[, 'p25'], col='gray', border='white', xlab='',
         ylab='# genes', main='25th percentile', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=quantile(median_para_prot_dist[targeted], 0.25, na.rm=TRUE))
    hist(rands[, 'median'], col='gray', border='white', xlab='',
         ylab='# genes', main='median', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=median(median_para_prot_dist[targeted], na.rm=TRUE))
    hist(rands[, 'mean'], col='gray', border='white', xlab='',
         ylab='# genes', main='mean', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=mean(median_para_prot_dist[targeted], na.rm=TRUE))
    hist(rands[, 'p75'], col='gray', border='white', xlab='median pairwise paralog protein distance',
         ylab='# genes', main='75th percentile', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=quantile(median_para_prot_dist[targeted], 0.75, na.rm=TRUE))
    par(op)


    plot(quantile(median_para_prot_dist[targeted], seq(0.1, 0.9, 0.1), na.rm=TRUE),
         pch=19, cex=1.5, col='black', ylim=c(0, max(deciles[, 9])*1.1), bty='n',
         xaxs='i', yaxs='i', xaxt='n', type='b', xlim=c(0, 10),
         ylab='median pairwise protein distance b/n duplications', xlab='decile',
         main='divergence between paralogs', xpd=NA)
    boxplot(deciles[, 1:9], add=TRUE, xaxt='n', frame=FALSE, yaxt='n', col='gray',
            border='gray60', xpd=NA)
    legend('topleft', legend=c('Targets', 'Random'), fill=c('black', 'gray'))
    axis(side=1, at=1:9)
    par(op)


    par(mfrow=c(2, 2))
    plot(median_prot_dist[!targeted],
         median_para_prot_dist[!targeted],
         xlim=c(0, M), ylim=c(0, M), col='gray',
         xlab='protein divergence', ylab='paralog protein divergence',
         main='Background genes')
    abline(a=0, b=1, lty=2)
    plot(median_prot_dist[targeted],
         median_para_prot_dist[targeted],
         xlab='protein divergence', ylab='paralog protein divergence',
         xlim=c(0, M), ylim=c(0, M), col='black',
         main='Target genes')
    abline(a=0, b=1, lty=2)
    par(op)

    dev.off()

    write.table(rbind(quantile(median_para_prot_dist[targeted], seq(0.1, 1, 0.1), na.rm=TRUE),
                      deciles),
                file=paste0(outpre, '.paralog_protein_divergence_comparison.deciles.tsv'),
                sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
}


if(opts[['options']][['compare-kaks']]) {
    ## Compare pairwise Ka/Ks values of target genes to the background.
    ## Pairs of sequences within a gene that are too diverged to get a Ks
    ## or Ka value are just ignored.
    pdf(paste0(outpre, '.kaks_comparison.pdf'), width=8, height=8, useDingbats=FALSE)
    op = par(no.readonly=TRUE)

    median_kaks = osets[, 'median_kaks']
    lmedian_kaks = log2(median_kaks)
    lmedian_kaks[is.infinite(lmedian_kaks)] =
    min(lmedian_kaks[!is.infinite(lmedian_kaks)], na.rm=TRUE) * 1.1
    M = ceiling(max(median_kaks, na.rm=TRUE))
    n = floor(min(median_kaks, na.rm=TRUE))
    lM = ceiling(max(lmedian_kaks, na.rm=TRUE))
    ln = floor(min(lmedian_kaks, na.rm=TRUE))

    ## Distribution of median pairwise protein distance
    ## Note that NA's are orthosets with only one gene or one strain
    par(mfrow=c(2, 1))
    hist(median_kaks,
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='All genes',
         ylab='# genes',
         xlab='', breaks=seq(0, M, 0.05),
         mar=c(2.1, 4.1, 3.1, 2.1))
    s = summary(median_kaks)
    text(M/2, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    hist(median_kaks[targeted],
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='Target genes',
         ylab='# genes',
         xlab='median Ka/Ks',
         breaks=seq(0, M, 0.05),
         mar=c(4.1, 4.1, 3.1, 2.1))
    s = summary(median_kaks[targeted])
    text(M/2, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    par(op)

    par(mfrow=c(2, 1))
    hist(lmedian_kaks,
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='All genes',
         ylab='# genes',
         xlab='', breaks=seq(ln, lM, 0.5),
         mar=c(2.1, 4.1, 3.1, 2.1))
    s = summary(lmedian_kaks)
    text(-5, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    hist(lmedian_kaks[targeted],
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='Target genes',
         ylab='# genes',
         xlab='log2(median Ka/Ks)',
         breaks=seq(ln, lM, 0.5),
         mar=c(4.1, 4.1, 3.1, 2.1))
    s = summary(lmedian_kaks[targeted])
    text(-5, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    par(op)

    ## Compare target genes to random draws.
    ## The matching method here does not try to avoid sampling the same background
    ## gene more than once per draw.
    par(mfrow=c(4, 2))
    tol = 0.25
    ngenes = sum(targeted)
    n_genes_draw = osets[targeted, 'n_genes']
    med_exp_draw = osets[targeted, 'expected_median_copynumber']
    min_exp_draw = osets[targeted, 'expected_min_copynumber']
    max_exp_draw = osets[targeted, 'expected_max_copynumber']
    rands = matrix(nrow=N, ncol=4)
    colnames(rands) = c('p25', 'median', 'mean', 'p75')
    lrands = rands
    deciles = matrix(nrow=N, ncol=10)
    ldeciles = deciles
    for(i in 1:N) {
        if(is.null(background)) {
            r = numeric(ngenes)
            lr = numeric(ngenes)
            for(j in seq_along(n_genes_draw)) {
                ng = n_genes_draw[j]
                de = med_exp_draw[j]
                ne = min_exp_draw[j]
                xe = max_exp_draw[j]
                choices =
                !targeted &
                osets[, 'expected_median_copynumber'] <= de * (1+prot_div_tol) &
                osets[, 'expected_median_copynumber'] >= de * (1-prot_div_tol) &
                osets[, 'expected_min_copynumber'] <= ne * (1+prot_div_tol) &
                osets[, 'expected_min_copynumber'] >= ne * (1-prot_div_tol) &
                osets[, 'expected_max_copynumber'] <= xe * (1+prot_div_tol) &
                osets[, 'expected_max_copynumber'] >= xe * (1-prot_div_tol) &
                osets[, 'n_genes'] <= ng * (1+prot_div_tol) &
                osets[, 'n_genes'] >= ng * (1-prot_div_tol)
                r[j] = sample(median_kaks[choices], 1)
                lr[j] = sample(lmedian_kaks[choices], 1)
            }
        } else {
            r = median_kaks[osets[, gene_column] %in% background[, i]]
            lr = lmedian_kaks[osets[, gene_column] %in% background[, i]]
        }
        rands[i, 'p25'] = quantile(r, 0.25, na.rm=TRUE)
        rands[i, 'p75'] = quantile(r, 0.75, na.rm=TRUE)
        rands[i, 'mean'] = mean(r, na.rm=TRUE)
        rands[i, 'median'] = median(r, na.rm=TRUE)
        deciles[i, ] = quantile(r, seq(0.1, 1, 0.1), na.rm=TRUE)
        lrands[i, 'p25'] = quantile(lr, 0.25, na.rm=TRUE)
        lrands[i, 'p75'] = quantile(lr, 0.75, na.rm=TRUE)
        lrands[i, 'mean'] = mean(lr, na.rm=TRUE)
        lrands[i, 'median'] = median(lr, na.rm=TRUE)
        ldeciles[i, ] = quantile(lr, seq(0.1, 1, 0.1), na.rm=TRUE)
    }

    par(mfcol=c(4, 2))
    hist(rands[, 'p25'], col='gray', border='white', xlab='',
         ylab='# genes', main='25th percentile', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=quantile(median_kaks[targeted], 0.25, na.rm=TRUE))
    hist(rands[, 'median'], col='gray', border='white', xlab='',
         ylab='# genes', main='median', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=median(median_kaks[targeted], na.rm=TRUE))
    hist(rands[, 'mean'], col='gray', border='white', xlab='',
         ylab='# genes', main='mean', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=mean(median_kaks[targeted], na.rm=TRUE))
    hist(rands[, 'p75'], col='gray', border='white', xlab='median Ka/Ks',
         ylab='# genes', main='75th percentile', breaks=seq(0, M, 0.1),
         xaxs='i', yaxs='i')
    abline(v=quantile(median_kaks[targeted], 0.75, na.rm=TRUE))
    hist(lrands[, 'p25'], col='gray', border='white', xlab='',
         ylab='# genes', main='25th percentile', breaks=seq(ln, lM, 0.1),
         xaxs='i', yaxs='i')
    abline(v=quantile(lmedian_kaks[targeted], 0.25, na.rm=TRUE))
    hist(lrands[, 'median'], col='gray', border='white', xlab='',
         ylab='# genes', main='median', breaks=seq(ln, lM, 0.1),
         xaxs='i', yaxs='i')
    abline(v=median(lmedian_kaks[targeted], na.rm=TRUE))
    hist(lrands[, 'mean'], col='gray', border='white', xlab='',
         ylab='# genes', main='mean', breaks=seq(ln, lM, 0.1),
         xaxs='i', yaxs='i')
    abline(v=mean(lmedian_kaks[targeted], na.rm=TRUE))
    hist(lrands[, 'p75'], col='gray', border='white', xlab='median Ka/Ks',
         ylab='# genes', main='75th percentile', breaks=seq(ln, lM, 0.1),
         xaxs='i', yaxs='i')
    abline(v=quantile(lmedian_kaks[targeted], 0.75, na.rm=TRUE))
    par(op)


    plot(quantile(median_kaks[targeted], seq(0.1, 0.9, 0.1), na.rm=TRUE),
         pch=19, cex=1.5, col='black', ylim=c(0, max(deciles[, 9])*1.1), bty='n',
         xaxs='i', yaxs='i', xaxt='n', type='b', xlim=c(0, 10),
         ylab='median Ka/Ks', xlab='decile',
         main='KaKs', xpd=NA)
    boxplot(deciles[, 1:9], add=TRUE, xaxt='n', frame=FALSE, yaxt='n', col='gray',
            border='gray60', xpd=NA)
    legend('topleft', legend=c('Targets', 'Random'), fill=c('black', 'gray'))
    axis(side=1, at=1:9)
    par(op)


    plot(quantile(lmedian_kaks[targeted], seq(0.1, 0.9, 0.1), na.rm=TRUE),
         pch=19, cex=1.5, col='black', ylim=c(ln, lM), bty='n',
         xaxs='i', yaxs='i', xaxt='n', type='b', xlim=c(0, 10),
         ylab='log2(median kaks)', xlab='decile',
         main='KaKs', xpd=NA)
    boxplot(ldeciles[, 1:9], add=TRUE, xaxt='n', frame=FALSE, yaxt='n', col='gray',
            border='gray60', xpd=NA)
    legend('topleft', legend=c('Targets', 'Random'), fill=c('black', 'gray'))
    axis(side=1, at=1:9)
    par(op)

    dev.off()

    write.table(rbind(quantile(median_kaks[targeted], seq(0.1, 1, 0.1), na.rm=TRUE),
                      deciles),
                file=paste0(outpre, '.kaks_comparison.deciles.tsv'),
                sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
}



if(opts[['options']][['compare-paralog-kaks']]) {
    ## Compare pairwise Ka/Ks values of paralogs in target genes to the background.
    ## Pairs of sequences within a gene that are too diverged to get a Ks
    ## or Ka value are just ignored.
    pdf(paste0(outpre, '.kaks_paralog_comparison.pdf'), width=8, height=8, useDingbats=FALSE)
    op = par(no.readonly=TRUE)

     median_para_kaks = osets[, 'median_dups_kaks']
     M = ceiling(max(median_para_kaks, na.rm=TRUE))
     n = floor(min(median_para_kaks, na.rm=TRUE))
     lmedian_para_kaks = log2(median_para_kaks)
     lmedian_para_kaks[is.infinite(lmedian_para_kaks)] =
     min(lmedian_para_kaks[!is.infinite(lmedian_para_kaks)], na.rm=TRUE) * 1.1
     lM = ceiling(max(lmedian_para_kaks, na.rm=TRUE))
     ln = floor(min(lmedian_para_kaks, na.rm=TRUE))

     ## Distribution of median pairwise Ka/Ks
     par(mfrow=c(2, 1))
     hist(median_para_kaks,
          col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
          cex.axis=1.5, main='All genes',
          ylab='# genes',
          xlab='', breaks=seq(0, M, 0.1),
          mar=c(2.1, 4.1, 3.1, 2.1))
     s = summary(median_para_kaks)
     text(M/2, par('usr')[4]*0.6,
          paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
     hist(median_para_kaks[targeted],
          col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
          cex.axis=1.5, main='Target genes',
          ylab='# genes',
          xlab='median pairwise prot. div. between paralogs',
          breaks=seq(0, M, 0.1),
          mar=c(4.1, 4.1, 3.1, 2.1))
     s = summary(median_para_kaks[targeted])
     text(M/2, par('usr')[4]*0.6,
          paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
     par(op)


     ## Compare target genes to random draws.
     ## The matching method here does not try to avoid sampling the same background
     ## gene more than once per draw.
     ## Match on copy number and total number of sequences
     ngenes = sum(targeted)
     n_genes_draw = osets[targeted, 'n_genes']
     cn_draw = copy_number[targeted]
     rands = matrix(nrow=N, ncol=4)
     colnames(rands) = c('p25', 'median', 'mean', 'p75')
     deciles = matrix(nrow=N, ncol=10)
     for(i in 1:N) {
         if(is.null(background)) {
             r = numeric(ngenes)
             for(j in seq_along(n_genes_draw)) {
                 ng = n_genes_draw[j]
                 cn = cn_draw[j]
                 choices =
                 !targeted &
                 copy_number <= cn * (1+para_cn_ng_tol) &
                 copy_number >= cn * (1-para_cn_ng_tol) &
                 osets[, 'n_genes'] <= ng * (1+para_cn_ng_tol) &
                 osets[, 'n_genes'] >= ng * (1-para_cn_ng_tol)
                 r[j] = sample(median_para_kaks[choices], 1)
             }
         } else {
             r = median_para_kaks[osets[, gene_column] %in% background[, i]]
         }
         rands[i, 'p25'] = quantile(r, 0.25, na.rm=TRUE)
         rands[i, 'p75'] = quantile(r, 0.75, na.rm=TRUE)
         rands[i, 'mean'] = mean(r, na.rm=TRUE)
         rands[i, 'median'] = median(r, na.rm=TRUE)
         deciles[i, ] = quantile(r, seq(0.1, 1, 0.1), na.rm=TRUE)
     }

     par(mfcol=c(4, 2))
     hist(rands[, 'p25'], col='gray', border='white', xlab='',
          ylab='# genes', main='25th percentile', breaks=seq(0, M, 0.1),
          xaxs='i', yaxs='i')
     abline(v=quantile(median_para_kaks[targeted], 0.25, na.rm=TRUE))
     hist(rands[, 'median'], col='gray', border='white', xlab='',
          ylab='# genes', main='median', breaks=seq(0, M, 0.1),
          xaxs='i', yaxs='i')
     abline(v=median(median_para_kaks[targeted], na.rm=TRUE))
     hist(rands[, 'mean'], col='gray', border='white', xlab='',
          ylab='# genes', main='mean', breaks=seq(0, M, 0.1),
          xaxs='i', yaxs='i')
     abline(v=mean(median_para_kaks[targeted], na.rm=TRUE))
     hist(rands[, 'p75'], col='gray', border='white', xlab='median pairwise paralog kaks',
          ylab='# genes', main='75th percentile', breaks=seq(0, M, 0.1),
          xaxs='i', yaxs='i')
     abline(v=quantile(median_para_kaks[targeted], 0.75, na.rm=TRUE))
     par(op)


     plot(quantile(median_para_kaks[targeted], seq(0.1, 0.9, 0.1), na.rm=TRUE),
          pch=19, cex=1.5, col='black', ylim=c(0, max(deciles[, 9])*1.1), bty='n',
          xaxs='i', yaxs='i', xaxt='n', type='b', xlim=c(0, 10),
          ylab='median pairwise KaKs b/n duplications', xlab='decile',
          main='Ka/Ks between paralogs', xpd=NA)
     boxplot(deciles[, 1:9], add=TRUE, xaxt='n', frame=FALSE, yaxt='n', col='gray',
             border='gray60', xpd=NA)
     legend('topleft', legend=c('Targets', 'Random'), fill=c('black', 'gray'))
     axis(side=1, at=1:9)
     par(op)


     par(mfrow=c(2, 2))
     plot(median_kaks[!targeted],
          median_para_kaks[!targeted],
          xlim=c(0, M), ylim=c(0, M), col='gray',
          xlab='median all pairs Ka/Ks', ylab='median paralog Ka/Ks',
          main='background genes')
     abline(a=0, b=1, lty=2)
     plot(median_kaks[targeted],
          median_para_kaks[targeted],
          xlab='median all pairs Ka/Ks', ylab='median paralog Ka/Ks',
          xlim=c(0, M), ylim=c(0, M), col='black',
          main='target genes')
     abline(a=0, b=1, lty=2)
     hist(lmedian_kaks[!targeted] -
          lmedian_para_kaks[!targeted],
          col='gray', border='white', breaks=seq(ln-0.1, -ln, 0.5),
          xaxs='i', yaxs='i', xpd=NA, ylab='# Orthosets',
          xlab='log(median all Ka/Ks) -\nlog2(median paralogs Ka/Ks)',
          main='')
     abline(v=0, lty=3, col='black')
     hist(lmedian_kaks[targeted] -
          lmedian_para_kaks[targeted],
          col='gray', border='white', breaks=seq(ln-0.1, -ln, 0.5),
          xaxs='i', yaxs='i', xpd=NA, ylab='# Orthosets',
          xlab='log(median all Ka/Ks) -\nlog2(median paralogs Ka/Ks)',
          main='')
     abline(v=0, lty=3, col='black')
     par(op)

     dev.off()

    write.table(rbind(quantile(median_para_kaks[targeted], seq(0.1, 1, 0.1), na.rm=TRUE),
                      deciles),
                file=paste0(outpre, '.kaks_paralog_comparison.deciles.tsv'),
                sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
}



if(opts[['options']][['compare-expected-divergence']]) {
    ## Compare expected pairwise distances in target genes to the background.
    ## This is a way of looking at the phylogenetic distribution of genes.
    ## Counts each strain once.
    pdf(paste0(outpre, '.expected_divergence_comparison.pdf'), width=8, height=8, useDingbats=FALSE)
    op = par(no.readonly=TRUE)
    op = par(no.readonly=TRUE)

    median_exp_div = osets[, 'expected_median']
    max_exp_div = osets[, 'expected_max']
    Md = ceiling(max(median_exp_div))
    Mx = ceiling(max(max_exp_div))

    ## Distribution of expected divergence
    par(mfrow=c(2, 1))
    hist(median_exp_div[!targeted],
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='Background genes',
         ylab='# orthosets',
         xlab='', breaks=seq(0, Md, 0.05),
         mar=c(2.1, 4.1, 3.1, 2.1))
    s = summary(median_exp_div[!targeted])
    text(Md * 0.7, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    hist(median_exp_div[targeted],
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='Target genes',
         ylab='# orthosets',
         xlab='median expected divergence',
         breaks=seq(0, Md, 0.05),
         mar=c(4.1, 4.1, 3.1, 2.1))
    s = summary(median_exp_div[targeted])
    text(Md * 0.7, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    par(op)

    par(mfrow=c(2, 1))
    hist(max_exp_div[!targeted],
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='Background genes',
         ylab='# orthosets',
         xlab='', breaks=seq(0, Mx, 0.05),
         mar=c(2.1, 4.1, 3.1, 2.1))
    s = summary(max_exp_div[!targeted])
    text(Mx * 0.7, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    hist(max_exp_div[targeted],
         col='gray', border='white', xaxs='i', yaxs='i', cex.lab=1.5,
         cex.axis=1.5, main='Target genes',
         ylab='# orthosets',
         xlab='max expected divergence',
         breaks=seq(0, Mx, 0.05),
         mar=c(4.1, 4.1, 3.1, 2.1))
    s = summary(max_exp_div[targeted])
    text(Mx * 0.7, par('usr')[4]*0.6,
         paste(mapply(function(x, y) paste(x, y, sep=':  '), names(s), round(s, 1)), collapse='\n'))
    par(op)


    ## Compare target genes to random draws
    ## Match the number of strains in the randomly chosen orthosets to the number
    ## of strains in the target orthosets. Do this in a way that avoids sampling
    ## the same orthoset more than once per random draw.
    par(mfrow=c(4, 2))
    ngenes = sum(targeted)
    n_strains_draw = table(osets[targeted, 'n_strains'])
    mrands = matrix(nrow=N, ncol=4)
    colnames(mrands) = c('p25', 'median', 'mean', 'p75')
    mdeciles = matrix(nrow=N, ncol=10)
    xrands = mrands
    xdeciles = mdeciles
    for(i in 1:N) {
        if(is.null(background)) {
            r_med = numeric(ngenes)
            r_max = numeric(ngenes)
            k = 1
            for(j in seq_along(n_strains_draw)) {
                ii = sample(which(!targeted & osets[, 'n_strains'] == as.numeric(names(n_strains_draw)[j])),
                            n_strains_draw[j],
                            FALSE)
                sm = median_exp_div[ii]
                sx = max_exp_div[ii]
                r_med[k:(k+length(ii)-1)] = sm
                r_max[k:(k+length(ii)-1)] = sx
                k = k + length(ii)
            }
        } else {
            r_med = median_exp_div[osets[, gene_column] %in% background[, i]]
            r_max = max_exp_div[osets[, gene_column] %in% background[, i]]
        }
        mrands[i, 'p25'] = quantile(r_med, 0.25)
        mrands[i, 'p75'] = quantile(r_med, 0.75)
        mrands[i, 'mean'] = mean(r_med)
        mrands[i, 'median'] = median(r_med)
        mdeciles[i, ] = quantile(r_med, seq(0.1, 1, 0.1))
        xrands[i, 'p25'] = quantile(r_max, 0.25)
        xrands[i, 'p75'] = quantile(r_max, 0.75)
        xrands[i, 'mean'] = mean(r_max)
        xrands[i, 'median'] = median(r_max)
        xdeciles[i, ] = quantile(r_max, seq(0.1, 1, 0.1))
    }

    par(mfcol=c(4, 2))
    hist(mrands[, 'p25'], col='gray', border='white', xlab='',
         ylab='# genes', main='25th percentile', breaks=seq(0, Md, 0.05),
         xaxs='i', yaxs='i')
    abline(v=quantile(median_exp_div[targeted], 0.25))
    hist(mrands[, 'median'], col='gray', border='white', xlab='',
         ylab='# genes', main='median', breaks=seq(0, Md, 0.05),
         xaxs='i', yaxs='i')
    abline(v=median(median_exp_div[targeted]))
    hist(mrands[, 'mean'], col='gray', border='white', xlab='',
         ylab='# genes', main='mean', breaks=seq(0, Md, 0.05),
         xaxs='i', yaxs='i')
    abline(v=mean(median_exp_div[targeted]))
    hist(mrands[, 'p75'], col='gray', border='white', xlab='median exp. prot. div.',
         ylab='# genes', main='75th percentile', breaks=seq(0, Md, 0.05),
         xaxs='i', yaxs='i')
    abline(v=quantile(median_exp_div[targeted], 0.75))

    hist(xrands[, 'p25'], col='gray', border='white', xlab='',
         ylab='# genes', main='25th percentile', breaks=seq(0, Md, 0.05),
         xaxs='i', yaxs='i')
    abline(v=quantile(max_exp_div[targeted], 0.25))
    hist(xrands[, 'median'], col='gray', border='white', xlab='',
         ylab='# genes', main='median', breaks=seq(0, Md, 0.05),
         xaxs='i', yaxs='i')
    abline(v=median(max_exp_div[targeted]))
    hist(xrands[, 'mean'], col='gray', border='white', xlab='',
         ylab='# genes', main='mean', breaks=seq(0, Md, 0.05),
         xaxs='i', yaxs='i')
    abline(v=mean(max_exp_div[targeted]))
    hist(xrands[, 'p75'], col='gray', border='white', xlab='max exp. prot. div.',
         ylab='# genes', main='75th percentile', breaks=seq(0, Md, 0.05),
         xaxs='i', yaxs='i')
    abline(v=quantile(max_exp_div[targeted], 0.75))
    par(op)


    plot(quantile(median_exp_div[targeted], seq(0.1, 0.9, 0.1)),
         pch=19, cex=1.5, col='black', ylim=c(0, Mx), bty='n',
         xaxs='i', yaxs='i', xaxt='n', type='b', xlim=c(0, 10),
         ylab='median exp. prot. div.', xlab='decile',
         main='expected divergence', xpd=NA)
    boxplot(mdeciles[, 1:9], add=TRUE, xaxt='n', frame=FALSE, yaxt='n', col='gray',
            border='gray60', xpd=NA)
    legend('topleft', legend=c('Targets', 'Random'), fill=c('black', 'gray'))
    axis(side=1, at=1:9)
    par(op)

    plot(quantile(max_exp_div[targeted], seq(0.1, 0.9, 0.1)),
         pch=19, cex=1.5, col='black', ylim=c(0, Mx), bty='n',
         xaxs='i', yaxs='i', xaxt='n', type='b', xlim=c(0, 10),
         ylab='max. exp. prot. div.', xlab='decile',
         main='expected divergence', xpd=NA)
    boxplot(xdeciles[, 1:9], add=TRUE, xaxt='n', frame=FALSE, yaxt='n', col='gray',
            border='gray60', xpd=NA)
    legend('topleft', legend=c('Targets', 'Random'), fill=c('black', 'gray'))
    axis(side=1, at=1:9)
    par(op)


    par(mfrow=c(2, 2))
    plot(median_exp_div[!targeted],
         osets[!targeted, 'n_strains'],
         xlab='median exp. prot. div.',
         ylab='number of strains',
         col='gray',
         yaxs='i', xaxs='i', bty='n', xpd=NA,
         main='Background genes')
    plot(median_exp_div[targeted],
         osets[targeted, 'n_strains'],
         xlab='median exp. prot. div.',
         ylab='number of strains',
         yaxs='i', xaxs='i', bty='n', xpd=NA,
         main='Target genes')
    plot(max_exp_div[!targeted],
         osets[!targeted, 'n_strains'],
         xlab='max. exp. prot. div.',
         ylab='number of strains',
         col='gray',
         yaxs='i', xaxs='i', bty='n', xpd=NA,
         main='Background genes')
    plot(max_exp_div[targeted],
         osets[targeted, 'n_strains'],
         xlab='max. exp. prot. div.',
         ylab='number of strains',
         yaxs='i', xaxs='i', bty='n', xpd=NA,
         main='Target genes')
    par(op)

    dev.off()

    write.table(rbind(quantile(max_exp_div[targeted], seq(0.1, 1, 0.1)),
                      xdeciles),
                file=paste0(outpre, '.expected_divergence_comparison.deciles.max.tsv'),
                sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
    write.table(rbind(quantile(median_exp_div[targeted], seq(0.1, 1, 0.1)),
                      mdeciles),
                file=paste0(outpre, '.expected_divergence_comparison.deciles.median.tsv'),
                sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
}

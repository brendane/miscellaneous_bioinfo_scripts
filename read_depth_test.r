#!/usr/bin/env Rscript
#
# Use a weighted linear regression to test for differences in read depth
# among treatments. The weights come from the total number of aligned
# reads.
#
#   read_depth_test.r <file with group 1 input files>
#       <file with group 2 input files>
#
#==============================================================================#

cargs = commandArgs(trailingOnly=TRUE)
group_1 = scan(cargs[1], what='character')
group_2 = scan(cargs[2], what='character')

trtmnt = c(rep('a', length(group_1)),
           rep('b', length(group_2)))

data_sets = lapply(c(group_1, group_2),
                   function(x) read.csv(x, sep='\t',
                                        stringsAsFactors=FALSE,
                                        check.names=FALSE,
                                        comment.char='#',
                                        header=TRUE))
n_genes = nrow(data_sets[[1]])

same_gene = sapply(1:n_genes,
                   function(i) {
                       genes = sapply(data_sets,
                                       function(x) x[i, 'name'])
                       all(genes == genes[1])
                    })
if(!all(same_gene))
    stop('Gene names do not match')

cat('contig\tstart\tend\tgene\tp\teffect\trpkm1\trpkm2\n', file=stdout())
for(i in 1:n_genes) {
    contig = data_sets[[1]][i, 'chrom']
    start_pos = data_sets[[1]][i, 'start']
    end_pos = data_sets[[1]][i, 'end']
    gene = data_sets[[1]][i, 'name']

    raw_counts = sapply(data_sets, function(x) x[i, 'count'])
    cpm = sapply(data_sets, function(x) x[i, 'count (CPM)'])
    rpkm = sapply(data_sets, function(x) x[i, 'RPKM'])

    if(sum(raw_counts[trtmnt == 'a'] != 0) > 0 &&
       sum(raw_counts[trtmnt == 'b'] != 0) > 0) {
        trts = trtmnt[raw_counts > 0]
        cpm = cpm[raw_counts > 0]
        rpkm = rpkm[raw_counts > 0]
        raw_counts = raw_counts[raw_counts > 0]
        wts = raw_counts / cpm

        test_results = lm(rpkm ~ trts, weights=wts)
        effect = coefficients(test_results)['trtsb']
        p = coefficients(summary(test_results))['trtsb', 'Pr(>|t|)']
    } else {
        p = NaN
        effect = NaN
    }
    cat(paste(contig, start_pos, end_pos, gene, p, effect, sep='\t'),
        '\t', paste(round(rpkm[trts == 'a'], 3), collapse=','),
        '\t', paste(round(rpkm[trts == 'b'], 3), collapse=','),
        '\n', sep='', file=stdout())
}

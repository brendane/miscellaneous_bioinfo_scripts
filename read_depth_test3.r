#!/usr/bin/env Rscript
#
# Use a weighted linear regression to test for differences in read depth
# among treatments. The weights come from the total number of aligned
# reads.
#
#   read_depth_test.r <file with group 1 input files>
#       <file with group 2 input files>
#
# In this version, the square root of the aligned reads is used.
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

# Get the weights - based on the number of aligned reads. Divide the
# count of aligned reads by the count of aligned reads million to get
# the number of aligned reads in millions. This only has to be done
# once per pool.
wts = sapply(data_sets, function(x) {
             r = which(x$count > 10)[1]
             sqrt(x[r, 'count'] / x[r, 'count (CPM)'])
                    })
cat(wts, '\n', file=stderr())

cat('contig\tstart\tend\tgene\tp\teffect\trpkm1\trpkm2\n', file=stdout())
for(i in 1:n_genes) {
    contig = data_sets[[1]][i, 'chrom']
    start_pos = data_sets[[1]][i, 'start']
    end_pos = data_sets[[1]][i, 'end']
    gene = data_sets[[1]][i, 'name']

    raw_counts = sapply(data_sets, function(x) x[i, 'count'])
    cpm = sapply(data_sets, function(x) x[i, 'count (CPM)'])
    rpkm = sapply(data_sets, function(x) x[i, 'RPKM'])

    # A previous version of this script had an "and" instead of an
    # "or", but I think it's okay to do a test if one treatment has
    # no reads, as long as the other one has reads.
    if(sum(raw_counts[trtmnt == 'a'] != 0) > 0 ||
       sum(raw_counts[trtmnt == 'b'] != 0) > 0) {
        # In a previous version of this script, I eliminated replicates
        # with zero reads. I have no idea why I did that - it seems
        # silly.
        test_results = lm(rpkm ~ trtmnt, weights=wts)
        effect = coefficients(test_results)['trtmntb']
        p = coefficients(summary(test_results))['trtmntb', 'Pr(>|t|)']
    } else {
        p = NaN
        effect = NaN
    }
    cat(paste(contig, start_pos, end_pos, gene, p, effect, sep='\t'),
        '\t', paste(round(rpkm[trtmnt == 'a'], 3), collapse=','),
        '\t', paste(round(rpkm[trtmnt == 'b'], 3), collapse=','),
        '\n', sep='', file=stdout())
}

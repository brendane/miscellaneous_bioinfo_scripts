#!/usr/bin/env Rscript
#
# Read a phylip distance matrix and create a neighbor-joining tree from
# the output.

library(ape)

read_phylip_full_dist_mat = function(fname) {
    dist_handle = file(fname, 'rb')
    n_samples = scan(dist_handle, nlines=1, quiet=TRUE)
    d = matrix(nrow=n_samples, ncol=n_samples)
    sample_names = character(n_samples)
    input_matrix = scan(dist_handle, quiet=TRUE, what='character')
    for(i in 1:n_samples) {
        row = input_matrix[((i-1) * (n_samples + 1) + 1):(i * (n_samples + 1))]
        sample_names[i] = row[1]
        d[i, i] = 0
        for(j in seq_along(row[-1])) {
            d[i, j] = as.numeric(row[j+1])
        }
    }   
    rownames(d) = sample_names
    colnames(d) = sample_names
    close(dist_handle)
    d
}

cargs = commandArgs(trailingOnly=TRUE)

dists = read_phylip_full_dist_mat(cargs[1])

tree = nj(as.dist(dists))
write.tree(tree, stdout())

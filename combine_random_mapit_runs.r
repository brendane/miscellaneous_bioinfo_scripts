#!/usr/bin/env Rscript
#
# combine_random_mapit_runs.r <files...>
#

library(data.table)

argv = commandArgs(trailingOnly=TRUE)

results = numeric(length(argv))
for(i in seq_along(argv)) {
    results[i] = sort(fread(argv[i], header=TRUE)[['p']], decreasing=FALSE)[1]
}

cat(sort(results, decreasing=FALSE), sep='\n', file=stdout())

#!/usr/bin/env Rscript

library(igraph)

argv = commandArgs(trailingOnly=TRUE)

## Read adjacency matrix
adjs = as.matrix(read.csv(argv[1], sep='\t', header=TRUE,
                          check.names=FALSE, comment.char=''))

## Make a graph
net = graph.adjacency(1-adjs, mode='undirected', weighted=TRUE, diag=FALSE)

## Edge names
e = attr(E(net), 'vnames')

## Calculate betweenness (takes a long time)
bet = estimate_betweenness(net, cutoff=5, directed=FALSE)
write.table(cbind(rep('vertex', length(bet)), names(bet), bet),
            file=stdout(), quote=FALSE, col.names=c('type', 'gene', 'betweenness'),
            row.names=FALSE)

ebet = estimate_edge_betweenness(net, cutoff=5, directed=FALSE)
write.table(cbind(rep('edge', length(ebet)), e, ebet),
            file=stdout(), quote=FALSE, col.names=FALSE, row.names=FALSE)

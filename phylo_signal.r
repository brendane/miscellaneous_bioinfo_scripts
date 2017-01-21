#!/usr/bin/env Rscript
#
# This function calculates a silly phylogenetic signal score that
# I made up for use on unrooted trees.
#
# For each split in the tree with >= N taxa on each side, it calculates
# the mean trait value on each side of the split, and takes the absolute
# value of the difference. I made this statistic up.
#
# It is important that there not be any duplicated names in the tree or
# in the traits.
#

silly_phy_signal = function(tree, traits, N) {
    #
    # tree should be a phylo object
    # traits should be a two column data.frame - 1st column is
    #   taxa names, 2nd column is trait values
    # N is the minimum number of leaves on each side of the
    #   split to allow for testing

    require(ape)
    require(picante)

    traits = traits[!is.na(as.numeric(traits[, 2])), ]
    traitv = as.numeric(traits[, 2])
    names(traitv) = traits[, 1]

    if(any(duplicated(traits[, 1])))
        stop('Duplicated taxa in traits file')
    if(any(duplicated(tree$tip.label)))
        stop('Duplicated names in tree')

    traitv = traitv[tree$tip.label]
    splits = prop.part(tree)
    ntip = Ntip(tree)
    maxval = 0
    for(s in splits) {
        n1 = sum(!is.na(traitv[s]))
        n2 = sum(!is.na(traitv[-s]))
        if(n1 < N || n2 < N)
            next
        m1 = mean(traitv[s], na.rm=TRUE)
        m2 = mean(traitv[-s], na.rm=TRUE)
        v  = abs(m1 - m2)
        if(v > maxval)
            maxval = v
    }
    maxval
}

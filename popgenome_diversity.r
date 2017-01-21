#!/usr/bin/env Rscript
.doc = '
    Calculate a variety of common diversity statistics on a single
    population using PopGenome.

    popgenome_diversity.r --output <output prefix> --window --step
        <input tabixed VCF> <gff file> <ref fasta file> <replicon>
        <text file with individuals in the population>

   Does not use an outgroup and only handles one population.
'

#==============================================================================#

library(magrittr)
library(optparse)
library(PopGenome)

#==============================================================================#

add_gene_names = function(obj, gff_file, replicon) {
    genes = get.feature.names(obj, gff_file, replicon)
    gene_ids    = gsub('.+ID=(.+?);.+', '\\1', genes)
    gene_names  = gsub('.+gene=(.+?);.+', '\\1', genes)
    gene_names[!grepl('gene=.+;', genes)] = ''
    obj@region.names = paste0(gene_ids, ': ', gene_names)
    obj
}

basic_stats = function(input_data, subsites) {
    output_list = structure(vector('list', length=length(subsites)),
                            names=subsites)
    for(subsites in names(output_list)) {
        sbs = subsites
        if(sbs == 'all') sbs = FALSE
        stats1= input_data %>%
            neutrality.stats(subsites=sbs) %>%
            get.neutrality(., theta=TRUE) %>%
            `[[`(1)
        Pi = input_data %>%
            F_ST.stats(mode='nucleotide') %>%
            `@`('nuc.diversity.within') %>%
            `[[`(1)
        if(!all(rownames(stats1) == rownames(stats2)))
            stop('Gene names do not match')
        output_list[[subsites]] = data.frame(stats1, stats2,
                                             stringsAsFactors=FALSE)
    }
    output_list
}

#==============================================================================#

optlist = list(make_option('--output'),
               make_option('--window', type='double'),
               make_option('--step', type='double'))
opts = parse_args(OptionParser(option_list=optlist, usage=.doc),
                  positional_arguments=TRUE)
outpref   = opts[['options']][['output']]
winsize   = opts[['options']][['window']]
stepsize  = opts[['options']][['step']]
vcffile   = opts[['args']][1]
gfffile   = opts[['args']][2]
reffile   = opts[['args']][3]
replicon  = opts[['args']][4]
popfile   = opts[['args']][5]
replen    = opts[['args']][6] %>% as.numeric()

# List of strains to include
strains = scan(popfile, what='character')

# Read in VCF file and add SNP type annotations; assume haploid
# To do this with a set of FASTA files, just put them all in a
# directory and point readData(path=directory, big.data=TRUE, ...)
# at the directory.
g = readVCF(filename=vcffile,
            gffpath=gfffile,
            frompos=1,
            topos=replen,
            numcols=10000,
            tid=replicon,
            approx=TRUE,
            include.unknown=TRUE) %>%
    set.populations(., new.populations=list(which(get.individuals(.) %in% strains))) %>%
    set.synnonsyn(ref.chr=reffile, save.codons=TRUE)


## Diversity statistics for the entire genome
whole_genome_stats = basic_stats(g, c('all', 'coding', 'syn',
                                      'nonsyn', 'intergenic'))


## Diversity statistics for genes
genic_data = g %>%
    splitting.data(subsites='coding') %>%
    add_gene_names(gfffile, replicon)
genic_stats = basic_stats(genic_data, c('all', 'syn', 'nonsyn'))


## Diversity statistics for sliding windows
window_data = sliding.window.transform(g, width=winsize, jump=stepsize,
                                       type=2, start.pos=FALSE,
                                       end.pos=FALSE, whole.data=FALSE)
window_stats = basic_stats(window_data, c('all', 'genes', 'syn',
                                          'nonsyn', 'intergenic'))

#!/usr/bin/env Rscript
.doc = '
    Calculate a variety of common diversity statistics on a single
    population using PopGenome.

    popgenome_diversity.r --output <output prefix> --window --step
        <input FASTA directory> <gff file> <ref fasta file> <replicon>
        <text file with individuals in the population>
        <nsites by genes> <nsites by 250 bp windows>

   Does not use an outgroup and only handles one population.

   Note that it is important that the fasta file and gff file have
   the same basename.
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
        stats1 = input_data %>%
            neutrality.stats(subsites=sbs) %>%
            get.neutrality(., theta=TRUE) %>%
            `[[`(1)
        Pi = input_data %>%
            F_ST.stats(mode='nucleotide', subsites=sbs) %>%
            `@`('nuc.diversity.within')
        output_list[[subsites]] = data.frame(stats1, Pi=Pi[, 1], 
                                             stringsAsFactors=FALSE)
    }
    output_list
}

adjust_basic_stats = function(stats, nsites) {
    for(i in seq_along(stats)) {
        if(grepl('theta|Pi', names(stats)[i])) {
            stats[i] = stats[i] / nsites
        }
    }
    stats
}

adjust_basic_stats_genes = function(stats, nsites) {
    gene_names = gsub(':.+', '', rownames(stats))
    nsites = nsites[match(gene_names, nsites$feature), 2]
    for(i in seq_along(stats)) {
        if(grepl('theta|Pi', names(stats)[i])) {
            stats[i] = stats[i] / nsites
        }
    }
    stats
}

write_output = function(stats_object, output_file) {
    for(i in seq_along(stats_object)) {
        subsites = names(stats_object)[i]
        write.table(data.frame('sites'=subsites,
                               region=rownames(stats_object[[i]]),
                               stats_object[[i]]),
                    file=output_file,
                    row.names=FALSE, quote=FALSE, col.names=(i == 1),
                    append=(i!=1), sep='\t')
    }
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
fastadir  = opts[['args']][1]
gfffile   = opts[['args']][2]
replicon  = opts[['args']][3]
popfile   = opts[['args']][4]
nsgfile   = opts[['args']][5]
nswfile   = opts[['args']][6]

# List of strains to include
strains = scan(popfile, what='character')

nsites_genes = read.csv(nsgfile, sep='\t', header=TRUE,
                        stringsAsFactors=FALSE,
                        colClasses=c('character', rep(NA, 7)))
nsites_wins = read.csv(nswfile, sep='\t', header=TRUE,
                       stringsAsFactors=FALSE)

# Read in the FASTA files and add SNP type annotations.
g = readData(path=fastadir,
            gffpath=dirname(gfffile),
            include.unknown=TRUE,
            big.data=TRUE) %>%
    set.populations(., new.populations=list(which(get.individuals(.)[[1]] %in% strains))) #%>%


## Diversity statistics for the entire genome
whole_genome_stats = basic_stats(g, c('all', 'genic', 'syn',
                                      'nonsyn', 'intergenic'))
whole_genome_nsites = colSums(nsites_wins[nsites_wins$replicon == replicon, -(1:3)])
whole_genome_stats[['all']] %<>% 
    adjust_basic_stats(whole_genome_nsites[['nsites']])
whole_genome_stats[['genic']] %<>% 
    adjust_basic_stats(whole_genome_nsites[['nsites_genic']])
whole_genome_stats[['nonsyn']] %<>% 
    adjust_basic_stats(whole_genome_nsites[['nsites_nonsyn']])
whole_genome_stats[['syn']] %<>% 
    adjust_basic_stats(whole_genome_nsites[['nsites_syn']])
whole_genome_stats[['intergenic']] %<>% 
    adjust_basic_stats(whole_genome_nsites[['nsites_intergenic']])
whole_genome_stats %>% write_output(paste0(outpref, '.genome.tsv'))


## Diversity statistics for genes
## Could add in g@region.data@Coding.matrix using
## list(parse_gff(gffRead(...), SNP.DATA=FALSE)$Coding)
genic_data = g %>%
    splitting.data(subsites='coding') %>%
    add_gene_names(gfffile, replicon)
genic_stats = basic_stats(genic_data, c('all', 'syn', 'nonsyn'))
x = adjust_basic_stats_genes(genic_stats[['all']],
                             nsites_genes[, c('feature', 'nsites')])
genic_stats[['all']] %<>%
    adjust_basic_stats_genes(nsites_genes[, c('feature', 'nsites')])
genic_stats[['syn']] %<>%
    adjust_basic_stats_genes(nsites_genes[, c('feature', 'nsites_syn')])
genic_stats[['nonsyn']] %<>%
    adjust_basic_stats_genes(nsites_genes[, c('feature', 'nsites_nonsyn')])
genic_stats %>% write_output(paste0(outpref, '.genes.tsv'))


## Diversity statistics for sliding windows
window_data = sliding.window.transform(g, width=winsize, jump=stepsize,
                                       type=2, start.pos=FALSE,
                                       end.pos=FALSE, whole.data=TRUE)
window_stats = basic_stats(window_data, c('all', 'syn', 'nonsyn'))
for(i in 1:nrow(window_stats[['all']])) {
    s = (i - 1) * stepsize
    e = s + winsize
    ns = nsites_wins[nsites_wins$replicon == replicon &
                     nsites_wins$pos >= s &
                     nsites_wins$pos + 250 <= e, 'nsites'] %>%
        sum()
    window_stats[['all']][i, ] %<>% adjust_basic_stats(ns)
}
for(i in 1:nrow(window_stats[['syn']])) {
    s = (i - 1) * stepsize
    e = s + winsize
    ns = nsites_wins[nsites_wins$replicon == replicon &
                     nsites_wins$pos >= s &
                     nsites_wins$pos + 250 <= e, 'nsites_syn'] %>%
        sum()
    window_stats[['syn']][i, ] %<>% adjust_basic_stats(ns)
}
for(i in 1:nrow(window_stats[['nonsyn']])) {
    s = (i - 1) * stepsize
    e = s + winsize
    ns = nsites_wins[nsites_wins$replicon == replicon &
                     nsites_wins$pos >= s &
                     nsites_wins$pos + 250 <= e, 'nsites_nonsyn'] %>%
        sum()
    window_stats[['nonsyn']][i, ] %<>% adjust_basic_stats(ns)
}
window_stats %>% write_output(paste0(outpref, '.windows.tsv'))

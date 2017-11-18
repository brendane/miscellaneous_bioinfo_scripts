#!/usr/bin/env Rscript
#
# Extract means and HPD intervals for PGE, PVE, and n_gamma
#

library(coda)
library(data.table)

argv = commandArgs(trailingOnly=TRUE)

thin = as.numeric(argv[2])
input_prefix = argv[1]
infiles = list.files(dirname(input_prefix),
                     pattern=paste0('^', basename(input_prefix)),
                     full.names=TRUE)
infiles = infiles[grepl('.hyp.txt', infiles)]

make_vec = function() numeric(length(infiles))
stats = data.frame('run'=make_vec(),
                   'pve_95_l'=make_vec(),
                   'pve_95_h'=make_vec(),
                   'pve_mean'=make_vec(),
                   'pge_95_l'=make_vec(),
                   'pge_95_h'=make_vec(),
                   'pge_mean'=make_vec(),
                   'ng_95_l'=make_vec(),
                   'ng_95_h'=make_vec(),
                   'ng_mean'=make_vec(),
                   'lmm_pve'=make_vec())
for(f in infiles) {
    run_number = as.numeric(gsub('\\.hyp\\.tx.+', '',
                                 gsub(input_prefix, '', f)))
    if(grepl('\\.gz$', f)) {
        x = fread(paste('zcat', f))
    } else {
        x = fread(f)
    }

    mcmc_data = as.mcmc(x, start=1, end=thin*nrow(x),
                        thin=thin)
    hpd_95 = HPDinterval(mcmc_data, 0.95)
    means = colMeans(x)

    stats[run_number+1, 'pve_95_l'] = hpd_95['pve', 'lower']
    stats[run_number+1, 'pve_95_h'] = hpd_95['pve', 'upper']
    stats[run_number+1, 'pge_95_l'] = hpd_95['pge', 'lower']
    stats[run_number+1, 'pge_95_h'] = hpd_95['pge', 'upper']
    stats[run_number+1, 'ng_95_l'] = hpd_95['n_gamma', 'lower']
    stats[run_number+1, 'ng_95_h'] = hpd_95['n_gamma', 'upper']
    stats[run_number+1, 'pve_mean'] = means['pve']
    stats[run_number+1, 'pge_mean'] = means['pge']
    stats[run_number+1, 'ng_mean'] = means['n_gamma']
    stats[run_number+1, 'run'] = run_number

    log_file = gsub('\\.hyp\\.tx.+', '.log.txt', f)
    l = scan(log_file, what='character', sep='\n', quiet=TRUE)
    stats[run_number+1, 'lmm_pve'] =
        as.numeric(gsub('.+= ', '', l[grepl('pve', l)][1]))


    cat('Done with', run_number, '\n', file=stderr())
}

stats = stats[order(stats[, 'run']), ]

write.table(stats, file=stdout(), sep='\t', col.names=TRUE,
            row.names=FALSE, quote=FALSE)

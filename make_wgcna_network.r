#!/usr/bin/env Rscript
#
# Use WGCNA to make a network.
#

library(data.table)
library(WGCNA)

options(stringsAsFactors=FALSE)

argv = commandArgs(trailingOnly=TRUE)

geno_file = argv[1]
ld_var_file = argv[2]
block_size = as.numeric(argv[3])
n_threads = as.numeric(argv[4])
trans_pow = as.numeric(argv[5])

enableWGCNAThreads(n_threads)

# Get list of variants
vars_to_use = scan(ld_var_file, what='character')

# Read in genotype data
geno_data = fread(geno_file, header=TRUE)

# Remove unwanted variants
geno_data = geno_data[geno_data$rs %in% vars_to_use, ]

# Transform into WGCNA format and eliminate variants with no
# variation
geno_mat = t(geno_data[, -(1:3)])
geno_sd = apply(geno_data, 2, sd, na.rm=TRUE)
geno_mat = geno_mat[, geno_sd > 0]

# Check for excessive missing data and remove samples and variants
# if necessary
gsg = goodSamplesGenes(geno_mat, verbose=3)
cat('No missing data issues?', gsg$allOK, '\n')
if(!gsg$allOK) {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
        printFlush(paste('Removing genes:', paste(names(geno_mat)[!gsg$goodGenes], collapse=', ')))
    if (sum(!gsg$goodSamples)>0)
        printFlush(paste('Removing samples:', paste(rownames(geno_mat)[!gsg$goodSamples], collapse=', ')))
    # Remove the offending genes and samples from the data:
    geno_mat = geno_mat[gsg$goodSamples, gsg$goodGenes]
}

# Make the network; there are lots of options to tinker with
# in this function
net = blockwiseModules(geno_mat, maxBlockSize=block_size,
                         power=trans_pow,
                         TOMType='unsigned', minModuleSize=10,
                         reassignThreshold=0, mergeCutHeight=0.25,
                         numericLabels=TRUE,
                         saveTOMs=TRUE,
                         saveTOMFileBase="output.wgcna",
                         verbose=3)

# Plotting
pdf('output.wgnca.dendrogram.pdf', width=12, height=9)
clrs = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], clrs[net$blockGenes[[1]]],
                    'Module colors',
                    dendroLabels=FALSE, hang=0.03,
                    addGuide=TRUE, guideHang=0.05)
dev.off()

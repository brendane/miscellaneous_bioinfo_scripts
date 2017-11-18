#!/usr/bin/env Rscript
#
# Using the methods described in the tutorial
#

library(data.table)
library(WGCNA)

options(stringsAsFactors=FALSE)

argv = commandArgs(trailingOnly=TRUE)

geno_file = argv[1]
ld_var_file = argv[2]
block_size = as.numeric(argv[3])
n_threads = as.numeric(argv[4])

enableWGCNAThreads(n_threads)

# Get list of variants
vars_to_use = scan(ld_var_file, what='character')

# Read in genotype data
geno_data = fread(geno_file, header=TRUE)

# Remove unwanted variants
geno_data = geno_data[geno_data$rs %in% vars_to_use, ]

# Transform into WGCNA format
geno_mat = t(geno_data[, -(1:3)])

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

# Make a sample dendrogram (the tutorial uses this step to
# check for outliers)
sample_dendro = hclust(dist(geno_mat), method='average')
pdf(file='sample_dendrogram.pdf', width=12, height=9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sample_dendro, main='Sample clustering to detect outliers',
     sub='', xlab='', cex.lab=1.5, cex.axis=1.5, cex.main=2)
dev.off()

# Run the threshold picking function
powers = c(1:10, seq(12, 20, 2))
sft = pickSoftThreshold(geno_mat, powerVector=powers,
                        verbose=5, blockSize=block_size)

# And make the plots
pdf('soft_threshold_plots.pdf', width=9, height=5)
par(mfrow=c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',
     ylab='Scale Free Topology Model Fit,signed R^2',
     type='n',
     main = paste('Scale independence'))
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     cex=cex1, col='red')
abline(h=0.90, col='red')
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab='Soft Threshold (power)',
     ylab='Mean Connectivity', type='n',
     main = paste('Mean connectivity'))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col='red')
dev.off()

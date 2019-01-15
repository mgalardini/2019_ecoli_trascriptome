#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
 
parser$add_argument("samples",
                    help="RNA samples file")
parser$add_argument("directory",
                    help="kallisto output directory")
parser$add_argument("straight",
                    help="VST output file")
parser$add_argument("corrected",
                    help="batch-corrected VST output file")

args <- parser$parse_args()

suppressPackageStartupMessages(library("tximport"))
suppressPackageStartupMessages(library('DESeq2'))
suppressPackageStartupMessages(library('limma'))


# prepare samples table
s2c <- read.table(file.path(args$samples),
                  header=TRUE,
                  stringsAsFactors=TRUE)
kal_dirs <- sapply(s2c$sample,
                   function(id) file.path(args$directory,
                                          strsplit(as.character(id), '_')[[1]][1],
					  strsplit(as.character(id), '_')[[1]][2],
                                          "abundance.h5"))
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# import kallisto's counts
txi.kallisto <- tximport(s2c$path,
                         type="kallisto",
                         txOut=TRUE,
                         countsCol='est_counts')
counts <- as.matrix(txi.kallisto$counts)
colnames(counts) <- s2c$sample
for (col in colnames(counts)){
    counts[, col] <- sapply(counts[, col], as.integer)
}
rownames(s2c) <- colnames(counts)

# apply variance stabilizing transformation (vst)
dds <- DESeqDataSetFromMatrix(counts,
			      s2c,
			      ~ 1)
tvst <- vst(dds)
write.table(assay(tvst),
	    file=args$straight,
	    sep="\t")

# correct for batch effects
# preserve treatment conditions
design <- model.matrix(~strain, s2c)
batch <- removeBatchEffect(assay(tvst),
                           batch=s2c$batch,
			   design=design)
write.table(batch,
	    file=args$corrected,
	    sep="\t")

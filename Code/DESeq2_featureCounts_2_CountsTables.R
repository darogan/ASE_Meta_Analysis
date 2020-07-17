#!/usr/local/bin/Rscript

# Analysis Performed by Russell S. Hamilton
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

# module load r-3.5.2-gcc-5.4.0-ordoh5p
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
#BiocManager::install("biomaRt")

#library("tidyr")
#library("dplyr")
library("methods")
library("utils")
library("Matrix")
library("matrixStats")
library("useful")
library("reshape")
library("reshape2")
library("DESeq2")
library("biomaRt")

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

Project <- args[1] 
#"Teichman_ESCs_NPCs"

#baseDir <- paste0("/home/rsh46/rds/rds-acf1004-afs-lab-rds/Imprinting/AFS_cae28_0001/", Project, "/PseudoReplicates/")
baseDir <- paste0("/home/rsh46/rds/rds-acf1004-afs-lab-rds/Imprinting/AFS_cae28_0001/", Project, "/")
setwd(baseDir)
print(baseDir)


message("+-------------------------------------------------------------------------------")
message("+ Use ensEMBL Annotations")
message("+-------------------------------------------------------------------------------")

ensembl    <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), mart = ensembl)  
head(ensEMBL2id)
nrow(ensEMBL2id)


message("+-------------------------------------------------------------------------------")
message("+ Set up a samples table")
message("+-------------------------------------------------------------------------------")

sampleFiles <- list.files(baseDir, pattern='*featureCounts_counts.txt$', recursive = TRUE)

sampleNames <- sampleFiles
sampleNames <- gsub(".bam_featureCounts_counts.txt", "", sampleNames)
sampleNames <- gsub("_1_val_1_CAST_EiJ_N-masked_hisat2_srtd", "", sampleNames)
sampleNames <- gsub("_trimmed_CAST_EiJ_N-masked_hisat2_srtd", "", sampleNames)
sampleGenome <- sampleNames
sampleGenome <- gsub(".*genome", "genome", sampleGenome)

sampleTable <- data.frame(sampleNames=sampleNames, fileNameDGE=sampleFiles, Genome=sampleGenome)
print(sampleTable)

write.csv(sampleTable, file=paste0(baseDir, Project, "_DESeq2_SampleTable.csv"), row.names=FALSE, quote=FALSE)

message("+-------------------------------------------------------------------------------")
message("+ Read in the featureCount gene count files per sample")
message("+-------------------------------------------------------------------------------")

DESeqDataSetFromFeatureCounts <- function (sampleTable, directory = ".", design, ignoreRank = FALSE, ...) 
{
  # From https://www.biostars.org/p/277316/
  if (missing(design)) 
    stop("design is missing")
  l <- lapply(as.character(sampleTable[, 2]), function(fn) read.table(file.path(directory, fn), skip=2))
  if (!all(sapply(l, function(a) all(a$V1 == l[[1]]$V1)))) 
    stop("Gene IDs (first column) differ between files.")
  tbl <- sapply(l, function(a) a$V7)
  colnames(tbl) <- sampleTable[, 1]
  rownames(tbl) <- l[[1]]$V1
  rownames(sampleTable) <- sampleTable[, 1]
  dds <- DESeqDataSetFromMatrix(countData = tbl, colData = sampleTable[, -(1:2), drop = FALSE], design = design, ignoreRank, ...)
  return(dds)
}


message("Read in featureCounts")

dds.Samples <- DESeqDataSetFromFeatureCounts(sampleTable=sampleTable, directory=baseDir, design= ~ Genome)


message("Run DESeq2 DDS")

dds.Samples <- DESeq(dds.Samples)


message("Counts")

Counts                    <- counts(dds.Samples, normalized=FALSE)
Counts.df                 <- as.data.frame(Counts)

write.csv(Counts.df, 
          file=paste0(baseDir, Project, "_DESeq2_ReadCounts.csv"),
          row.names=T, quote=FALSE)

Counts.df$ensembl_gene_id <- rownames(Counts.df)
Counts.df.annot           <- merge(ensEMBL2id, Counts.df, by="ensembl_gene_id" )

write.csv(Counts.df.annot, 
          file=paste0(baseDir, Project, "_DESeq2_ReadCounts.annot.csv"), 
          row.names=FALSE, quote=FALSE)


message("Normalised Counts")

normCounts                    <- counts(dds.Samples, normalized=TRUE)
normCounts.df                 <- as.data.frame(normCounts)

write.csv(normCounts.df, 
          file=paste0(baseDir, Project, "_DESeq2_NormalisedReadCounts.csv"), 
          row.names=T, quote=FALSE)

normCounts.df$ensembl_gene_id <- rownames(normCounts.df)
normCounts.df.annot           <- merge(ensEMBL2id, normCounts.df, by="ensembl_gene_id" )

write.csv(normCounts.df.annot, 
          file=paste0(baseDir, Project, "_DESeq2_NormalisedReadCounts.annot.csv"), 
          row.names=FALSE, quote=FALSE)

message("+-------------------------------------------------------------------------------")
message("+ END OF SCRIPT")
message("+-------------------------------------------------------------------------------")

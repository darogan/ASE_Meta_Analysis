#!/usr/local/bin/Rscript
#
#
#
#
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

library("tidyr")
library("dplyr")
library("methods")
library("utils")
library("ggplot2")
library("ggrepel")
library("cowplot")
library("Matrix")
library("matrixStats")
library("useful")
library("reshape")
library("reshape2")
library("DESeq2")
library("biomaRt")
library("ggforce")
library("Cairo")
library("pheatmap")
library('RColorBrewer')
library("eulerr")
library("scales")

baseDir <- "/Users/rhamilto/Documents/AFS-Projects"
setwd(baseDir)
print(baseDir)



table           <- read.csv(paste0(baseDir, "/", "babak.read.filter.txt"), header=F)
colnames(table) <- c("Sample", "Stage", "ReadCount")
table$Stage     <- gsub("genome1.bam",    "genome1",    table$Stage)
table$Stage     <- gsub("genome2.bam",    "genome2",    table$Stage)
table$Stage     <- gsub("unassigned.bam", "unassigned", table$Stage)
table$Stage     <- gsub("bam",            "aligned",    table$Stage)
head(table)

if(table$Stage == "fastq"){ table$ReadCount = table$ReadCount/4}
head(table)

table$Tissue <- table$Sample
table$Tissue <- gsub("SRR823449", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR823450", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR823469", "Liver", table$Tissue)
table$Tissue <- gsub("SRR823474", "Liver", table$Tissue)
table$Tissue <- gsub("SRR823485", "Hypothalamus", table$Tissue)
table$Tissue <- gsub("SRR823478", "Hypothalamus", table$Tissue)
table$Tissue <- gsub("SRR823461", "Cerebellum", table$Tissue)
table$Tissue <- gsub("SRR823458", "Cerebellum", table$Tissue)
table$Tissue <- gsub("SRR823472", "AdultWholeBrain", table$Tissue)
table$Tissue <- gsub("SRR823473", "AdultWholeBrain", table$Tissue)
head(table)


pdf(paste0(baseDir, "/ExploringReadFiltering_Babak",".pdf"),width=5,height=7, onefile=FALSE)
par(bg=NA)
ggplot(data=table, aes(x=Stage, y=ReadCount, colour=Tissue, group=Sample)) +
  geom_line(alpha=0.5, show.legend=FALSE) +
  geom_point() +
  scale_x_discrete(name="Stage", limits=c("fastq", "aligned", "unassigned", "genome1", "genome2")) +
  scale_y_continuous(name="Read Count", labels = comma) +
  ggtitle("Babak_2015") +
  theme_bw() +
  theme(legend.position = "top")
dev.off()



table           <- read.csv(paste0(baseDir, "/", "bonthuis.read.filter.txt"), header=F)
colnames(table) <- c("Sample", "Stage", "ReadCount")
table$Stage     <- gsub("genome1.bam",    "genome1",    table$Stage)
table$Stage     <- gsub("genome2.bam",    "genome2",    table$Stage)
table$Stage     <- gsub("unassigned.bam", "unassigned", table$Stage)
table$Stage     <- gsub("bam",            "aligned",    table$Stage)
head(table)

table$Tissue <- table$Sample
table$Tissue <- gsub("SRR823449", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086215", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086216", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086217", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086218", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086219", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086220", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086221", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086222", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086223", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086224", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086225", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086226", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086227", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086228", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086229", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086230", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086231", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086232", "ARN", table$Tissue)
table$Tissue <- gsub("SRR2086233", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086234", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086235", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086236", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086237", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086238", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086239", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086240", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086241", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086242", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086243", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086244", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086245", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086246", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086247", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086248", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086249", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086250", "DRN", table$Tissue)
table$Tissue <- gsub("SRR2086251", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086252", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086253", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086254", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086255", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086256", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086257", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086258", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086259", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086260", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086261", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086262", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086263", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086264", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086265", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086266", "Liver", table$Tissue)
table$Tissue <- gsub("SRR2086267", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086268", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086269", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086270", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086271", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086272", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086273", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086274", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086275", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086276", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086277", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086278", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086279", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086280", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086281", "Muscle", table$Tissue)
table$Tissue <- gsub("SRR2086282", "Muscle", table$Tissue)


pdf(paste0(baseDir, "/ExploringReadFiltering_Bonthuis",".pdf"),width=5,height=7, onefile=FALSE)
par(bg=NA)
ggplot(data=table, aes(x=Stage, y=ReadCount, colour=Tissue, group=Sample)) +
  geom_line(alpha=0.5, show.legend=FALSE) +
  geom_point() +
  scale_x_discrete(name="Stage", limits=c("fastq", "aligned", "unassigned", "genome1", "genome2")) +
  scale_y_continuous(name="Read Count", labels = comma) +
  ggtitle("Bonthuis_2015") +
  theme_bw() +
  theme(legend.position = "top")
dev.off()



message("+-------------------------------------------------------------------------------")
message("+ END OF SCRIPT ")
message("+-------------------------------------------------------------------------------")
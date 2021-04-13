library(BiocManager)
BiocManager::install("maEndToEnd")
library(gplots)
library(ggplot2)

#General Bioconductor packages
library(Biobase)
library(oligoClasses)

#Quality control and pre-processing packages
library(oligo)

#Analysis through linear models
BiocManager::install("limma")
library(limma)

#Annotation libraries
BiocManager::install("pd.hugene.2.1.st")
library(pd.hugene.2.1.st)
BiocManager::install("hugene21sttranscriptcluster.db")
library(hugene21sttranscriptcluster.db)

#Helper libraries
library(tidyverse)
library(ggrepel)

### Load in proper libraries for analysis
source("libraries.R")

#Read in cel files and store into a dataset
data <- ReadAffy()
eset <- rma(data)
celFiles <- list.celfiles()
affyRaw <- read.celfiles(celFiles)

#Normalize and log 2 transformed
eset <- oligo::rma(affyRaw)

### At this point we have a normalized expression set
### We have to add an annotation dataframe using a human transcript
### GPL17692-38203.txt was the annotation set that was used in this experiment

class(eset)
exprs(eset)

hist(eset)
boxplot(eset)

### Create correct sequences to extract fasted and fed samples from data
even_seq <- seq(2,22,2)
odd_seq <- seq(1,21,2)

fasted_data <- eset[,even_seq]
fed_data <- eset[,odd_seq]

### I want to look at the differentially expressed genes
### using the fed states as a reference

fed_means <- rowMeans(exprs(fed_data))
fasted_means <- rowMeans(exprs(fasted_data))

subtracted <- fed_means - fasted_means
sum(abs(subtracted)>1)

genes_different <- subtracted[abs(subtracted) > 1]
hist(genes_different)

genes_different <- as.data.frame(genes_different)

### Indices of genes that with gene expression levels
### Now must refer numbers to annotation matrix to 
### infer the gene that is expressed differently

gene_id <- rownames(genes_different)

annot <- data.frame(ACCNUM=sapply(contents(hugene21sttranscriptclusterACCNUM),
                                  paste, collapse=", "), 
                    SYMBOL=sapply(contents(hugene21sttranscriptclusterSYMBOL), paste, collapse=", "), 
                    DESC=sapply(contents(hugene21sttranscriptclusterGENENAME), paste, collapse=", "))
eset_data <- data.frame(exprs(eset))

symbol_data_combined <- merge(annot, eset_data, by.x=0, by.y=0, all=T)

### Removing all entry points that do not contain 
### a relevant gene symbol

complete_data <- symbol_data_combined[symbol_data_combined["SYMBOL"] != "NA", ]
rownames(complete_data) <- complete_data$Row.names
complete_data <- subset(complete_data, select = -c(Row.names))

class(complete_data)
dim(complete_data)

### Checking for symbol repeats to remove

dim(unique(complete_data['SYMBOL']))[1] == dim(complete_data)[1]

final_data <- complete_data[!duplicated(complete_data['SYMBOL']),]

### We now have a data set that is normalized and log2 transformed
### We have removed points that did not have GENE assignements
### We have also removed points that had duplicates for efficacy3 


################################################################################################
### Filtering lowly expressed genes
### Thru median soft intensity based filtering
### Using the median putatively as a cutoff for background intensities
### However once we look at the median (~3.6) we see that it 
### cuts out too much data therefore we will arbitrarily use 2 as a cutoff
### The logic beind using 2 is that we will get rid of a lot of 
### "Maintenance" genes that are expressed equally between all samples
### Therefore they would be seen to have the highest density in a hist plot

final_expressions <- as.matrix(final_data[,4:25])

final_medians <- rowMedians(final_expressions)

hist_res <- hist(final_medians, 100, col = "orange", freq = FALSE, 
                 main = "Histogram of the median intensities", 
                 border = "orange",
                 xlab = "Median intensities")

putative_bg_cutoff <- median(final_medians)

abline(v = putative_bg_cutoff, col = "coral4", lwd = 2)
abline(v = 2, col = "green", lwd = 2)

final_expressions$Medians = final_medians

final_expressions <- final_expressions[which(final_expressions$Medians > 2),]
final_expressions <- subset(final_expressions, select = -Medians)

### MA plot. 
### M <- log ratio (logx-logy)
### A <- mean average 1/2 (logx + logy)
### Positive log-ratio means higher expression
### Negative log-ratio means lower expression

fed_final_mean <- rowMeans(final_expressions[, odd_seq])
fasted_final_mean <- rowMeans(final_expressions[, even_seq])

mean_df <- data.frame(fed = fed_final_mean, fasted = fasted_final_mean)

plotMA(mean_df, main = "Fed vs Fasted MA Plot")
abline(h = c(-0.5,0.5), col = 'red', lwd = 1.5)

### subsetting eset to only have 
### expressions > 2

rownames <- rownames(final_expressions)
eset <- eset[which(rownames(eset) %in% rownames),]
eset

### Using linear model from limma to find DE genes
### Volcano plot to visualize values of the genes 
### Visualizes the DE genes based on their position in the plot
### log2FC >= abs(0.5) means DE gene

design <- data.frame(Fed = replicate(22,0))
design$Fasted = replicate(22,0)
design[even_seq, "Fasted"] = 1
design[odd_seq, "Fed"] = 1

fit <- lmFit(eset, design2)
cont.matrix <- makeContrasts(FastedvsFed=Fasted-Fed, levels=design2)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
table <- topTable(fit2, number = 200, adjust="BH", p.value = 0.05)

###Using AnnotationDbi to match probeIDs to gene symbols
###Allowing for annotation of genes in the summary table obtained from a linear model

table.10000 <- topTable(fit2, number = 20000, adjust="BH", p.value = 1)
table10000keys <- rownames(table.10000)

annot.table10000 <- AnnotationDbi::select(x=hugene21sttranscriptcluster.db, keys=table10000keys, 
                                          columns = c("SYMBOL","ENSEMBL", "GENENAME"))

annot.top10000 <- annot.table10000[which(!is.na(annot.table10000$ENSEMBL)),]
annot.top10000 <- subset(annot.table10000, select=-c(ENSEMBL))
table_DEgenes <- dplyr::distinct(annot.table10000)

probeID.10000 <- rownames(table.10000)
temp.10000 <- cbind(probeID.10000, table.10000)
tracemem(temp.10000) == tracemem(table.10000)
names(temp.10000)[1] <- "PROBEID"
colnames(temp.10000)

dplyr::inner_join(temp.10000, table_DEgenes, x1 = PROBEID)

temp.10000 <- dplyr::inner_join(temp.10000, table_DEgenes, x1 = PROBEID)
rownames(temp.10000) <- temp.10000[,1]
temp.10000 <- subset(temp.10000, select = -c(PROBEID))
temp.10000

temp.10000 <- subset(temp.10000, select=c(SYMBOL,logFC:GENENAME))
table.10000 <- temp.10000

### We now have a summary table that has indices that correlate
### to gene symbols instead of probeIDs
### We now will plot a volcano plot visualizing our results
### First must create labels allowing for proper categorization of gene expression

table.10000$expressed <- "NO"
table.10000$expressed[table.10000$logFC > 0.5 & table.10000$adj.P.Val < 0.05] <- "UP"
table.10000$expressed[table.10000$logFC < -0.5 & table.10000$adj.P.Val < 0.05] <- "DOWN"

colors <- c("blue", "red", "black")
names(colors) <- c("DOWN", "UP", "NO")
threshold2 <- threshold + scale_colour_manual(values = colors)

table.10000$delabel <- NA
table.10000$delabel[table.10000$expressed != 'NO'] <- 
  table.10000$SYMBOL[table.10000$expressed != 'NO']



###Filter summary table based on log fold change values
###Then create a legend for the top 5 genes that have their expressions
###decreased and the top 5 genes that have their expressions increased

table.10000 <- table.10000[order(table.10000$logFC),]
lowestfold5 <- table.10000$delabel[1:5]
highestfold5 <- table.10000[order(table.10000$logFC, decreasing = TRUE),]$delabel[1:5]


bottom_rownames <- rownames(subset(table.10000, SYMBOL %in% lowestfold5))
bottom_symbols <- (as.list(subset(table.10000, SYMBOL %in% lowestfold5)$SYMBOL))

top_rownames <- rownames(subset(table.10000, SYMBOL %in% highestfold5))
top_symbols <- (as.list(subset(table.10000, SYMBOL %in% highestfold5)$SYMBOL))

table.10000$delabel <- ''
table.10000[top_rownames, ]$delabel <- top_symbols
table.10000[bottom_rownames, ]$delabel <- bottom_symbols

ggplot(data=table.10000, aes(x=logFC, y=-log10(adj.P.Val), col=expressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  ggtitle("Volcano Plot for Gene Expression") +
  geom_text_repel(max.overlaps = 100000)

### GO analysis
### Purpose of this analysis is to extract
### Biologically important data from differentially
### expressed genes already found and creating a table to succinctly summarize the findings

top200keys <- rownames(table)

annot.top200 <- AnnotationDbi::select(x=hugene21sttranscriptcluster.db, keys=top200keys, 
                                      columns = c("SYMBOL","ENSEMBL", "GENENAME"))

annot.top200$PROBEID

annot.top200 <- annot.top200[which(!is.na(annot.top200$ENSEMBL)),]
annot.top200 <- subset(annot.top200, select=-c(ENSEMBL))
distinct_DEgenes <- dplyr::distinct(annot.top200)

probeID <- rownames(table)
temp <- cbind(probeID, table)
tracemem(temp) == tracemem(table)
names(temp)[1] <- "PROBEID"
colnames(temp)

dplyr::inner_join(temp, distinct_DEgenes, x1 = PROBEID)

temp <- dplyr::inner_join(temp, distinct_DEgenes, x1 = PROBEID)
rownames(temp) <- temp[,1]
temp <- subset(temp, select = -c(PROBEID))
temp

temp <- subset(temp, select=c(SYMBOL,logFC:GENENAME))
table <- temp


### Perform GO analysis using the Maayan Lab's Enrichr web-based software
### Going to perform different GO analysis on two groups
### Genes that were enhanced and genes that were repressed
### Genes are selected to be differentially expressed iff abs(logFC) > 0.5
### This is chosen because of the experiment design consisting of only 24 hour difference in feeding times

table_inc <- table[order(table$logFC),]
table_dec <- table[order(table$logFC,decreasing = TRUE),]
increased_DEgenes <- table_dec[which(table_dec['logFC'] > 0.5),]
decreased_DEgenes <- table_inc[which(table_inc['logFC'] < -0.5),]

gene_symbol_increased_cut<- as.list(as.character(increased_DEgenes$SYMBOL))
gene_symbol_decreased_cut<- as.list(as.character(decreased_DEgenes$SYMBOL))
lapply(gene_symbol_increased_cut, write, "CutIncreasedSymbols.txt", append=TRUE, ncolumns=1000)
lapply(gene_symbol_decreased_cut, write, "CutDecreasedSymbols.txt", append=TRUE, ncolumns=1000)

### Read in GO table
### Set correct column names
### Remove unnecessary columns

GO_increased <- read.table("./CutGO_Increased.txt", sep='\t')
new_names <- as.character(unlist(GO_increased[1, ], use.names = FALSE))
names(GO_increased) <- new_names
GO_increased <- GO_increased[-1,]
drop_columns <-  c("Old P-value", "Old Adjusted P-value", "Overlap", "Combined Score", "Adjusted P-value")
GO_increased <- GO_increased[,!(names(GO_increased) %in% drop_columns)]

GO_decreased <- read.table("./CutGO_Decreased.txt", sep='\t')
new_names <- as.character(unlist(GO_decreased[1,]), use.names = FALSE)
names(GO_decreased) <- new_names

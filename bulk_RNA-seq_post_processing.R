library(DESeq2)
library(tidyverse)
library(edgeR)
library(biomaRt)

setwd('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/bulk-rna-seq/')

count_table <- read.table('./output_pipeline/salmon.merged.gene_counts.tsv', header = T)

count_table <- count_table[,c(1,3,4,5,6,7,8,9,10)]

count_table$gene_id <- gsub('\\.[0-9]*$','',count_table$gene_id)

rownames(count_table) <- count_table$gene_id

count_table$gene_id <- NULL

colSums(count_table)

sampleinfo <- data.frame(readxl::read_excel("sampledata.xlsx"))

colnames(count_table) <- c("VKR_001", "VKR_002", "VKR_003", "VKR_004","VKR_005","VKR_006","VKR_007","VKR_008")

# filter only protein coding genes

listEnsembl()          # Lists current Ensembl versions
listEnsemblArchives()  # Lists archived versions of Ensembl
listMarts()  

# Connect to Ensembl and select the "hsapiens_gene_ensembl" dataset for mouse gene annotation
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get and store a list of filters (search constraints) available in the dataset
filters <- listFilters(ensembl) %>% as.data.frame()

# Use the row names (gene IDs) from the countdata as the list of genes to be annotated
genes <- rownames(count_table)
# 78,724 genes id in the count table

# Get and store a list of attributes (gene properties) that can be retrieved from Ensembl
attributes <- listAttributes(ensembl) %>% as.data.frame()

# Select attributes of interest from the attributes data frame (these columns hold the information you need for annotation)
att <- attributes$name[c(1, 31)]

# Select the appropriate filter for ensembl_gene_id_version, which filters genes based on Ensembl gene ID
filt <- filters$name[53]

# Retrieve gene information from Ensembl using the selected attributes and filter for the given list of genes
gene_info <- getBM(attributes = att, 
                   filters = filt,
                   values = genes,
                   mart = ensembl)

length(unique(gene_info$ensembl_gene_id))
# 78,724 genes id for which an annotation has been found 

all(gene_info$ensembl_gene_id %in% rownames(count_table))
# yes

all(rownames(count_table) %in% gene_info$ensembl_gene_id)
# yes

### Filter to keep only protein-coding genes
protein_coding <- gene_info[which(gene_info$gene_biotype == "protein_coding"),'ensembl_gene_id']

length(unique(protein_coding))
# 20,092 genes id in the count table are protein coding

# Filter countdata to keep only rows (genes) that match protein-coding genes in gene_info
count_table <- count_table[which(rownames(count_table) %in% protein_coding),]

colnames(count_table) <- c("microglia_noTGFB_rep1", "microglia_TGFB_rep1", "microglia_TGFB_rep2", "microglia_noTGFB_rep2", "EMP_tomato_rep1","EMP_tomato_rep2",
                           "EMP_no_tomato_rep1","EMP_no_tomato_rep2")

size <- colSums(count_table)
size
#microglia_noTGFB_rep1   microglia_TGFB_rep1   microglia_TGFB_rep2 microglia_noTGFB_rep2       EMP_tomato_rep1       EMP_tomato_rep2 
#6509883               8480255               8469296               6945713               9998810               3734707 
#EMP_no_tomato_rep1    EMP_no_tomato_rep2 
#8783354               9264558 

# These are the number of reads of each replicate that will be used for quantification. 
# We can see that the sizes are not so much variable.

count_table_only_microglia <- count_table[,c(1,4,2,3)]

y <- DGEList(counts=count_table_only_microglia)
y

# Assign the group label to each replicate, that is the tissue from which it derives
group <- as.factor(c("microglia_noTGFB", "microglia_noTGFB","microglia_TGFB", "microglia_TGFB"))
y$samples$group <- group
y

# Remove all the genes with low or zero counts (expression): first check how many genes do not appear in any of the 9 replicates
table(rowSums(y$counts==0)==4)
#FALSE  TRUE 
#16248  3844 

# keep.exprs function keeps all the genes that are expressed in all the 2 replicates of a given condition by removing
# those with zero or low expression since they have to be ignored during normalization and parameter estimation
keep.exprs <- filterByExpr(y, group=group)
# # default min count = 10
y <- y[keep.exprs, keep.lib.sizes=FALSE]
dim(y)
# 13591 genes that remain after the filtering and are those on which the DE analysis will be based on.

# Extract and store in a vector the log2 of the counts per million before normalization and plot their distribution
logcpm_before <- cpm(y, log=TRUE)
myColors <- ifelse(colnames(logcpm_before) %in% c("microglia_noTGFB_rep1","microglia_noTGFB_rep2") , 'mediumseagreen' ,'lightskyblue')
boxplot(logcpm_before,notch=T,xlab='Replicates',ylab='Log(CPM)', main='Log(CPM) before TMM normalization',col=myColors, varwidth=T)
legend("top", inset=.01, title="Condition", c("microglia_noTGFB","microglia_TGFB"), fill=c('mediumseagreen','lightskyblue'), horiz=TRUE, cex=0.35)

for (i in 1:9){
  print(median(logcpm_before[,i]))
}
#4.445544
# [1] 4.821608
# [1] 4.79559
# [1] 4.821738

# edgeR normalizes the counts using TMM normalization
y <- calcNormFactors(y, method = "TMM")
y
# TMM normalization is based on multiplying each count value of each sample for a constant, that is the normalization factor defined at
# sample level, trying to shift the count values gene by gene in order to have no change of expression across each pair of samples. The
# assumption of TMM normalization is that most of the genes do not change their expression in a significant way.

# Extract and store in a vector the log2 of the counts per million after normalization and plot their distribution
logcpm_after <- cpm(y, log=TRUE)
myColors <- ifelse(colnames(logcpm_before) %in% c("microglia_noTGFB_rep1","microglia_noTGFB_rep2") , 'mediumseagreen' ,'lightskyblue')
boxplot(logcpm_after,notch=T,xlab='Replicates',ylab='Log(CPM)', main='Log(CPM) after TMM normalization',col=myColors, varwidth=T)
legend("top", inset=.01, title="Condition", c("microglia_noTGFB","microglia_TGFB"), fill=c('mediumseagreen','lightskyblue'), horiz=TRUE, cex=0.35)

# The effect of TMM normalization is that now the medians are very similar
for (i in 1:9){
  print(median(logcpm_after[,i]))
}
#  4.642516
# [1] 4.750089
# [1] 4.742974
# [1] 4.748033

# Design the linear model without the intercept since it makes no sense to choose a type of tissue as reference being the 2 samples independent from one another
design <- model.matrix(~0 + group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design

pdf('./results_R/edgeR/MDS_only_microglia.pdf')
# MDS plotting the samples labeled by group
Shape <- c(15,15,17,17)
Color <- c("blue","blue","red","red")
plotMDS(logcpm_after, col=Color, pch=Shape,cex=2, main = 'Multidimensional scaling plot of distances\nbetween gene expression profiles') 
legend("bottomright", col=c("blue","red"), pch=c(15,17), legend=c("microglia_noTGFB1","microglia_TGFB1"))
# After normalization the different samples are projected in bidimensional space. The distance between points is the overall
# fold ratio gene by gene between two samples and it is computed only on the top 500 variable genes.
# The closer two points are the more similar are the expression values of these two samples: 
# the  replicates of the same condition should be close between each other and far from the replicates of the other condition
dev.off()

#plotMDS(logcpm_after, labels = y$samples$lib.size, main = 'Multidimensional scaling plot of distances\nbetween gene expression profiles') 


# Estimate the Negative Binomial dispersion and plot the BCV (square root of the dispersion)
y <- estimateDisp(y, design)
plotBCV(y, main = 'Biological Coefficient of Variation with respect to the average logCPM', ylim = c(0.15,1.75))
# The bigger is the BCV the higher is the dispersion, and so the variance, of the corresponding gene. The Common red
# line is the global NB dispersion φ: it is estimated from all the genes. But, one single dispersion value does not fit well all the
# genes and, on the other hand, the assumption is that we have too few replicates to have a reliable estimate of the dispersion of each gene.
# The Trend blue line is the estimated trend that tries to model the dependence between the mean expression and the
# dispersion.

# The estimated trend is used to shrinkage the dispersion of each gene towards the trend line itself so that the final gene-wise
# dispersion estimate is no more the observed dispersion of that gene, but  the original dispersion value of that gene modified by pulling it
# towards the estimated trend.

# the common BCV is higher than 0.5: it is quite high!!!

y$common.dispersion
# 0.32

head(y$trended.dispersion)
# 0.1279455 0.1428689 0.3656991 0.3189173 0.3537538 0.1707703

head(y$tagwise.dispersion)
# 0.15983811 0.09759181 0.31744505 0.31012429 0.43409925 0.11858597

# Fit our data to the “generalized linear” model we designed
fit <- glmQLFit(y, design)

# microglia TGFB vs microglia no TGFB 
micro_TGFB_micro_no <- glmQLFTest(fit, contrast = c(-1,1))

summary(decideTests(micro_TGFB_micro_no, p.value=0.01, lfc=1))
# -1*microglia_noTGFB 1*microglia_TGFB

#This is the table in which gene by gene the log2FC, the log2CPM of the same gene between the two samples and the p-value (not adjusted!!!)
#are reported.
micro_TGFB_micro_no <- topTags(micro_TGFB_micro_no, n=10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
micro_TGFB_micro_no <- micro_TGFB_micro_no$table

# library("org.Hs.eg.db")
# columns(org.Hs.eg.db)
# 
# # Extract the gene identifiers (Ensembl gene IDs) from the count matrix
# genes <- rownames(micro_TGFB_micro_no)
# 
# ann <- select(org.Hs.eg.db,keys=genes,columns='SYMBOL', keytype="ENSEMBL")
# 
# gene_names <- unlist(lapply(as.list(rownames(micro_TGFB_micro_no)), function(x) {
#   if (length(ann[ann$ENSEMBL == x, 'SYMBOL']) > 1) {
#     print(x)
#   }
#   
# }))

# Connect to Ensembl database and select the mouse dataset for annotation
ensembl <- useMart("ensembl",
                   dataset = "hsapiens_gene_ensembl")

# List available filters and attributes in the selected Ensembl dataset
listFilters(ensembl)
attributes <- listAttributes(ensembl) %>% as.data.frame()
filters <- listFilters(ensembl) %>% as.data.frame()


# Select attributes for annotation: gene ID, gene name, biotype, and other relevant fields
att <- attributes$name[c(1, 25,31)]

# Select the appropriate filter for gene annotation (Ensembl gene ID version)
filt <- filters$name[53]

genes <- rownames(micro_TGFB_micro_no)

# Retrieve gene information from Ensembl using the selected attributes and filter
gene_info <- getBM(attributes = att, 
                   filters = filt,
                   values = genes,
                   mart = ensembl)

micro_TGFB_micro_no$gene_name <- unlist(lapply(as.list(rownames(micro_TGFB_micro_no)), function(x) {
  gene_info[gene_info$ensembl_gene_id == x, 'external_gene_name']
}))

# Remove the genes with low expression, log2CPM < 0, since they are more likely false positives
micro_TGFB_micro_no <- micro_TGFB_micro_no[-(which(micro_TGFB_micro_no$logCPM<0)),]
# 4 genes

### Filter to keep only protein-coding genes
protein_coding <- gene_info[which(gene_info$gene_biotype == "protein_coding"),'ensembl_gene_id']

length(unique(protein_coding))
# 13591 genes id in the count table are protein coding

micro_TGFB_micro_no <- micro_TGFB_micro_no[rownames(micro_TGFB_micro_no) %in% protein_coding,]
# 13591 genes

nrow(micro_TGFB_micro_no[micro_TGFB_micro_no$PValue <=0.05 & micro_TGFB_micro_no$logFC >=1,])
nrow(micro_TGFB_micro_no[micro_TGFB_micro_no$PValue <=0.05 & micro_TGFB_micro_no$logFC <=-1,])
#          -1*microglia_noTGFB 1*microglia_TGFB
#Down                                    546
#Up                                      604


micro_TGFB_micro_no$diff <- 'Not'
micro_TGFB_micro_no[micro_TGFB_micro_no$PValue <= 0.05 & micro_TGFB_micro_no$logFC <= -1,'diff'] <- 'Microglia TGFB1-'
micro_TGFB_micro_no[micro_TGFB_micro_no$PValue <= 0.05 & micro_TGFB_micro_no$logFC >= 1,'diff'] <- 'Microglia TGFB1+'

write.table(micro_TGFB_micro_no, '/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/bulk-rna-seq/results_R/edgeR/DEGs_edgeR_microglia_TGFB+_micro_TGFB-.txt', quote = F, sep = '\t')

library(patchwork)
pdf('./results_R/edgeR/volcano_edgeR_microglia_TGFB+_micro_TGFB-.pdf', height = 7, width = 14)
ggplot2::ggplot(data=micro_TGFB_micro_no, aes(x=logFC, y=-log10(PValue), fill=diff)) +
  geom_point(shape=21, size=1.5) + 
  scale_fill_manual(values=c('Microglia TGFB1-'="#28658f", 'Microglia TGFB1+'="#ef857c")) +
  geom_point(data = subset(micro_TGFB_micro_no, diff != 'Not'), col = "black", shape=21, size=1.5) + 
  theme_classic() +
  ggrepel::geom_text_repel(aes(label=ifelse(PValue<0.005 & abs(logFC)>=1.5,gene_name,'')), max.overlaps = 50,xlim = c(-25,25), ylim = c(-Inf, Inf), min.segment.length = 0) +
  geom_vline(xintercept=c(-1, 1), col="red", linewidth=0.3) +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between Microglia TGFB1+ samples and Microglia TGFB1- samples')

dev.off()

## functional enrichment 
library('enrichR')
setEnrichrSite("Enrichr")
websiteLive <- TRUE

listEnrichrDbs()

dbs_pathway <- c("GO_Biological_Process_2021", "KEGG_2021_Human")
if (websiteLive) {
  enriched_pathway <- enrichr(micro_TGFB_micro_no[micro_TGFB_micro_no$diff == 'Microglia TGFB1+','gene_name'], dbs_pathway)
}

pdf('./results_R/edgeR/BP_microglia_TGFB.pdf', height = 4, width = 10)
if (websiteLive) plotEnrich(title = "Enriched terms of\nGO Biological Process 2021 db - Microglia TGFB1+", enriched_pathway[[1]], showTerms
                            = 15, numChar = 100, y = "Ratio", orderBy = "P.value")
dev.off()

pdf('./results_R/edgeR/KEGG_microglia_TGFB.pdf', height = 4, width = 10)
if (websiteLive) plotEnrich(title = "Enriched terms of\nKEGG 2021 Human db - Microglia TGFB1+", enriched_pathway[[2]], showTerms
                            = 15, numChar = 100, y = "Ratio", orderBy = "P.value")
dev.off()

dbs_pathway <- c("GO_Biological_Process_2021", "KEGG_2021_Human")
if (websiteLive) {
  enriched_pathway <- enrichr(micro_TGFB_micro_no[micro_TGFB_micro_no$diff == 'Microglia TGFB1-','gene_name'], dbs_pathway)
}

pdf('./results_R/edgeR/BP_microglia_no_TGFB.pdf', height = 4,width = 10)
if (websiteLive) plotEnrich(title = "Enriched terms of\nGO Biological Process 2021 db - Microglia TGFB1-", enriched_pathway[[1]], showTerms
                            = 15, numChar = 100, y = "Ratio", orderBy = "P.value")
dev.off()

pdf('./results_R/edgeR/KEGG_microglia_no_TGFB.pdf', height = 4,width = 10)
if (websiteLive) plotEnrich(title = "Enriched terms of\nKEGG 2021 Human db - Microglia TGFB1-", enriched_pathway[[2]], showTerms
                            = 15, numChar = 100, y = "Ratio", orderBy = "P.value")
dev.off()



######################################################

## All Microglia samples vs all EMP samples

count_table <- count_table[,c(1,4,2,3,5,6,7,8)]

y <- DGEList(counts=count_table)
y

# Assign the group label to each replicate, that is the tissue from which it derives
group <- as.factor(c("microglia", "microglia","microglia", "microglia","EMP","EMP","EMP","EMP"))
y$samples$group <- group
y

# Remove all the genes with low or zero counts (expression): first check how many genes do not appear in any of the 9 replicates
table(rowSums(y$counts==0)==8)
#FALSE  TRUE 
#17248  2844 

# keep.exprs function keeps all the genes that are expressed in all the 2 replicates of a given condition by removing
# those with zero or low expression since they have to be ignored during normalization and parameter estimation
keep.exprs <- filterByExpr(y, group=group)
# # default min count = 10
y <- y[keep.exprs, keep.lib.sizes=FALSE]
dim(y)
# 13417 genes that remain after the filtering and are those on which the DE analysis will be based on.

# Extract and store in a vector the log2 of the counts per million before normalization and plot their distribution
logcpm_before <- cpm(y, log=TRUE)
myColors <- ifelse(colnames(logcpm_before) %in% c("microglia_noTGFB_rep1", "microglia_noTGFB_rep2", "microglia_TGFB_rep1",  "microglia_TGFB_rep2") , 'mediumseagreen' ,'lightskyblue')
boxplot(logcpm_before,notch=T,xlab='Replicates',ylab='Log(CPM)', main='Log(CPM) before TMM normalization',col=myColors, varwidth=T)
legend("top", inset=.01, title="Condition", c("microglia","EMP"), fill=c('mediumseagreen','lightskyblue'), horiz=TRUE, cex=0.35)

for (i in 1:9){
  print(median(logcpm_before[,i]))
}
#4.476944
# [1] 4.853364
# [1] 4.827894
# [1] 4.865465
# [1] 4.899253
# [1] 4.754357
# [1] 4.642516
# [1] 4.662656

# edgeR normalizes the counts using TMM normalization
y <- calcNormFactors(y, method = "TMM")
y
# TMM normalization is based on multiplying each count value of each sample for a constant, that is the normalization factor defined at
# sample level, trying to shift the count values gene by gene in order to have no change of expression across each pair of samples. The
# assumption of TMM normalization is that most of the genes do not change their expression in a significant way.

# Extract and store in a vector the log2 of the counts per million after normalization and plot their distribution
logcpm_after <- cpm(y, log=TRUE)
myColors <- ifelse(colnames(logcpm_before) %in% c("microglia_noTGFB_rep1", "microglia_noTGFB_rep2", "microglia_TGFB_rep1",  "microglia_TGFB_rep2") , 'mediumseagreen' ,'lightskyblue')
boxplot(logcpm_after,notch=T,xlab='Replicates',ylab='Log(CPM)', main='Log(CPM) after TMM normalization',col=myColors, varwidth=T)
legend("top", inset=.01, title="Condition", c("microglia_noTGFB","microglia_TGFB"), fill=c('mediumseagreen','lightskyblue'), horiz=TRUE, cex=0.35)

# The effect of TMM normalization is that now the medians are very similar
for (i in 1:9){
  print(median(logcpm_after[,i]))
}
# 4.709685
# [1] 4.814195
# [1] 4.81222
# [1] 4.8259
# [1] 4.834737
# [1] 4.57894
# [1] 4.692785
# [1] 4.713924

# Design the linear model without the intercept since it makes no sense to choose a type of tissue as reference being the 2 samples independent from one another
design <- model.matrix(~0 + group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design

pdf('./results_R/edgeR/MDS_all_samples.pdf')
# MDS plotting the samples labeled by group
Shape <- c(15,15,16,16,17,17,18,18)
Color <- c("blue","blue","red","red","green","green",'pink','pink')
plotMDS(logcpm_after, col=Color, pch=Shape,cex=2, main = 'Multidimensional scaling plot of distances\nbetween gene expression profiles') 
legend("topleft", col=c("blue","red","green",'pink'), pch=c(15,16,17,18), legend=c("microglia_noTGFB1","microglia_TGFB1",'EMP_tomato','EMP_no_tomato'))
# After normalization the different samples are projected in bidimensional space. The distance between points is the overall
# fold ratio gene by gene between two samples and it is computed only on the top 500 variable genes.
# The closer two points are the more similar are the expression values of these two samples: 
# the  replicates of the same condition should be close between each other and far from the replicates of the other condition
dev.off()

#plotMDS(logcpm_after, labels = y$samples$lib.size, main = 'Multidimensional scaling plot of distances\nbetween gene expression profiles') 


# Estimate the Negative Binomial dispersion and plot the BCV (square root of the dispersion)
y <- estimateDisp(y, design)
plotBCV(y, main = 'Biological Coefficient of Variation with respect to the average logCPM', ylim = c(0.15,1.75))
# The bigger is the BCV the higher is the dispersion, and so the variance, of the corresponding gene. The Common red
# line is the global NB dispersion φ: it is estimated from all the genes. But, one single dispersion value does not fit well all the
# genes and, on the other hand, the assumption is that we have too few replicates to have a reliable estimate of the dispersion of each gene.
# The Trend blue line is the estimated trend that tries to model the dependence between the mean expression and the
# dispersion.

# The estimated trend is used to shrinkage the dispersion of each gene towards the trend line itself so that the final gene-wise
# dispersion estimate is no more the observed dispersion of that gene, but  the original dispersion value of that gene modified by pulling it
# towards the estimated trend.

# the common BCV is higher than 0.5: it is quite high!!!


# Fit our data to the “generalized linear” model we designed
fit <- glmQLFit(y, design)

# microglia TGFB vs microglia no TGFB 
micro_TGFB_micro_EMP <- glmQLFTest(fit, contrast = c(-1,1))

summary(decideTests(micro_TGFB_micro_EMP, p.value=0.01, lfc=1))
#-1*EMP 1*microglia

#This is the table in which gene by gene the log2FC, the log2CPM of the same gene between the two samples and the p-value (not adjusted!!!)
#are reported.
micro_TGFB_micro_EMP <- topTags(micro_TGFB_micro_EMP, n=10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
micro_TGFB_micro_EMP <- micro_TGFB_micro_EMP$table

# Connect to Ensembl database and select the mouse dataset for annotation
ensembl <- useMart("ensembl",
                   dataset = "hsapiens_gene_ensembl")

# List available filters and attributes in the selected Ensembl dataset
listFilters(ensembl)
attributes <- listAttributes(ensembl) %>% as.data.frame()
filters <- listFilters(ensembl) %>% as.data.frame()


# Select attributes for annotation: gene ID, gene name, biotype, and other relevant fields
att <- attributes$name[c(1, 25,31)]

# Select the appropriate filter for gene annotation (Ensembl gene ID version)
filt <- filters$name[53]

genes <- rownames(micro_TGFB_micro_EMP)

# Retrieve gene information from Ensembl using the selected attributes and filter
gene_info <- getBM(attributes = att, 
                   filters = filt,
                   values = genes,
                   mart = ensembl)

micro_TGFB_micro_EMP$gene_name <- unlist(lapply(as.list(rownames(micro_TGFB_micro_EMP)), function(x) {
  gene_info[gene_info$ensembl_gene_id == x, 'external_gene_name']
}))

# Remove the genes with low expression, log2CPM < 0, since they are more likely false positives
micro_TGFB_micro_EMP <- micro_TGFB_micro_EMP[-micro_TGFB_micro_EMP$logCPM<0,]
# 0 genes removed

### Filter to keep only protein-coding genes
protein_coding <- gene_info[which(gene_info$gene_biotype == "protein_coding"),'ensembl_gene_id']

length(unique(protein_coding))
# 13417 genes id in the count table are protein coding

micro_TGFB_micro_EMP <- micro_TGFB_micro_EMP[rownames(micro_TGFB_micro_EMP) %in% protein_coding,]
# 13417 genes

nrow(micro_TGFB_micro_EMP[micro_TGFB_micro_EMP$PValue <=0.05 & micro_TGFB_micro_EMP$logFC >=1,])
nrow(micro_TGFB_micro_EMP[micro_TGFB_micro_EMP$PValue <=0.05 & micro_TGFB_micro_EMP$logFC <=-1,])
#                         -1*EMP 1*microglia
#Down                                    2089
#Up                                      2665

nrow(micro_TGFB_micro_EMP[micro_TGFB_micro_EMP$FDR <=0.05 & micro_TGFB_micro_EMP$logFC >=1,])
nrow(micro_TGFB_micro_EMP[micro_TGFB_micro_EMP$FDR <=0.05 & micro_TGFB_micro_EMP$logFC <=-1,])
#                         -1*EMP 1*microglia
#Down                                    1551
#Up                                      2098

micro_TGFB_micro_EMP$diff <- 'Not'
micro_TGFB_micro_EMP[micro_TGFB_micro_EMP$PValue <= 0.05 & micro_TGFB_micro_EMP$logFC <= -1,'diff'] <- 'EMP'
micro_TGFB_micro_EMP[micro_TGFB_micro_EMP$PValue <= 0.05 & micro_TGFB_micro_EMP$logFC >= 1,'diff'] <- 'Microglia'

write.table(micro_TGFB_micro_EMP, '/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/bulk-rna-seq/results_R/edgeR/DEGs_edgeR_microglia_EMP.txt', quote = F, sep = '\t')

library(patchwork)
pdf('./results_R/edgeR/volcano_edgeR_microglia_EMP.pdf', height = 7, width = 14)
ggplot2::ggplot(data=micro_TGFB_micro_EMP, aes(x=logFC, y=-log10(PValue), fill=diff)) +
  geom_point(shape=21, size=1.5) + 
  scale_fill_manual(values=c('Microglia'="#28658f", 'EMP'="#ef857c")) +
  geom_point(data = subset(micro_TGFB_micro_EMP, diff != 'Not'), col = "black", shape=21, size=1.5) + 
  theme_classic() +
  ggrepel::geom_text_repel(aes(label=ifelse(FDR<0.00005 & abs(logFC)>=1.5,gene_name,'')), max.overlaps = 50,xlim = c(-40,40), ylim = c(-Inf, Inf), min.segment.length = 0) +
  geom_vline(xintercept=c(-1, 1), col="red", linewidth=0.3) +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between all Microglia samples and all EMP samples')

dev.off()

## functional enrichment 
library('enrichR')
setEnrichrSite("Enrichr")
websiteLive <- TRUE

listEnrichrDbs()

dbs_pathway <- c("GO_Biological_Process_2021", "KEGG_2021_Human")
if (websiteLive) {
  enriched_pathway <- enrichr(micro_TGFB_micro_EMP[micro_TGFB_micro_EMP$diff == 'Microglia','gene_name'], dbs_pathway)
}

pdf('./results_R/edgeR/BP_microglia_all.pdf', height = 4, width = 10)
if (websiteLive) plotEnrich(title = "Enriched terms of\nGO Biological Process 2021 db - Microglia", enriched_pathway[[1]], showTerms
                            = 15, numChar = 100, y = "Ratio", orderBy = "P.value")
dev.off()

pdf('./results_R/edgeR/KEGG_microglia_all.pdf', height = 4, width = 10)
if (websiteLive) plotEnrich(title = "Enriched terms of\nKEGG 2021 Human db - Microglia", enriched_pathway[[2]], showTerms
                            = 15, numChar = 100, y = "Ratio", orderBy = "P.value")
dev.off()

dbs_pathway <- c("GO_Biological_Process_2021", "KEGG_2021_Human")
if (websiteLive) {
  enriched_pathway <- enrichr(micro_TGFB_micro_EMP[micro_TGFB_micro_EMP$diff == 'EMP','gene_name'], dbs_pathway)
}

pdf('./results_R/edgeR/BP_EMP_all.pdf', height = 4,width = 10)
if (websiteLive) plotEnrich(title = "Enriched terms of\nGO Biological Process 2021 db - EMP", enriched_pathway[[1]], showTerms
                            = 15, numChar = 100, y = "Ratio", orderBy = "P.value")
dev.off()

pdf('./results_R/edgeR/KEGG_EMP_all.pdf', height = 4,width = 10)
if (websiteLive) plotEnrich(title = "Enriched terms of\nKEGG 2021 Human db - EMP", enriched_pathway[[2]], showTerms
                            = 15, numChar = 100, y = "Ratio", orderBy = "P.value")
dev.off()


######################################################

## EMP tomato samples vs EMP non tomato samples

count_table_only_emp <- count_table[,c(5,6,7,8)]
y <- DGEList(counts=count_table_only_emp)
y

# Assign the group label to each replicate, that is the tissue from which it derives
group <- as.factor(c("EMP tomato","EMP tomato","EMP no tomato","EMP no tomato"))
y$samples$group <- group
y

# Remove all the genes with low or zero counts (expression): first check how many genes do not appear in any of the 9 replicates
table(rowSums(y$counts==0)==4)
#FALSE  TRUE 
#16038  4054 

# keep.exprs function keeps all the genes that are expressed in all the 2 replicates of a given condition by removing
# those with zero or low expression since they have to be ignored during normalization and parameter estimation
keep.exprs <- filterByExpr(y, group=group)
# # default min count = 10
y <- y[keep.exprs, keep.lib.sizes=FALSE]
dim(y)
# 12967 genes that remain after the filtering and are those on which the DE analysis will be based on.

# Extract and store in a vector the log2 of the counts per million before normalization and plot their distribution
logcpm_before <- cpm(y, log=TRUE)
myColors <- ifelse(colnames(logcpm_before) %in% c("EMP_tomato_rep1", "EMP_tomato_rep2") , 'mediumseagreen' ,'lightskyblue')
boxplot(logcpm_before,notch=T,xlab='Replicates',ylab='Log(CPM)', main='Log(CPM) before TMM normalization',col=myColors, varwidth=T)
legend("top", inset=.01, title="Condition", c("EMP tomato","EMP no tomato"), fill=c('mediumseagreen','lightskyblue'), horiz=TRUE, cex=0.35)

for (i in 1:9){
  print(median(logcpm_before[,i]))
}
#4.993815
# [1] 4.882759
# [1] 4.738292
# [1] 4.762984

# edgeR normalizes the counts using TMM normalization
y <- calcNormFactors(y, method = "TMM")
y
# TMM normalization is based on multiplying each count value of each sample for a constant, that is the normalization factor defined at
# sample level, trying to shift the count values gene by gene in order to have no change of expression across each pair of samples. The
# assumption of TMM normalization is that most of the genes do not change their expression in a significant way.

# Extract and store in a vector the log2 of the counts per million after normalization and plot their distribution
logcpm_after <- cpm(y, log=TRUE)
myColors <- ifelse(colnames(logcpm_before) %in% c("EMP_tomato_rep1", "EMP_tomato_rep2") , 'mediumseagreen' ,'lightskyblue')
boxplot(logcpm_after,notch=T,xlab='Replicates',ylab='Log(CPM)', main='Log(CPM) after TMM normalization',col=myColors, varwidth=T)
legend("top", inset=.01, title="Condition", c("EMP tomato","EMP no tomato"), fill=c('mediumseagreen','lightskyblue'), horiz=TRUE, cex=0.35)

# The effect of TMM normalization is that now the medians are very similar
for (i in 1:9){
  print(median(logcpm_after[,i]))
}
# 4.941909
# [1] 4.822387
# [1] 4.796948
# [1] 4.816784

# Design the linear model without the intercept since it makes no sense to choose a type of tissue as reference being the 2 samples independent from one another
design <- model.matrix(~0 + group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design

pdf('./results_R/edgeR/MDS_only_EMP.pdf')
# MDS plotting the samples labeled by group
Shape <- c(16,16,18,18)
Color <- c("green","green",'pink','pink')
plotMDS(logcpm_after, col=Color, pch=Shape,cex=2, main = 'Multidimensional scaling plot of distances\nbetween gene expression profiles') 
legend("bottomleft", col=c("green",'pink'), pch=c(16,18), legend=c('EMP_tomato','EMP_no_tomato'))
# After normalization the different samples are projected in bidimensional space. The distance between points is the overall
# fold ratio gene by gene between two samples and it is computed only on the top 500 variable genes.
# The closer two points are the more similar are the expression values of these two samples: 
# the  replicates of the same condition should be close between each other and far from the replicates of the other condition
dev.off()

#plotMDS(logcpm_after, labels = y$samples$lib.size, main = 'Multidimensional scaling plot of distances\nbetween gene expression profiles') 


# Estimate the Negative Binomial dispersion and plot the BCV (square root of the dispersion)
y <- estimateDisp(y, design)
plotBCV(y, main = 'Biological Coefficient of Variation with respect to the average logCPM', ylim = c(0.15,1.75))
# The bigger is the BCV the higher is the dispersion, and so the variance, of the corresponding gene. The Common red
# line is the global NB dispersion φ: it is estimated from all the genes. But, one single dispersion value does not fit well all the
# genes and, on the other hand, the assumption is that we have too few replicates to have a reliable estimate of the dispersion of each gene.
# The Trend blue line is the estimated trend that tries to model the dependence between the mean expression and the
# dispersion.

# The estimated trend is used to shrinkage the dispersion of each gene towards the trend line itself so that the final gene-wise
# dispersion estimate is no more the observed dispersion of that gene, but  the original dispersion value of that gene modified by pulling it
# towards the estimated trend.

# the common BCV is higher than 0.5: it is quite high!!!


# Fit our data to the “generalized linear” model we designed
fit <- glmQLFit(y, design)

# microglia TGFB vs microglia no TGFB 
micro_TGFB_EMP_tomato_EMP_no <- glmQLFTest(fit, contrast = c(-1,1))

summary(decideTests(micro_TGFB_EMP_tomato_EMP_no, p.value=0.01, lfc=1))
# -1*EMP no tomato 1*EMP tomato

#This is the table in which gene by gene the log2FC, the log2CPM of the same gene between the two samples and the p-value (not adjusted!!!)
#are reported.
micro_TGFB_EMP_tomato_EMP_no <- topTags(micro_TGFB_EMP_tomato_EMP_no, n=10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
micro_TGFB_EMP_tomato_EMP_no <- micro_TGFB_EMP_tomato_EMP_no$table

# Connect to Ensembl database and select the mouse dataset for annotation
ensembl <- useMart("ensembl",
                   dataset = "hsapiens_gene_ensembl")

# List available filters and attributes in the selected Ensembl dataset
listFilters(ensembl)
attributes <- listAttributes(ensembl) %>% as.data.frame()
filters <- listFilters(ensembl) %>% as.data.frame()


# Select attributes for annotation: gene ID, gene name, biotype, and other relevant fields
att <- attributes$name[c(1, 25,31)]

# Select the appropriate filter for gene annotation (Ensembl gene ID version)
filt <- filters$name[53]

genes <- rownames(micro_TGFB_EMP_tomato_EMP_no)

# Retrieve gene information from Ensembl using the selected attributes and filter
gene_info <- getBM(attributes = att, 
                   filters = filt,
                   values = genes,
                   mart = ensembl)

micro_TGFB_EMP_tomato_EMP_no$gene_name <- unlist(lapply(as.list(rownames(micro_TGFB_EMP_tomato_EMP_no)), function(x) {
  gene_info[gene_info$ensembl_gene_id == x, 'external_gene_name']
}))

# Remove the genes with low expression, log2CPM < 0, since they are more likely false positives
micro_TGFB_EMP_tomato_EMP_no <- micro_TGFB_EMP_tomato_EMP_no[-micro_TGFB_EMP_tomato_EMP_no$logCPM<0,]
# 9 genes removed

### Filter to keep only protein-coding genes
protein_coding <- gene_info[which(gene_info$gene_biotype == "protein_coding"),'ensembl_gene_id']

length(unique(protein_coding))
# 12967 genes id in the count table are protein coding

micro_TGFB_EMP_tomato_EMP_no <- micro_TGFB_EMP_tomato_EMP_no[rownames(micro_TGFB_EMP_tomato_EMP_no) %in% protein_coding,]
# 12967 genes

nrow(micro_TGFB_EMP_tomato_EMP_no[micro_TGFB_EMP_tomato_EMP_no$PValue <=0.05 & micro_TGFB_EMP_tomato_EMP_no$logFC >=1,])
nrow(micro_TGFB_EMP_tomato_EMP_no[micro_TGFB_EMP_tomato_EMP_no$PValue <=0.05 & micro_TGFB_EMP_tomato_EMP_no$logFC <=-1,])
#                    -1*EMP no tomato 1*EMP tomato
#Down                                    909
#Up                                      1774

nrow(micro_TGFB_EMP_tomato_EMP_no[micro_TGFB_EMP_tomato_EMP_no$FDR <=0.05 & micro_TGFB_EMP_tomato_EMP_no$logFC >=1,])
nrow(micro_TGFB_EMP_tomato_EMP_no[micro_TGFB_EMP_tomato_EMP_no$FDR <=0.05 & micro_TGFB_EMP_tomato_EMP_no$logFC <=-1,])
#                         -1*EMP no tomato 1*EMP tomato
#Down                                    74
#Up                                      318

micro_TGFB_EMP_tomato_EMP_no$diff <- 'Not'
micro_TGFB_EMP_tomato_EMP_no[micro_TGFB_EMP_tomato_EMP_no$PValue <= 0.05 & micro_TGFB_EMP_tomato_EMP_no$logFC <= -1,'diff'] <- 'EMP tomato-'
micro_TGFB_EMP_tomato_EMP_no[micro_TGFB_EMP_tomato_EMP_no$PValue <= 0.05 & micro_TGFB_EMP_tomato_EMP_no$logFC >= 1,'diff'] <- 'EMP tomato+'

write.table(micro_TGFB_EMP_tomato_EMP_no, '/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/bulk-rna-seq/results_R/edgeR/DEGs_edgeR_EMP_tomato_VS_EMP_no_tomato.txt', quote = F, sep='\t')

library(patchwork)
pdf('./results_R/edgeR/volcano_edgeR_EMP_tomato_VS_EMP_no_tomato.pdf', height = 7, width = 14)
ggplot2::ggplot(data=micro_TGFB_EMP_tomato_EMP_no, aes(x=logFC, y=-log10(PValue), fill=diff)) +
  geom_point(shape=21, size=1.5) + 
  scale_fill_manual(values=c('EMP tomato-'="#28658f", 'EMP tomato+'="#ef857c")) +
  geom_point(data = subset(micro_TGFB_EMP_tomato_EMP_no, diff != 'Not'), col = "black", shape=21, size=1.5) + 
  theme_classic() +
  ggrepel::geom_text_repel(aes(label=ifelse(FDR<0.02 & abs(logFC)>=1.5,gene_name,'')), max.overlaps = 50,xlim = c(-40,40), ylim = c(-Inf, Inf), min.segment.length = 0) +
  geom_vline(xintercept=c(-1, 1), col="red", linewidth=0.3) +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between all EMP tomato+ samples and EMP tomato- samples')

dev.off()

## functional enrichment 
library('enrichR')
setEnrichrSite("Enrichr")
websiteLive <- TRUE

listEnrichrDbs()

dbs_pathway <- c("GO_Biological_Process_2021", "KEGG_2021_Human")
if (websiteLive) {
  enriched_pathway <- enrichr(micro_TGFB_EMP_tomato_EMP_no[micro_TGFB_EMP_tomato_EMP_no$diff == 'EMP tomato+','gene_name'], dbs_pathway)
}

pdf('./results_R/edgeR/BP_EMP_tomato+.pdf', height = 4, width = 10)
if (websiteLive) plotEnrich(title = "Enriched terms of\nGO Biological Process 2021 db - EMP tomato+", enriched_pathway[[1]], showTerms
                            = 15, numChar = 100, y = "Ratio", orderBy = "P.value")
dev.off()

pdf('./results_R/edgeR/KEGG_EMP_tomato+.pdf', height = 4, width = 10)
if (websiteLive) plotEnrich(title = "Enriched terms of\nKEGG 2021 Human db - EMP tomato+", enriched_pathway[[2]], showTerms
                            = 15, numChar = 100, y = "Ratio", orderBy = "P.value")
dev.off()

dbs_pathway <- c("GO_Biological_Process_2021", "KEGG_2021_Human")
if (websiteLive) {
  enriched_pathway <- enrichr(micro_TGFB_EMP_tomato_EMP_no[micro_TGFB_EMP_tomato_EMP_no$diff == 'EMP tomato-','gene_name'], dbs_pathway)
}

pdf('./results_R/edgeR/BP_EMP_tomato-.pdf', height = 4,width = 10)
if (websiteLive) plotEnrich(title = "Enriched terms of GO Biological\nProcess 2021 db - EMP tomato-", enriched_pathway[[1]], showTerms
                            = 15, numChar = 100, y = "Ratio", orderBy = "P.value")
dev.off()

pdf('./results_R/edgeR/KEGG_EMP_tomato-.pdf', height = 4,width = 10)
if (websiteLive) plotEnrich(title = "Enriched terms of\nKEGG 2021 Human db - EMP tomato-", enriched_pathway[[2]], showTerms
                            = 15, numChar = 100, y = "Ratio", orderBy = "P.value")
dev.off()

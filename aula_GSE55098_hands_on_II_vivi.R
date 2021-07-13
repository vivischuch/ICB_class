
# This is a comment

dias <- c(1, 5, 10, 15, 20, 24)
mortes <- c(1, 2, 3, 4, 7, 9)

class(dias)
class(mortes)

help(plot)
?plot

plot(dias, mortes) 
plot(dias, mortes, main="Scatterplot", xlab="Dias", ylab="Mortes")
plot(dias, mortes, main="Scatterplot", xlab="Dias", ylab="Mortes", las=1)
plot(dias, mortes, main="Scatterplot", xlab="Dias", ylab="Mortes", las=1, pch=8)
plot(dias, mortes, main="Scatterplot", xlab="Dias", ylab="Mortes", las=1, pch=8, col=2)
plot(dias, mortes, main="Scatterplot", xlab="Dias", ylab="Mortes", las=1, pch=8, col=2, frame.plot = FALSE)
plot(dias, mortes, main="Scatterplot", xlab="Dias", ylab="Mortes", las=1, pch=8, col=2, frame.plot = FALSE, xlim=c(0,25), ylim=c(0,10))
plot(dias, mortes, main="Scatterplot", xlab="Dias", ylab="Mortes", las=1, pch=8, col=2, frame.plot = FALSE, xlim=c(0,25), ylim=c(0,10), yaxs="i", xaxs="i")

plot(dias, mortes, # plot a scatter chart of dias and mortes
     main="Scatterplot", # add a title label
     xlab="Dias",  # add x-axis label
     ylab="Mortes", # add y-axis label
     las=1, # rotate the values on the y-axis
     pch=8, # change the plotting character
     col=2, # change the color
     frame.plot = FALSE, # remove top and right borders
     xlim=c(0,25), # change x-axis limits
     ylim=c(0,10), # change y-axis limits
     yaxs="i", # set the y-axis limits to exact values
     xaxs="i") # set the x-axis limits to exact values

# ++++++++

tabela <- GSE55098_toptable

str(tabela) # display the internal structure of a data frame
tabela$Gene.symbol # $ operator to address a particular column
head(tabela) # display the first rows

mean(tabela$logFC) # calculate the average of numbers of logFC column
median(tabela$logFC) # calculate the median of numbers of logFC column

tabela$log10P<- log10(tabela$P.Value) # calculate the log10 of the P-value
tabela$minuslog10P <- -log10(tabela$P.Value) # calculate the -log10 of the P-value

degs <- tabela[tabela$P.Value < 0.05,]  # get the DEGs based on statistical significance
up <- degs[degs$logFC > 1,] # Based on their fold change, select up regulated genes
down <- degs[degs$logFC < -1,] #and down regulated genes


write.table(degs, "degs_DEGs_GSE55098.tsv", sep='\t', row.names = F, quote = F)
write.table(up, "up_DEGs_GSE55098.tsv", sep='\t', row.names = F, quote = F)
write.table(down, "down_DEGs_GSE55098.tsv", sep='\t', row.names = F, quote = F)


plot(tabela$logFC, tabela$minuslog10P)
plot(tabela$logFC, tabela$minuslog10P, main="Volcano Plot", xlab="log2FC", ylab="-log10P")
plot(tabela$logFC, tabela$minuslog10P, main="Volcano Plot", xlab="log2FC", ylab="-log10P", las=1)
plot(tabela$logFC, tabela$minuslog10P, main="Volcano Plot", xlab="log2FC", ylab="-log10P", las=1, pch=20)
plot(tabela$logFC, tabela$minuslog10P, main="Volcano Plot", xlab="log2FC", ylab="-log10P", las=1, pch=20, xlim=c(-2.5,2))

abline(h=1.3, col="blue") # add horizontal straight line at the cut-off
abline(v=c(-1,1), col="blue") # add vertical straight lines

# Add grey points to unperturbed genes
with(subset(tabela, minuslog10P < 1.3), points(logFC, minuslog10P, pch=20, col="gray")) 
# Add green points to downregulated genes
with(subset(tabela, logFC < -1 & minuslog10P > 1.3), points(logFC, minuslog10P, pch=20, col="blue")) 
# Add red points to upregulated genes
with(subset(tabela, logFC > 1 & minuslog10P > 1.3), points(logFC, minuslog10P, pch=20, col="red")) 

# install.packages("calibrate")
library(calibrate)

with(subset(tabela, logFC < -1 & minuslog10P >1.3), textxy(logFC, minuslog10P, labs=Gene.symbol, cex=.8, col="blue"))
with(subset(tabela, logFC > 1 & minuslog10P >1.3), textxy(logFC, minuslog10P, labs=Gene.symbol, cex=.8, col="red"))

#+++++ FIG

png(file = "GSE55098_volcano_plot.png", units = "in", width = 12, height = 8, pointsize = 10, res = 100)
plot(tabela$logFC, tabela$minuslog10P, main="GSE55098 Volcano Plot", xlab="log2FC", ylab="-log10P", las=1, pch=20, xlim=c(-2.5,2))
abline(h=1.3, col="blue") # add horizontal straight line at the cut-off
abline(v=c(-1,1), col="blue") # add vertical straight lines
with(subset(tabela, minuslog10P < 1.3), points(logFC, minuslog10P, pch=20, col="gray")) 
with(subset(tabela, logFC < -1 & minuslog10P > 1.3), points(logFC, minuslog10P, pch=20, col="blue")) 
with(subset(tabela, logFC > 1 & minuslog10P > 1.3), points(logFC, minuslog10P, pch=20, col="red")) 
with(subset(tabela, logFC < -1 & minuslog10P >1.3), textxy(logFC, minuslog10P, labs=Gene.symbol, cex=.8, col="blue"))
with(subset(tabela, logFC > 1 & minuslog10P >1.3), textxy(logFC, minuslog10P, labs=Gene.symbol, cex=.8, col="red"))
dev.off()
# ++++++++

# if (!requireNamespace("BiocManager", quietly = TRUE))
#         install.packages("BiocManager")
# BiocManager::install("GEOquery")

library(GEOquery)

gse <- getGEO("GSE55098", GSEMatrix = TRUE)
eset <- gse[[1]] # Primeiro ExpressionSet
expr <- as.data.frame(Biobase::exprs(eset))
feature_data <- Biobase::fData(eset)
str(feature_data)
expr$Symbol <- feature_data$`Gene Symbol`
head(expr)

rownames(expr) <- NULL
head(expr)

pheno_data <- Biobase::pData(eset)
head(pheno_data)
write.table(pheno_data, "phenodata_GSE55098.tsv", sep='\t', row.names = F, quote = F)
write.table(expr, "GSE55098_expr_forMDP_entrez.txt", sep='\t', row.names = F, quote = F)

# +++++++ MDP

expressionFile <- "GSE55098_expr_forMDP_entrez.txt"
annotationFile <- "GSE55098_f.csv"

install.packages("RCurl")
library(RCurl)

exp1<- read.table( expressionFile, header = T, na.strings = "NA", sep = "\t", row.names = NULL )
samplesinfo <- read.csv(text = getURL("https://raw.githubusercontent.com/vivischuch/ICB_class/main/GSE55098_f.csv"), header = T, na.strings = "NA", sep = "\t")

controlGroup <- "control"
treatedGroup <- "test"

#Remove samples that are not in comparison
samplesinfo <- samplesinfo[ c( which( samplesinfo$Class == controlGroup ),
                               which( samplesinfo$Class == treatedGroup ) ), ]

samplesinfo <- samplesinfo[ samplesinfo$Sample %in% colnames( exp1 ), ]

exp1 <- exp1[ ,c( "Symbol", as.character( samplesinfo$Sample ) ) ]
samplesinfo <- samplesinfo[ samplesinfo$Sample %in% colnames( exp1 ), ]

#Colapse symbols that are duplicated by taking the one with higest expression
exp1$meanG <- apply( exp1[ , 2:ncol( exp1 ) ], 1, mean )
exp1 <- exp1[ order( exp1[ , 'Symbol' ],
                     exp1[ , 'meanG' ] ), ]


exp1$Symbol <- sapply(exp1$Symbol, function(Symbol) gsub("\\s?\\/\\/\\/\\s?.*$", "", Symbol))
exp1 <- exp1[ !duplicated( exp1$Symbol ), ]
rownames(exp1) <- NULL

# Get all sample without mean column (last column)
exp2 <- exp1[ , 2:( ncol( exp1 ) - 1 ) ]
NumberOfSamplesWithZeroAllowed <- round( 0.25 * ncol( exp2 ) )
exp2 <- exp2[rowSums( exp2 == 0 ) <= NumberOfSamplesWithZeroAllowed, ]


# Rodando MDP
# BiocManager::install("mdp")
library(mdp)
mdp.results <- mdp( data = exp2, pdata = samplesinfo, control_lab = controlGroup,
                    file_name = paste( "GSE55098_MDP_", treatedGroup, "_", controlGroup, "_", sep = "" ) )

# ++++++ Cemitool

## ----style, echo=FALSE, results="asis", message=FALSE-------------------------
knitr::opts_chunk$set(tidy = FALSE,
                      warning = FALSE,
                      message = FALSE,
                      cache=TRUE)

## -----------------------------------------------------------------------------
# BiocManager::install("CEMiTool")
library("CEMiTool")

# load expression data
expr0<- exp2
head(expr0[,1:4])

## -----------------------------------------------------------------------------
# load your sample annotation data
sample_annot <- read.table("sample_annotation_GSE55098.csv", sep = "\t", header = T, stringsAsFactors = FALSE)
head(sample_annot)
sample_annot <- sample_annot[,c(1, 2)]
colnames(sample_annot)[1] <- "SampleName"
head(sample_annot)

# Put the expression matrix in same order and with same samples of the sample annotation file
index <- match(sample_annot$SampleName, colnames(expr0))
exprs_matrix_reordered <- expr0[,index]
table(sample_annot$SampleName == colnames(exprs_matrix_reordered))

# Load the gene sets file
gene_sets <- read_gmt("BTM_for_GSEA_20131008.gmt")

# gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
# gmt_in <- read_gmt(gmt_fname)

# Load the CEMiTool gene interactions file
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_ppi <- read.delim(int_fname)

# Run CEMiTool and save the results
cem <- cemitool(expr = exprs_matrix_reordered, annot = sample_annot, gmt = gene_sets, interactions = int_ppi, filter = TRUE, plot = TRUE, verbose = TRUE)

show_plot(cem, "gsea")

# Save all results in reports, tables, and figures
# Generate an html report with all CEMiTool results
# It is saved in the directory called "Reports/Report"
generate_report(cem, force = TRUE)
# Generate several files with the results in tables
# It is saved in the directory called "Tables"
write_files(cem, force = TRUE)
# Generate several plots in pdf
# It is saved in the directory called "Plots"
save_plots(cem, "all", force = TRUE)

# Generates the diagnostic report
# It is saved in the directory called "Reports/Diagnostics"
diagnostic_report(cem, force = TRUE)

                      
# ++++++ plot number of samples using ggplot

sample_annot <- read.table("sample_annotation_GSE55098.csv", sep = "\t", header = T, stringsAsFactors = FALSE)
head(sample_annot)
sample_annot <- sample_annot[,c(1, 2)]
colnames(sample_annot)[1] <- "SampleName"
head(sample_annot)

a <- as.data.frame(table(sample_annot$Class))
head(a)
a <- a[order(a$Freq),]

library(tidyverse)
# BiocManager::install("tidyverse")
a <- a %>% mutate(row = row_number())
write.table(a, "GSE55098_number_samples.tsv", sep='\t', row.names = F, quote = F)

# ++++ FIG ++++
pdf(file = "GSE55098_Number_of_samples_barplot.pdf", width=6,height=4)
p<-ggplot(data=a, aes(x=reorder(Var1, -row), y=Freq, fill=Var1)) +
        geom_bar(stat="identity", aes(fill = Var1))+
        scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
        geom_text(aes(label=Freq), hjust=1.6, color="white", size=8)+
        ggtitle("GSE55098 Number of Samples") +
        labs(x = "", y = "") +
        scale_y_continuous(limits = c(0, 26), breaks = c(0, 5, 10, 15, 20, 25))+
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_blank(),
                plot.title = element_text(color = "black", size = 18),
                axis.text.x =  element_text(colour = "#333333", size = 14),
                axis.text.y =  element_text(color = "#333333", size = 18),
                legend.position="none",
                axis.ticks.y = element_blank()
        )+
        coord_flip()
p
dev.off()

png(file = "GSE55098_Number_of_samples_barplot.png", units = "in", width = 6, height = 4, pointsize = 10, res = 100)
p<-ggplot(data=a, aes(x=reorder(Var1, -row), y=Freq, fill=Var1)) +
        geom_bar(stat="identity", aes(fill = Var1))+
        scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
        geom_text(aes(label=Freq), hjust=1.6, color="white", size=8)+
        ggtitle("GSE55098 Number of Samples") +
        labs(x = "", y = "") +
        scale_y_continuous(limits = c(0, 26), breaks = c(0, 5, 10, 15, 20, 25))+
        theme_bw() +
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border = element_blank(),
                plot.title = element_text(color = "black", size = 18),
                axis.text.x =  element_text(colour = "#333333", size = 14),
                axis.text.y =  element_text(color = "#333333", size = 18),
                legend.position="none",
                axis.ticks.y = element_blank()
        )+
        coord_flip()
p
dev.off()





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

# tabela <- read.csv(sprintf("https://docs.google.com/uc?id=%s&export=download", "1B3Zvti2KqCt15LopJxPmB1yZbTARUW95"), sep = "\t")

tabela <- GSE122737_top_table

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
with(subset(tabela, logFC < -1 & minuslog10P > 1.3), points(logFC, minuslog10P, pch=20, col="green")) 
# Add red points to upregulated genes
with(subset(tabela, logFC > 1 & minuslog10P > 1.3), points(logFC, minuslog10P, pch=20, col="red")) 

# install.packages("calibrate")
library(calibrate)
#label down genes
with(subset(tabela, logFC < -1 & minuslog10P >1.3), textxy(logFC, minuslog10P, labs=Gene.symbol, cex=.8, col="green")) 
#label up genes
with(subset(tabela, logFC > 1 & minuslog10P >1.3), textxy(logFC, minuslog10P, labs=Gene.symbol, cex=.8, col="red")) 


# ++++++++


if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
BiocManager::install("GEOquery")

library(GEOquery)

gse <- getGEO("GSE122737", GSEMatrix = TRUE)
eset <- gse[[1]] # Primeiro ExpressionSet
expr <- as.data.frame(Biobase::exprs(eset))
feature_data <- Biobase::fData(eset)
str(feature_data)
expr$Symbol <- feature_data$Symbol
head(expr)

rownames(expr) <- NULL
head(expr)

pheno_data <- Biobase::pData(eset)
head(pheno_data)
write.table(pheno_data, "phenodata_GSE122737.tsv", sep='\t', row.names = F, quote = F)

# ++++++++++

dir.create(paste0(workingdir,"/MDP"))
write.table(expr, paste0(workingdir,"/MDP/","GSE1402_expr_forMDP_entrez.txt"), sep='\t', row.names = F, quote = F)

       
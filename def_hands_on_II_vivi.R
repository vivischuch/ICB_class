
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

tabela <- GSE122737_toptable

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


# if (!requireNamespace("BiocManager", quietly = TRUE))
#         install.packages("BiocManager")
# BiocManager::install("GEOquery")

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
# expr$Symbol <- as.character((expr$Symbol))


# ++++++


pheno_data <- Biobase::pData(eset)
head(pheno_data)
write.table(pheno_data, "phenodata_GSE1402.tsv", sep='\t', row.names = F, quote = F)


dir.create("number_samples")
setwd("number_samples")

targets <- GSE122737_f

a <- as.data.frame(table(targets$Class))
head(a)
# a$Freq <- as.numeric(a$Freq)
# head(a)
a <- a[order(a$Freq),]

library(tidyverse)
# BiocManager::install("tidyverse")
a <- a %>% mutate(row = row_number())
write.table(a, "GSE1402_number_samples.tsv", sep='\t', row.names = F, quote = F)

# ++++ FIG ++++
pdf(file = "GSE1402_Number_of_samples_barplot.pdf", width=6,height=4)
p<-ggplot(data=a, aes(x=reorder(Var1, -row), y=Freq, fill=Var1)) +
        geom_bar(stat="identity", aes(fill = Var1))+
        scale_fill_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
        geom_text(aes(label=Freq), hjust=1.6, color="white", size=8)+
        ggtitle("GSE1402 Number of Samples") +
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

png(file = "GSE1402_Number_of_samples_barplot.png", units = "in", width = 6, height = 4, pointsize = 10, res = 100)
p<-ggplot(data=a, aes(x=reorder(Var1, -row), y=Freq, fill=Var1)) +
        geom_bar(stat="identity", aes(fill = Var1))+
        scale_fill_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
        geom_text(aes(label=Freq), hjust=1.6, color="white", size=8)+
        ggtitle("GSE1402 Number of Samples") +
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

# +++++++++++

stringsAsFactors=FALSE
workingdir <- "/home/user/Documents/hands_on_II/ICB_class/"
setwd(workingdir)

dir.create("MDP")
setwd("MDP")

exp1<- expr
samplesinfo <- targets

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
exp1 <- exp1[ !duplicated( exp1$Symbol ), ]

table(is.na(exp1))
exp1 <- na.omit(exp1) 
rownames(exp1) <- as.character(exp1$Symbol)

# Get all sample without mean column (last column)
exp2 <- exp1[ , 2:( ncol( exp1 ) - 1 ) ]
NumberOfSamplesWithZeroAllowed <- round( 0.25 * ncol( exp2 ) )
exp2 <- exp2[rowSums( exp2 == 0 ) <= NumberOfSamplesWithZeroAllowed, ]




# Rodando MDP
BiocManager::install("mdp")
library(mdp)
mdp.results <- mdp( data = exp2, pdata = samplesinfo, control_lab = controlGroup,
                    file_name = paste( "GSE4806_MDP_", treatedGroup, "_", controlGroup, "_", sep = "" ) )






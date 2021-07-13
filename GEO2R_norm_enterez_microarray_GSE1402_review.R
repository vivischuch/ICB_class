#@vivischuch
#03/10/2018
#
library(RCurl)

y <- read.csv(text = getURL("https://raw.githubusercontent.com/vivischuch/classes/main/GSE23832.top.table.csv"), sep = "\t")
#
# install.packages("calibrate")
library(calibrate)
library(data.table)
library(GEOquery)
library(limma)
library(minfi)
library(RColorBrewer)
library(matrixStats)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(tidyverse)
library(data.table)
library(mdp)
library(reactome.db)
library(fgsea)

# ++++++++++

rm(list = ls())
stringsAsFactors=FALSE
workingdir <- "/home/user/Documents/hands_on_II/ICB/"
setwd(workingdir)
folder <- ""

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

library(GEOquery)

gse <- getGEO("GSE1402", GSEMatrix = TRUE)
eset <- gse[[1]] # Primeiro ExpressionSet
expr <- as.data.frame(Biobase::exprs(eset))
feature_data <- Biobase::fData(eset)
str(feature_data)
expr$Symbol <- feature_data$ENTREZ_GENE_ID
head(expr)

rownames(expr) <- NULL
head(expr)
expr$Symbol <- as.character((expr$Symbol))

pheno_data <- Biobase::pData(eset)
head(pheno_data)
write.table(pheno_data, "phenodata_GSE1402.tsv", sep='\t', row.names = F, quote = F)

# ++++++++++

dir.create(paste0(workingdir,"/MDP"))
write.table(expr, paste0(workingdir,"/MDP/","GSE1402_expr_forMDP_entrez.txt"), sep='\t', row.names = F, quote = F)

dirWant <- "/home/user/Documents/hands_on_II/ICB/MDP/"
expressionFile <- "GSE1402_expr_forMDP_entrez.txt"
annotationFile <- "/home/user/Documents/hands_on_II/ICB/phenodata_GSE1402_ok.tsv"

setwd(dirWant)

exp1<- read.table( expressionFile, header = T, na.strings = "NA", sep = "\t", row.names = NULL )
samplesinfo <- read.table( annotationFile, header = T, na.strings = "NA", sep = "\t" )

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

rownames( exp1 ) <- as.character( exp1$Symbol )

# Get all sample without mean column (last column)
exp2 <- exp1[ , 2:( ncol( exp1 ) - 1 ) ]
NumberOfSamplesWithZeroAllowed <- round( 0.25 * ncol( exp2 ) )
exp2 <- exp2[rowSums( exp2 == 0 ) <= NumberOfSamplesWithZeroAllowed, ]


# Rodando MDP
BiocManager::install("mdp")
library(mdp)
mdp.results <- mdp( data = exp2, pdata = samplesinfo, control_lab = controlGroup,
                    file_name = paste( "GSE4806_MDP_", treatedGroup, "_", controlGroup, "_", sep = "" ) )


# ++++++++++

dir.create(paste0(workingdir,"/number_samples"))
setwd(paste0(workingdir,"/number_samples"))

targets <- samplesinfo

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
# ++++ Fim FIG ++++

# ++++++++++

dir.create(paste0(workingdir,"/expression_files"))
setwd(paste0(workingdir,"/expression_files"))

expression <- exp2
expression <- as.matrix(expression)
# densityPlot(expression, xlab = "", main = "densityplot", las = 2, legend=FALSE)
boxplot(expression, main = "boxplot", border="#1B9E77", las = 2, cex.axis = 0.6)

# pal = brewer.pal(8,"Dark2")
# "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"
# ++++ FIG ++++
pdf(file = "GSE1402_expression_densityplot_boxplot.pdf", width=8,height=4)
# par(mfrow=c(1,2), cex.main=2)
# densityPlot(expression, xlab = "", main = "GSE1402 densityplot", las = 2, legend=FALSE)
boxplot(expression, main = "GSE1402 boxplot", border="#1B9E77", las = 2, cex.axis = 0.6)
dev.off()

png(file = "GSE1402_expression_densityplot_boxplot.png", units = "in", width = 8, height = 4, pointsize = 10, res = 100)
par(mfrow=c(1,2), cex.main=2)
densityPlot(expression, xlab = "", main = "GSE1402 densityplot", las = 2, legend=FALSE)
boxplot(expression, main = "GSE1402 boxplot", border="#1B9E77", las = 2, cex.axis = 0.6)
dev.off()
# ++++ Fim FIG ++++

Symbol2 <- rownames(exp3)
exprv2 <- cbind(Symbol, exp3)
rownames(exprv2) <- NULL
head(exprv2)
write.table(exprv2, "GSE1402_expression_MDP_enterez.tsv", sep='\t', row.names = F, quote = F)

summary(expression)
min(expression)
max(expression)

expr0 <- normalizeBetweenArrays(expression)

# ++++ FIG ++++
pdf(file = "GSE1402_normalized_expression_densityplot_boxplot.pdf", width=8,height=4)
par(mfrow=c(1,2), cex.main=2)
densityPlot(expr0, xlab = "", main = "GSE1402 densityplot", las = 2, legend=FALSE)
boxplot(expr0, main = "GSE1402 boxplot", border="#1B9E77", las = 2, cex.axis = 0.6)
dev.off()

png(file = "GSE1402_normalized_expression_densityplot_boxplot.png", units = "in", width = 8, height = 4, pointsize = 10, res = 100)
par(mfrow=c(1,2), cex.main=2)
densityPlot(expr0, xlab = "", main = "GSE1402 densityplot", las = 2, legend=FALSE)
boxplot(expr0, main = "GSE1402 boxplot", border="#1B9E77", las = 2, cex.axis = 0.6)
dev.off()
# ++++ Fim FIG ++++

min(expr0)
max(expr0)
summary(expr0)

ex <- log2(expr0)
min(ex)
max(ex)
summary(ex)

# ++++ FIG ++++
pdf(file = "GSE1402_log2_normalized_expression_densityplot_boxplot.pdf", width=8,height=4)
par(mfrow=c(1,2), cex.main=2)
densityPlot(ex, xlab = "", main = "GSE1402 densityplot", las = 2, legend=FALSE)
boxplot(ex, main = "GSE1402 boxplot", border="#1B9E77", las = 2, cex.axis = 0.6)
dev.off()

png(file = "GSE1402_log2_normalized_expression_densityplot_boxplot.png", units = "in", width = 8, height = 4, pointsize = 10, res = 100)
par(mfrow=c(1,2), cex.main=2)
densityPlot(ex, xlab = "", main = "GSE1402 densityplot", las = 2, legend=FALSE)
boxplot(ex, main = "GSE1402 boxplot", border="#1B9E77", las = 2, cex.axis = 0.6)
dev.off()
# ++++ Fim FIG ++++

exprzz <- 5.9 +(ex)
min(exprzz)
max(exprzz)
summary(exprzz)
dim(exprzz)

# ++++ FIG ++++
pdf(file = "GSE1402_5.9+(expr)_log2_normalized_expression_densityplot_boxplot.pdf", width=8,height=4)
par(mfrow=c(1,2), cex.main=2)
densityPlot(exprzz, xlab = "", main = "GSE1402 densityplot", las = 2, legend=FALSE)
boxplot(exprzz, main = "GSE1402 boxplot", border="#1B9E77", las = 2, cex.axis = 0.6)
dev.off()

png(file = "GSE1402_5.9+(expr)_log2_normalized_expression_densityplot_boxplot.png", units = "in", width = 8, height = 4, pointsize = 10, res = 100)
par(mfrow=c(1,2), cex.main=2)
densityPlot(exprzz, xlab = "", main = "GSE1402 densityplot", las = 2, legend=FALSE)
boxplot(exprzz, main = "GSE1402 boxplot", border="#1B9E77", las = 2, cex.axis = 0.6)
dev.off()
# ++++ Fim FIG ++++

Symbol <- rownames(exprzz)
exprv <- cbind(Symbol, exprzz)
rownames(exprv) <- NULL
head(exprv)
write.table(exprv, "GSE1402_5.9+(expr)_log2_normalized_expression_MDP_enterez.tsv", sep='\t', row.names = F, quote = F)

Symbol1 <- rownames(ex)
exprv1 <- cbind(Symbol, ex)
rownames(exprv1) <- NULL
head(exprv1)
write.table(exprv1, "GSE1402_log2_normalized_expression_MDP_enterez.tsv", sep='\t', row.names = F, quote = F)


# ++++++++++

dir.create(paste0(workingdir,"/t-test"))
setwd(paste0(workingdir,"/t-test"))

exp <- as.data.frame(expr0)

#Run t-test between Groups BEFORE MDP
Ps <- sapply(1:nrow(exp), function(x) { 
  if(all(is.na(exp[x,]))) NULL 
  else try(t.test(as.numeric(exp[x, which(as.character(samplesinfo$Class) == controlGroup) ]),
                  as.numeric(exp[x, which(as.character(samplesinfo$Class) == treatedGroup) ]))$p.value, 
           silent=TRUE) } )
Ps <- unlist(Ps)
AdjPs <- p.adjust(Ps, method = "BH", n = length(Ps))

exp$meanControl <- apply( exp[, which(as.character(samplesinfo$Class) == controlGroup) ], 1, mean )
exp$meanTreated <- apply( exp[, which(as.character(samplesinfo$Class) == treatedGroup) ], 1, mean )

exp$FC <- exp$meanTreated / exp$meanControl
exp$Log2FC <- log(exp$meanTreated / exp$meanControl,2)

FCs <- sapply(1:nrow(exp), function(x) { 
  if(all(is.na(exp[x,]))) NULL 
  else try(exp[x, "FC" ], 
           silent=TRUE) } )
FCs <- unlist(FCs)

Log2FCs <- sapply(1:nrow(exp), function(x) { 
  if(all(is.na(exp[x,]))) NULL 
  else try(exp[x, "Log2FC" ], 
           silent=TRUE) } )
Log2FCs <- unlist(Log2FCs)

DEG_table <- as.data.frame( cbind( Ps, AdjPs, FCs, Log2FCs ) )
rownames(DEG_table) <- rownames(exp)

AdjPcut <- 0.05
DEG_table1 <- DEG_table[DEG_table[,"AdjPs"] < AdjPcut,]
DEG_table2 <- DEG_table1[DEG_table1[,"Ps"] < AdjPcut,]
DEG_table3 <- DEG_table2[abs(DEG_table2[,"Log2FCs"]) > 1,]

DEG_table$Symbol <- rownames(DEG_table)
DEG_table <- DEG_table[order(DEG_table$Log2FCs),] 

genesFC <- as.vector <- DEG_table[, 4]
names( genesFC ) <- rownames( DEG_table )
genesFC <- genesFC[ order( genesFC ) ]

dim(DEG_table)
summary(DEG_table)

DEG_table$minuslog10P <- -log10(DEG_table$Ps)

# ++++ FIG ++++
pdf(file = "GSE1402_volcanoplot.pdf", width=8,height=4)
plot(DEG_table$Log2FCs, DEG_table$minuslog10P, main="GSE1402 Volcano Plot", xlab="log2FC", ylab="-log10P", las=1, pch=20)
abline(h=1.3, col="blue") # add horizontal straight line at the cut-off
abline(v=c(-1,1), col="blue") # add vertical straight lines
# Add grey points to unperturbed genes
with(subset(DEG_table, minuslog10P < 1.3), points(Log2FCs, minuslog10P, pch=20, col="gray")) 
# Add green points to downregulated genes
with(subset(DEG_table, Log2FCs < -1 & minuslog10P > 1.3), points(Log2FCs, minuslog10P, pch=20, col="blue")) 
# Add red points to upregulated genes
with(subset(DEG_table, Log2FCs > 1 & minuslog10P > 1.3), points(Log2FCs, minuslog10P, pch=20, col="red")) 
#label down genes
with(subset(DEG_table, Log2FCs < -2 & minuslog10P >1.3), textxy(Log2FCs, minuslog10P, labs=Symbol, cex=.8, col="blue")) 
#label up genes
with(subset(DEG_table, Log2FCs > 2 & minuslog10P >1.3), textxy(Log2FCs, minuslog10P, labs=Symbol, cex=.8, col="red")) 
dev.off()


png(file = "GSE1402_volcanoplot.png", units = "in", width = 8, height = 4, pointsize = 10, res = 100)
plot(DEG_table$Log2FCs, DEG_table$minuslog10P, main="GSE1402 Volcano Plot", xlab="log2FC", ylab="-log10P", las=1, pch=20)
abline(h=1.3, col="blue") # add horizontal straight line at the cut-off
abline(v=c(-1,1), col="blue") # add vertical straight lines
# Add grey points to unperturbed genes
with(subset(DEG_table, minuslog10P < 1.3), points(Log2FCs, minuslog10P, pch=20, col="gray")) 
# Add green points to downregulated genes
with(subset(DEG_table, Log2FCs < -1 & minuslog10P > 1.3), points(Log2FCs, minuslog10P, pch=20, col="blue")) 
# Add red points to upregulated genes
with(subset(DEG_table, Log2FCs > 1 & minuslog10P > 1.3), points(Log2FCs, minuslog10P, pch=20, col="red")) 
#label down genes
with(subset(DEG_table, Log2FCs < -2 & minuslog10P >1.3), textxy(Log2FCs, minuslog10P, labs=Symbol, cex=.8, col="blue")) 
#label up genes
with(subset(DEG_table, Log2FCs > 2 & minuslog10P >1.3), textxy(Log2FCs, minuslog10P, labs=Symbol, cex=.8, col="red")) 
dev.off()
# ++++ Fim FIG ++++

Symbol3 <- rownames(DEG_table)
exprv3 <- cbind(Symbol3, DEG_table)
rownames(exprv3) <- NULL
head(exprv3)
write.table(exprv3, "GSE1402_DEG_table_enterez.tsv", sep='\t', row.names = F, quote = F)

# ++++++++++

dir.create(paste0(workingdir,"/fgsea"))
setwd(paste0(workingdir,"/fgsea"))

ranks <- genesFC

head(ranks)
class(ranks)
str(ranks)

# ++++ FIG ++++
pdf(file = "GSE1402_fGSEA_barplot_rank.pdf", width=8,height=4)
barplot(sort(ranks, decreasing = T), main = "GSE1402 barplot rank")
dev.off()

png(file = "GSE1402_fGSEA_barplot_rank.png", units = "in", width = 8, height = 4, pointsize = 10, res = 100)
barplot(sort(ranks, decreasing = T), main = "GSE1402 barplot rank")
dev.off()
# ++++ Fim FIG ++++

library(reactome.db)
pathways <- reactomePathways(names(ranks))
str(pathways)

fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)

# fwrite(fgseaRes, file="GSE51402_fgseaRes.txt", sep="\t", sep2=c("", " ", ""))

head(fgseaRes[order(padj, -abs(NES)), ], n=10)
head(fgseaRes[order(pval), ], 10)

sum(fgseaRes[, padj < 0.05])
sum(fgseaRes[, pval < 0.05])

# ++++ FIG ++++
pdf(file = "GSE1402_fGSEA_plotEnrichment.pdf", width=8,height=4)
plotEnrichment(pathways, ranks)
dev.off()

png(file = "GSE1402_fGSEA_plotEnrichment.png", units = "in", width = 8, height = 4, pointsize = 10, res = 100)
plotEnrichment(pathways, ranks)
dev.off()
# ++++ Fim FIG ++++

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fwrite(fgseaResTidy, file="GSE1402_fGSEA_fgseaResTidy.txt", sep="\t", sep2=c("", " ", ""))

# ++++ FIG ++++
pdf(file = "GSE1402_fGSEA_Pathways_NES_padj.pdf", width=8,height=4)
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  theme(axis.text.y = element_blank())+
  labs(x="Pathway", y="NES",
       title=" GSE1402 Pathways NES from GSEA")
dev.off()

png(file = "GSE1402_fGSEA_Pathways_NES_padj.png", units = "in", width = 8, height = 4, pointsize = 10, res = 100)
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  theme(axis.text.y = element_blank())+
  labs(x="Pathway", y="NES",
       title=" GSE1402 Pathways NES from GSEA")
dev.off()
# ++++ Fim FIG ++++


topcfs <- fgseaResTidy
topcfs <- subset(topcfs, padj < 0.05)
fwrite(topcfs, file="GSE1402_fGSEA_subset_padj_0.05.txt", sep="\t", sep2=c("", " ", ""))

topUp <- fgseaRes %>% 
  filter(ES > 0) %>% 
  top_n(15, wt=-pval)
topDown <- fgseaRes %>% 
  filter(ES < 0) %>% 
  top_n(15, wt=-pval)
topPathways <- bind_rows(topUp, topDown) %>% 
  arrange(-NES)

fwrite(topPathways, file="GSE1402_fGSEA_subset_pval_NES.txt", sep="\t", sep2=c("", " ", ""))

data <- topPathways

# ++++ FIG ++++
pdf(file = "GSE1402_fGSEA_ggplot_subset_pval_NES.pdf", width=12,height=8)
p <- ggplot(data, aes(NES, pathway))

p + geom_point(aes(colour=pval, size=size)) +
  scale_color_gradientn(colours=rainbow(4)) +
  geom_vline(xintercept=0, size=0.5, colour="gray50") +
  theme(panel.background=element_rect(fill="gray95", colour="gray95"),
        panel.grid.major=element_line(size=0.25,linetype='solid', colour="gray90"), 
        panel.grid.minor=element_line(size=0.25,linetype='solid', colour="gray90"),
        axis.title.y=element_blank()) +
  expand_limits(x=c(-3,3)) +
  scale_x_continuous(breaks=c(-3,-2,-1,0,1,2,3)) +
  scale_y_discrete(limits=rev(data$pathway)) +
  labs(title="GSE1402 Pathways")
dev.off()

png(file = "GSE1402_fGSEA_ggplot_subset_pval_NES.pdf.png",units = "in", width = 12, height = 8, pointsize = 10, res = 100)
p <- ggplot(data, aes(NES, pathway))

p + geom_point(aes(colour=pval, size=size)) +
  scale_color_gradientn(colours=rainbow(4)) +
  geom_vline(xintercept=0, size=0.5, colour="gray50") +
  theme(panel.background=element_rect(fill="gray95", colour="gray95"),
        panel.grid.major=element_line(size=0.25,linetype='solid', colour="gray90"), 
        panel.grid.minor=element_line(size=0.25,linetype='solid', colour="gray90"),
        axis.title.y=element_blank()) +
  expand_limits(x=c(-3,3)) +
  scale_x_continuous(breaks=c(-3,-2,-1,0,1,2,3)) +
  scale_y_discrete(limits=rev(data$pathway)) +
  labs(title="GSE1402 Pathways")
dev.off()
# ++++ Fim FIG ++++






#Run t-test between Groups AFTER MDP
outliersHC <- c("control1", "control18", "control12", "control13",
                "SLE3", "SLE12", "SLE15", "SLE16", "SLE17", "SLE20", "SLE29", "SLE77")
#outliersHC <- c("control1", "control18", "control12", "control13")
#outliersT <- c("SLE3", "SLE12", "SLE15", "SLE16", "SLE17", "SLE20", "SLE29", "SLE77")

exp2 <- exp[,!(names(exp) %in% outliersHC) ]
samplesinfo2 <- samplesinfo[!(samplesinfo$Sample %in% outliersHC),]

#exp2 <- exp[,!(names(exp) %in% outliersT) ]
#samplesinfo2 <- samplesinfo[!(samplesinfo$Sample %in% outliersT),]

Ps2 <- sapply(1:nrow(exp2), function(x) { 
  if(all(is.na(exp2[x,]))) NULL 
  else try(t.test(as.numeric(exp2[x, which(as.character(samplesinfo2$Class) == controlGroup) ]),
                  as.numeric(exp2[x, which(as.character(samplesinfo2$Class) == treatedGroup) ]))$p.value, 
           silent=TRUE) } )
Ps2 <- unlist(Ps2)
AdjPs2 <- p.adjust(Ps2, method = "BH", n = length(Ps2))


exp2$meanControl <- apply( exp2[, which(as.character(samplesinfo2$Class) == controlGroup) ], 1, mean )
exp2$meanTreated <- apply( exp2[, which(as.character(samplesinfo2$Class) == treatedGroup) ], 1, mean )

exp2$FC <- exp2$meanTreated / exp2$meanControl
exp2$Log2FC <- log(exp2$meanTreated / exp2$meanControl,2)

FCs2 <- sapply(1:nrow(exp2), function(x) { 
  if(all(is.na(exp2[x,]))) NULL 
  else try(exp2[x, "FC" ], 
           silent=TRUE) } )
FCs2 <- unlist(FCs2)

Log2FCs2 <- sapply(1:nrow(exp2), function(x) { 
  if(all(is.na(exp2[x,]))) NULL 
  else try(exp2[x, "Log2FC" ], 
           silent=TRUE) } )
Log2FCs2 <- unlist(Log2FCs2)


DEG_table2 <- as.data.frame( cbind( Ps2, AdjPs2, FCs2, Log2FCs2 ) )
rownames(DEG_table2) <- rownames(exp2)
DEG_table2 <- DEG_table2[DEG_table2[,"AdjPs2"] < AdjPcut,]
#DEG_table <- DEG_table[DEG_table[,"Ps"] < AdjPcut,]
#DEG_table <- DEG_table[abs(DEG_table[,"FCs"]) > 0.1,]


DEG_table2$Symbol <- rownames(DEG_table2)
DEG_table2 <- DEG_table2[order(DEG_table2$Log2FCs),] 

genesFC2 <- as.vector <- DEG_table2[, 3]
names( genesFC2 ) <- rownames( DEG_table2 )
genesFC2 <- genesFC2[ order( genesFC2 ) ]

dim(DEG_table)
dim(DEG_table2)


#fgseaResultAFTERRemoveOutlier <- run_fgsea( "~/../Downloads/c2.cp.reactome.v6.1.symbols.gmt", ranksOfGenes = genesFC2, minSizeGroup = 15, maxSizeGroup = 500, npermGroup = 10000, filterPathways = TRUE, pFilter = AdjPcut)



mdp.results <- mdp( data = exp2, pdata = samplesinfo2, control_lab = controlGroup,
                    file_name = paste( "MDP_AFTER_OUTLIERS_", treatedGroup, "_", controlGroup, "_", sep = "" ) )


##Remove the same number of samples, randomly and see how many DEGs



source("https://raw.githubusercontent.com/nicolau/code-R/master/funcoes_para_diagrama_venn.R")
plot.pairwise.venn(a1 = fgseaResultAFTERRemoveOutlier$pathway, a2 = fgseaResultBEFORERemoveOutlier$pathway, labels = c("AFTER", "BEFORE"), file = "GSE72509_venn_diagram.pdf", graphicType = "pdf")





dim(DEG_table)
dim(DEG_table2)


















#### Teste para criacao do metodo estatistico para remocao de outliers
xControl <- mdp.results$sample_scores$allgenes$Score[ mdp.results$sample_scores$allgenes$Class == controlGroup ]
meanControl <- mean( mdp.results$sample_scores$allgenes$Score[ mdp.results$sample_scores$allgenes$Class == controlGroup ] )
sdControl <- sd( mdp.results$sample_scores$allgenes$Score[ mdp.results$sample_scores$allgenes$Class == controlGroup ] )

(xControl - meanControl) / sdControl


xTreated <- mdp.results$sample_scores$allgenes$Score[ mdp.results$sample_scores$allgenes$Class == treatedGroup ]
meanTreated <- mean( mdp.results$sample_scores$allgenes$Score[ mdp.results$sample_scores$allgenes$Class == treatedGroup ] )
sdTreated <- sd( mdp.results$sample_scores$allgenes$Score[ mdp.results$sample_scores$allgenes$Class == treatedGroup ] )

hist((xTreated - meanTreated) / sdTreated)

sort(xControl)
min(xTreated)
max(xControl)

mdp.results$sample_scores$allgenes


newdata <- mtcars[order(mpg, -cyl),]
mdp.results$sample_scores$allgenes[order(mdp.results$sample_scores$allgenes$Score),]



sample_scores_list <- mdp.results$sample_scores
# select sample scores calculated using the perturbed genes
sample_scores <- sample_scores_list[["perturbedgenes"]]
head(sample_scores)
##### FInal do teste









































#### Comparar as amostras outlier do controle em rela??o ao controle remanescentes ####
expressionFile <- "GSE72509_SLE_RPKMs_expression.tsv"
annotationFile <- "GSE72509_SLE_RPKMs_pheno_outlier_comparison.tsv"
controlGroup <- "control"
treatedGroup <- "controlOutlier"

controlGroup <- "SLEOutlier"
treatedGroup <- "SLE"

AdjPcut <- 0.05

library(mdp)
setwd(dirWant)

exp         <- read.table( expressionFile, header = T, na.strings = "NA", sep = "\t", row.names = NULL )
samplesinfo <- read.table( annotationFile, header = T, na.strings = "NA", sep = "\t" )

#Remove samples that are not in comparison
samplesinfo <- samplesinfo[ c( which( samplesinfo$Class == controlGroup ),
                               which( samplesinfo$Class == treatedGroup ) ), ]

samplesinfo <- samplesinfo[ samplesinfo$Sample %in% colnames( exp ), ]

exp <- exp[ ,c( "Symbol", as.character( samplesinfo$Sample ) ) ]
samplesinfo <- samplesinfo[ samplesinfo$Sample %in% colnames( exp ), ]


#Colapse symbols that are duplicated by taking the one with higest expression
exp$meanG <- apply( exp[ , 2:ncol( exp ) ], 1, mean )
exp <- exp[ order( exp[ , 'Symbol' ],
                   exp[ , 'meanG' ] ), ]
exp <- exp[ !duplicated( exp$Symbol ), ]

rownames( exp ) <- as.character( exp$Symbol )
# Get all sample without mean column (last column)
exp2 <- exp[ , 2:( ncol( exp ) - 1 ) ]
NumberOfSamplesWithZeroAllowed <- round( 0.25 * ncol( exp2 ) )
exp2 <- exp2[rowSums( exp2 == 0 ) <= NumberOfSamplesWithZeroAllowed, ]
exp <- exp2

#Exporte o novo (e final) expression
#write.table(exp,file="GSE72509_SLE_RPKMs_expression_final.tsv", col.names = NA, sep="\t",row.names = TRUE)


# Rodando MDP antes de remover os adaptadores
mdp.results <- mdp( data = exp, pdata = samplesinfo, control_lab = controlGroup,
                    file_name = paste( "MDP_", treatedGroup, "_", controlGroup, "_", sep = "" ))#,
                    #pathways = fgsea::gmtPathways("~/../Downloads/c2.cp.reactome.v6.0.symbols.gmt"))


#Run t-test between Groups BEFORE MDP
Ps <- sapply(1:nrow(exp), function(x) { 
  if(all(is.na(exp[x,]))) NULL 
  else try(t.test(as.numeric(exp[x, which(as.character(samplesinfo$Class) == controlGroup) ]),
                  as.numeric(exp[x, which(as.character(samplesinfo$Class) == treatedGroup) ]))$p.value, 
           silent=TRUE) } )
Ps <- unlist(Ps)
AdjPs <- p.adjust(Ps, method = "BH", n = length(Ps))

exp$meanControl <- apply( exp[, which(as.character(samplesinfo$Class) == controlGroup) ], 1, mean )
exp$meanTreated <- apply( exp[, which(as.character(samplesinfo$Class) == treatedGroup) ], 1, mean )

exp$FC <- exp$meanTreated / exp$meanControl
exp$Log2FC <- log(exp$meanTreated / exp$meanControl,2)

FCs <- sapply(1:nrow(exp), function(x) { 
  if(all(is.na(exp[x,]))) NULL 
  else try(exp[x, "FC" ], 
           silent=TRUE) } )
FCs <- unlist(FCs)

Log2FCs <- sapply(1:nrow(exp), function(x) { 
  if(all(is.na(exp[x,]))) NULL 
  else try(exp[x, "Log2FC" ], 
           silent=TRUE) } )
Log2FCs <- unlist(Log2FCs)


DEG_table <- as.data.frame( cbind( Ps, AdjPs, FCs, Log2FCs ) )
rownames(DEG_table) <- rownames(exp)
DEG_table <- DEG_table[DEG_table[,"AdjPs"] < AdjPcut,]
#View(DEG_table)
#DEG_table <- DEG_table[DEG_table[,"Ps"] < AdjPcut,]
#DEG_table <- DEG_table[abs(DEG_table[,"FCs"]) > 0.1,]


DEG_table$Symbol <- rownames(DEG_table)
DEG_table <- DEG_table[order(DEG_table$Log2FCs),] 

genesFC <- as.vector <- DEG_table[, 3]
names( genesFC ) <- rownames( DEG_table )
genesFC <- genesFC[ order( genesFC ) ]


dim(DEG_table)

































### Calcula o arranjo combinatorio para evitar os mesmos outliers e consequentemente o mesmo numero de DEGs

print("tste")
set.seed(2018)
write.table( x = "", file = "random_outlier_removing_DEGs.txt", append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE)
for(i in 1:1000) {
  print(paste("############ Interaction ", i, "############", sep = "" ))
  #### Random outlier removing####
  numberOfHealthySamples <- 4
  numberOfTreatedSamples <- 8
  
  print("Random outlier samples from healthy and treated conditions...")
  outlierHealthy <- as.vector(sample(x = samplesinfo[ which( samplesinfo$Class == controlGroup ), 1], numberOfHealthySamples, replace = TRUE))
  outlierTreated <- as.vector(sample(x = samplesinfo[ which( samplesinfo$Class == treatedGroup ), 1], numberOfTreatedSamples, replace = TRUE))
  outlierJoined <- c(outlierHealthy, outlierTreated)

  exp2 <- exp[,!(names(exp) %in% outlierJoined) ]
  samplesinfo2 <- samplesinfo[!(samplesinfo$Sample %in% outlierJoined),]
  
  print("Doing the statistical analysis...")
  Ps2 <- sapply(1:nrow(exp2), function(x) { 
    if(all(is.na(exp2[x,]))) NULL 
    else try(t.test(as.numeric(exp2[x, which(as.character(samplesinfo2$Class) == controlGroup) ]),
                    as.numeric(exp2[x, which(as.character(samplesinfo2$Class) == treatedGroup) ]))$p.value, 
             silent=TRUE) } )
  Ps2 <- unlist(Ps2)
  AdjPs2 <- p.adjust(Ps2, method = "BH", n = length(Ps2))
  
  DEG_table2 <- as.data.frame( cbind( Ps2, AdjPs2 ) )
  rownames(DEG_table2) <- rownames(exp2)
  DEG_table2 <- DEG_table2[DEG_table2[,"AdjPs2"] < AdjPcut,]

  print("Saving results in file...")
  write.table( x = paste( "interaction: ", i, " - DEGs: ", dim(DEG_table2)[1], " - OutlierSamples: ", paste( outlierJoined, collapse = "|" ), sep = "" ), file = "random_outlier_removing_DEGs.txt", append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)

  print("############ Done ############")
  print("##############################")
}

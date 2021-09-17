####################################################
#####################################################
#### Paper Name: Disrupted PGR-B and ESR1 signaling underlies defective decidualization linked to severe preeclampsia
#### eLife, 2021


# Data Processing ---------------------------------------------------------
# - Data (fastq files) can be downloaded from GEO: Series GSE172381
# - Preprocessinfg data can be found at script 01_preProcessingSamples.sh
# - After the script "01", you obtain sample_htseq_counts.txt


# 1. Obtain raw counts from htseq files -----------------------------------

# a) This scripts requires the following paramaters as input:
# analysisPath--> Folder where fastq files for all samples are saved 
#  analysisPath/sample/counts/sample_htseq_counts.txt
#
# b) Or you can load directly from GEO/github raw counts, if you want to avoid step a)


samples<- dir(analysisPath)

cov <- sapply(samples, function(x) {
  sample_cov_file <- file.path(analysisPath, x, "counts", as.character(list.files(file.path(analysisPath, x, "counts"), pattern = "_htseq_counts.txt")[1]))
  res <- get.coverage.htseqcounts(sample_cov_file)
  res$sample <- x
  return(list(res))
})
cov_df <- do.call(rbind, cov)
raw_counts.temp <- acast(cov_df, gene ~ sample, value.var="nreads")

non_genes <- c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
raw_counts<- raw_counts.temp[!rownames(raw_counts.temp) %in% non_genes,]

dim(raw_counts) #  56638    40

final_gcts<- raw_counts


# 2. Load necessary objects -----------------------------------------

## Functions and libraries

source(file=paste0(pipeRoot, "functionsPREECLAMPSIA.R"))


## This script requires the following files:
#  ed.rev.cli: Experimental design you can load from github
#  final_gcts: Raw counts, obtained previously or you can load from GEO/github

ed.rev.cli<- read.csv(file=paste0(pipeRoot, "edRevPaper.csv"))
final_gcts<- read.csv(file=paste0(pipeRoot, "rawData.csv"), row.names = 1,check.names = F)
dim(final_gcts) ## 56638    40


## Colors of the figures

colores<- c("#1f78b4","#ff7f00") 
clasesCol<- c("Control", "sPE")
names(colores)<- clasesCol


## Table with genes with Entrez Id

genesEntrezId<- read.table(file=paste0(pipeRoot, "genesEntrezId_paper.txt"))



# 3. Exploratory analysis: All Samples ---------------------------------------------

## Normalization all samples

final_norm <- get.normalization(final_gcts, ed.rev.cli, round(mean(table(ed.rev.cli$group)))/2)
names(final_norm) # "counts"  "samples"


## Dispersion

final_norm_disp <- estimate.disp.pairw(final_norm) 
final_norm_disp_cpm <- cpm(final_norm_disp) 
final_logcpm <- cpm(final_norm, prior.count=2, log=TRUE)
dim(final_logcpm) # 18301    40

## PCA

final_fit_norm <- prcomp(t(na.omit(final_logcpm)))


#  3.1 Balanced split data: TRAIN y TEST ----------------------------------------------------

table(ed.rev.cli[ed.rev.cli$setNew=="train", "group"])
# Control     sPE 
# 12      17 

table(ed.rev.cli[ed.rev.cli$setNew=="test", "group"])
# Control     sPE 
# 4       7 

### PCA plot with Train y Test set

final_pci_norm1 <- data.frame(final_fit_norm$x, 
                              color_group=ed.rev.cli[match(colnames(final_logcpm), ed.rev.cli$sample), "group"],
                              shape_set=ed.rev.cli[match(colnames(final_logcpm), ed.rev.cli$sample), "setNew"],
                              paperName= ed.rev.cli[match(colnames(final_logcpm), ed.rev.cli$sample), "paperName"])


PCA.set.norm<- ggplot(final_pci_norm1, aes(x=PC1, y=PC2, color=color_group, shape=shape_set)) + 
  geom_point(size=3) +
  scale_x_continuous(name=paste("PC1 (", abs(round(summary(final_fit_norm)$importance[2, "PC1"]*100)), "%)", sep="")) +
  scale_y_continuous(name=paste("PC2 (", abs(round(summary(final_fit_norm)$importance[2, "PC2"]*100)), "%)", sep="")) +
  labs(title="PCA norm counts")+
  scale_color_manual(name = "Group", values = colores) +
  theme_light() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust=0.5))+
  geom_text_repel(aes(x=PC1, y=PC2, color=color_group, label = paperName))


# 3.2 Preparing intermediate files for Train and Test data set ---------------------------------------

ed.train<- ed.rev.cli %>%
  filter(setNew=="train")

ed.train$group<- factor(ed.train$group)

ed.test<- ed.rev.cli %>%
  filter(setNew=="test")

ed.test$group<- factor(ed.test$group)

final_gcts.train <- final_gcts[, colnames(final_gcts) %in% ed.train$sample]
final_gcts.test <- final_gcts[, colnames(final_gcts) %in% ed.test$sample]

rownames(ed.train)<- ed.train$sample
rownames(ed.test)<- ed.test$sample

ed.train<- ed.train[colnames(final_gcts.train),]
ed.test<- ed.test[colnames(final_gcts.test),]

# 4.Normalization Train data set ----------------------------------------------

group<- as.factor(ed.train$group)
table(group)

y.train<- DGEList(counts=final_gcts.train, group=group)

## Normalization

min_samples <- round(mean(table(ed.train$group)))/2 

keep <- rowSums(cpm(y.train) > 1) >= min_samples 

y.train <- y.train[keep, keep.lib.sizes=FALSE]
y.train.norm <- calcNormFactors(y.train, method="TMM") 


# 5. Testing for DE genes --------------------------------------------------------------

# 5.1 Design Matrix -------------------------------------------------------

design<- model.matrix(~0+group) 
rownames(design)<- colnames(y.train.norm)

# 5.2 Estimating Dispersion -----------------------------------------------

## Estimating the dispersion

y.train.norm<-estimateDisp(y.train.norm, design, robust=TRUE)

# 5.3 Data Exploration  ---------------------------------

## Multi-dimensional scaling (MDS) plot

points <- c(16,16)
colors <- rep(c("blue", "orange"), 2)

plotMDS(y.train.norm, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)


## CPM Normalized Data:Train

final_logcpm.train <- cpm(y.train.norm, prior.count=2, log=TRUE)
final_logcpm.train<- as.data.frame(final_logcpm.train)

## ## CPM Normalized Data:Test 

final_logcpm.test <- cpm(final_gcts.test, prior.count=2, log=TRUE)

## PCA plot

final_qc_norm <- get.df2plot(data.frame(final_logcpm.train, check.names = F), ed.train)
final_qc_norm$group <- factor(final_qc_norm$group, levels=c("Control","sPE"))
final_qc_norm$sample <- factor(final_qc_norm$sample, levels=unique(final_qc_norm$sample[order(final_qc_norm$group, final_qc_norm$run_id)] ), ordered=TRUE)

final_fit_norm <- prcomp(t(na.omit(final_logcpm.train)))
final_pci_norm <- data.frame(final_fit_norm$x, 
                             color_group=ed.train[match(colnames(final_logcpm.train), ed.train$sample), "group"],
                             paperName= ed.train[match(colnames(final_logcpm.train), ed.train$sample), "paperName"])

PCA12.norm.train<- ggplot(final_pci_norm, aes(x=PC1, y=PC2, color=color_group)) + 
  geom_point(size=3) +
  scale_x_continuous(name=paste("PC1 (", abs(round(summary(final_fit_norm)$importance[2, "PC1"]*100)), "%)", sep="")) +
  scale_y_continuous(name=paste("PC2 (", abs(round(summary(final_fit_norm)$importance[2, "PC2"]*100)), "%)", sep="")) +
  labs(title="PCA norm counts - Train Set")+
  scale_color_manual(name = "Group", values = colores) +
  theme_light() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust=0.5))+
  geom_text_repel(aes(x=PC1, y=PC2, color=color_group, label = paperName))

# 5.4 Fit model: GLM ---------------------------------------------------------------

fit<-glmQLFit(y.train.norm, design,robust=TRUE)

## Contrast to test

con<- makeContrasts(groupsPE-groupControl, levels = design)
qlf<- glmQLFTest(fit, contrast=con)


# 6. TREAT model ---------------------------------------------------------

# Given a fold-change threshold, the TREAT method is applied using the DGEGLM object produced previously.

############################################
### fc=1.2
############################################

tr.FC12 <- glmTreat(fit, contrast=con, lfc=log2(1.2))

sPEctrl.treat.FC12<- topTags(tr.FC12, n=nrow(y.train.norm), p.value=1, adjust.method = "BH" )
sPEctrl.treat.FC12.deg<- sPEctrl.treat.FC12$table

# filter significant p-value
# FDR: False Discovery Rate

sPEctrl.treat.FC12.deg<- sPEctrl.treat.FC12.deg %>%
  transform(hgnc = rownames(sPEctrl.treat.FC12.deg)) %>%
  filter(FDR<=0.05)

topTags(tr.FC12)
CvsP.tab.F12<-summary(decideTests(tr.FC12))

## Volcano plot with DE genes

dataPlot<- sPEctrl.treat.FC12$table %>%
  mutate(threshold=ifelse(FDR<=0.05 & logFC>=log2(1.2), "UP",
                          ifelse(FDR<=0.05 & logFC<log(1.2), "DOWN", "Non-Significant")))


VP.FC12<- ggplot(data=dataPlot, aes(x=logFC, y=-log10(FDR))) + 
  geom_point(aes(colour = threshold),alpha=0.8, size=2.5) + 
  geom_hline(yintercept = -log10(0.05), linetype="solid", colour="darkgrey") + 
  theme_bw() + theme(legend.position = "right") + 
  scale_colour_manual(values = c("UP"= "red", "DOWN"= "blue","Non-Significant"= "grey")) + 
  xlab("log2 fold change") + 
  ylab("-log10 FDR") +
  labs(title="Control vs sPE: log2(FC=1.2) & FDR< 0.05 \n  ")


## Subset of genes with Entrez Id asociated

sPEctrl.treat.FC12.deg.EG<- merge(x=sPEctrl.treat.FC12.deg, y=genesEntrezId, by.x="hgnc", by.y = "hgnc", all.y = F)
sPEctrl.treat.FC12.deg.EG<- subset(sPEctrl.treat.FC12.deg.EG, !is.na(EntrezGene)) ## 509 genes


################################################
### fc=1.4
################################################


tr.FC14 <- glmTreat(fit, contrast=con, lfc=log2(1.4))

CvsP.tab.F14<- summary(decideTests(tr.FC14))

sPEctrl.treat.FC14<- topTags(tr.FC14, n=nrow(y.train.norm), p.value=1, adjust.method = "BH" )
sPEctrl.treat.FC14.deg<- sPEctrl.treat.FC14$table

sPEctrl.treat.FC14.deg<- sPEctrl.treat.FC14.deg %>%
  transform(hgnc = rownames(sPEctrl.treat.FC14.deg)) %>%
  filter(FDR<=0.05)



sPEctrl.treat.FC14.deg.EG<- merge(x=sPEctrl.treat.FC14.deg, y=genesEntrezId, by.x="hgnc", by.y = "hgnc", all.y = F)
sPEctrl.treat.FC14.deg.EG<- subset(sPEctrl.treat.FC14.deg.EG, !is.na(EntrezGene)) ## 120 genes

## Volcano plot with some important genes


mygenes<- c("MMP11","IHH","IL6","TNF","FOSL1",
            "FILIP1","GRAMD1C","LIPI","LPIN3","MTND1P23")

dataPlot.FC14<- sPEctrl.treat.FC14$table %>%
  mutate(hgnc=rownames(sPEctrl.treat.FC14$table))

dataPlot.FC14.Entrez<- merge(dataPlot.FC14, genesEntrezId, by="hgnc")

dataPlot.FC14.Entrez<- dataPlot.FC14.Entrez %>%
  filter(!is.na(EntrezGene)) %>%
  mutate( threshold=ifelse(FDR<=0.05 & logFC>=log2(1.4), "UP",
                           ifelse(FDR<=0.05 & logFC<log(1.4), "DOWN", "Non-Significant")))

VP.FC14.genes<- ggplot(data=dataPlot.FC14.Entrez, aes(x=logFC, y=-log10(FDR))) + 
  geom_point(aes(colour = threshold),alpha=0.8, size=2.5) + 
  geom_hline(yintercept = -log10(0.05), linetype="solid", colour="darkgrey") + 
  theme_bw() + theme(legend.position = "right") + 
  scale_colour_manual(values = c("UP"= "red", "DOWN"= "blue","Non-Significant"= "grey")) + 
  xlab("log2 fold change") + 
  ylab("-log10 FDR") +
  labs(title="Control vs sPE: log2(FC=1.4) & FDR<0.05 \n  ") +
  annotate(geom="text", x= -4.9, y= 2.3, label="106 DEGs", color="blue", size=5.5)+
  annotate(geom="text", x= 4, y= 2.3, label="14 DEGs", color="red", size=5.5)+
  geom_text_repel(
    data=subset(dataPlot.FC14, hgnc %in% mygenes),
    aes(label=hgnc),
    size = 4,
    box.padding = unit(0.7, "lines"),
    point.padding = unit(0.5, "lines")
  )


# 6.1 Heatmap and PCA plot with DE genes ----------------------------------------------------

## TRAIN

dataNorm.CvsP.FC12<-final_logcpm.train[rownames(final_logcpm.train) %in% sPEctrl.treat.FC12.deg.EG$hgnc,]
dataNorm.CvsP.FC14<-final_logcpm.train[rownames(final_logcpm.train) %in% sPEctrl.treat.FC14.deg.EG$hgnc,]


final_ed1<- ed.train
final_ed1<- final_ed1[colnames(final_logcpm.train),]

colnames(dataNorm.CvsP.FC12)<- final_ed1$paperName
colnames(dataNorm.CvsP.FC14)<- final_ed1$paperName

rownames(final_ed1)<- final_ed1$paperName

dist<- "canberra" ## distance used for plotting Heatmap

plotHM.CvsP.FC12.dist<- makeHeatMap.dist(dataNormDE=dataNorm.CvsP.FC12, ed=final_ed1, group="group",  
                                         mainTit=paste0("log2(FC=1.2), Control vs sPE - ",  nrow(dataNorm.CvsP.FC12), " genes ", dist,  " dist."),
                                         colores = colores, sampleId = "paperName", distance = dist) 


plotHM.CvsP.FC14.dist<- makeHeatMap.dist(dataNormDE=dataNorm.CvsP.FC14, ed=final_ed1, group="group",  
                                         mainTit=paste0("log2(FC=1.4), Control vs sPE - ",  nrow(dataNorm.CvsP.FC14), " genes ", dist,  " dist."),
                                         colores = colores, sampleId = "paperName", distance = dist) 




## PCA

pca.CvsP.FC12<- pcaShape(dataNormDE=dataNorm.CvsP.FC12, ed=final_ed1, group="group",
                         mainTit=paste0("log2(FC=1.2), Control vs sPE - ", nrow(dataNorm.CvsP.FC12), " genes"),
                         shape="none", col=colores)
pca.CvsP.FC12$plotPCA


pca.CvsP.FC14<- pcaShape(dataNormDE=dataNorm.CvsP.FC14, ed=final_ed1, group="group",
                         mainTit=paste0("log2(FC=1.4), Control vs sPE - ", nrow(dataNorm.CvsP.FC14), " genes"),
                         shape="none", col=colores)
pca.CvsP.FC14$plotPCA


### TEST SET

dataNorm.CvsP.FC12.test<-final_logcpm.test[rownames(final_logcpm.test) %in%  sPEctrl.treat.FC12.deg.EG$hgnc,]
dataNorm.CvsP.FC14.test<-final_logcpm.test[rownames(final_logcpm.test) %in%  sPEctrl.treat.FC14.deg.EG$hgnc,]



final_ed2<- ed.test
final_ed2<- final_ed2[colnames(final_logcpm.test),]

colnames(dataNorm.CvsP.FC12.test)<- final_ed2$paperName
colnames(dataNorm.CvsP.FC14.test)<- final_ed2$paperName

rownames(final_ed2)<- final_ed2$paperName

## Heatmap

plotHM.CvsP.FC12.dist.test<- makeHeatMap.dist(dataNormDE=dataNorm.CvsP.FC12.test, ed=final_ed2, group="group",  
                                              mainTit=paste0("log2(FC=1.2), Control vs sPE - ",  nrow(dataNorm.CvsP.FC12), " genes ", dist,  " dist."),
                                              colores = colores, sampleId = "paperName", distance = dist) 



plotHM.CvsP.FC14.dist.test<- makeHeatMap.dist(dataNormDE=dataNorm.CvsP.FC14.test, ed=final_ed2, group="group",  
                                              mainTit=paste0("log2(FC=1.4), Control vs sPE - ",  nrow(dataNorm.CvsP.FC14), " genes ", dist,  " dist."),
                                              colores = colores, sampleId = "paperName", distance = dist) 


## PCA

pca.CvsP.FC12.test<- pcaShape(dataNormDE=dataNorm.CvsP.FC12.test, ed=final_ed2, group="group",
                              mainTit=paste0("log2(FC=1.2), Control vs sPE - ", nrow(dataNorm.CvsP.FC12), " genes"),
                              shape="none", col=colores)

pca.CvsP.FC12.test$plotPCA

pca.CvsP.FC14.test<- pcaShape(dataNormDE=dataNorm.CvsP.FC14.test, ed=final_ed2, group="group",
                              mainTit=paste0("log2(FC=1.4), Control vs sPE - ", nrow(dataNorm.CvsP.FC14), " genes"),
                              shape="none", col=colores)
pca.CvsP.FC14.test$plotPCA




# 7. Analysis GO ----------------------------------------------------------

## GO tables

mydataFC.UP<-  subset(sPEctrl.treat.FC14.deg.EG, logFC>(log2(1.4))) 
mydataFC.DOWN<- subset(sPEctrl.treat.FC14.deg.EG, logFC<(log2(1.4))) 

## UP

goUP.temp<- goana(mydataFC.UP$EntrezGene,species="Hs", geneid = "ENTREZID")
goUP.top<- topGO(goUP.temp, sort="DE", n=Inf)

## Adjusted p-values

myPvalues.UP<- goUP.top$P.DE
names(myPvalues.UP)<- rownames(goUP.top)
myPvalues.UP.adj<- p.adjust(myPvalues.UP,method = "BH")
myPvalues.UP.adj<- data.frame(myPvalues.UP.adj)

goUP.top.sig<- cbind(goUP.top, myPvalues.UP.adj[rownames(goUP.top),])
colnames(goUP.top.sig)[6]<- "FDR"
goUP.top.sig<- goUP.top.sig[goUP.top.sig$FDR<0.05,] ## 0 GO


## DOWN

goDOWN.temp<- goana(mydataFC.DOWN$EntrezGene,species="Hs", geneid = "ENTREZID")
goDOWN.top<- topGO(goDOWN.temp, sort="DE", n=Inf)

myPvalues.DOWN<- goDOWN.top$P.DE
names(myPvalues.DOWN)<- rownames(goDOWN.top)
myPvalues.DOWN.adj<- p.adjust(myPvalues.DOWN,method = "BH")
myPvalues.DOWN.adj<- data.frame(myPvalues.DOWN.adj)

goDOWN.top.sig<- cbind(goDOWN.top, myPvalues.DOWN.adj[rownames(goDOWN.top),])
colnames(goDOWN.top.sig)[6]<- "FDR"
goDOWN.top.sig<- goDOWN.top.sig[goDOWN.top.sig$FDR<0.05,] ## 168 GO


# 7.1 Maps Between Entrez Gene ID and GO --------------------------------------

# org.Hs.egGO2ALLEGS --> Maps between Entrez Gene IDs and Gene Ontology (GO) IDs

EG.GO<- AnnotationDbi::toTable(get("org.Hs.egGO2ALLEGS"))
dim(EG.GO) # 3430396       4

d<- duplicated(EG.GO[,c("gene_id", "go_id", "Ontology")])

EG.GO <- EG.GO[!d, ]

# Now choosing DE genes for each GO ID/Ontology combination:
de.by.go <- split(EG.GO$gene_id,EG.GO$go_id)

# 7.2 Genes UP ------------------------------------------------------------

## There are not GO for UP genes


# 7.3 Genes DOWN ----------------------------------------------------------

## Genes DE DOWN for each GO

deDOWN<-mydataFC.DOWN$EntrezGene ## 
deDOWN<-deDOWN[!is.na(deDOWN)] 

de.by.goDOWN <- lapply(de.by.go, FUN=function(x) { x[x %in% deDOWN] })
de.by.goDOWN<- de.by.goDOWN[names(de.by.goDOWN)  %in% rownames(goDOWN.top.sig)] 

goDOWN.top1<- goDOWN.top.sig
goDOWN.top1$ENTREZID<- "kk"
goDOWN.top1$SYMBOL<-"kk"

for (i in 1:nrow(goDOWN.top1)) {
  # i=1 
  myGO<-rownames(goDOWN.top1)[i]
  myEId<- de.by.goDOWN[[myGO]] ## EntrezId
  mySId<- mydataFC.DOWN[mydataFC.DOWN$EntrezGene %in% myEId, "hgnc"]
  goDOWN.top1[i,"ENTREZID"]<- paste0(myEId, collapse = "|")
  goDOWN.top1[i,"SYMBOL"]<- paste0(mySId, collapse = "|")
  
}

## Tables in Excel

hs <- openxlsx::createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize=12,
                            fontName="Arial Narrow", fgFill = "#4F80BD")

openxlsx::write.xlsx(list("goDOWN"=goDOWN.top1),
                     file = paste0(pipeRoot, "genesSymbolUpDownFC14.xlsx"),
                     rowNames=T, headerStyle = hs)

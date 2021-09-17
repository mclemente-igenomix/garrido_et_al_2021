##############################################################
##############################################################
## Clinical Variables Analysis
## DE analysis: term vs preterm


# Libraries ---------------------------------------------------------------

## Functions and libraries

source(file=paste0(pipeRoot, "functionsPREECLAMPSIA.R"))



# Reading files -----------------------------------------------------------

## This script requires the following files:
#  ed.rev.cli: Experimental design you can load from github
#  final_gcts: Raw counts, obtained previously or you can load from GEO/github

ed.rev.cli<- read.csv(file=paste0(pipeRoot, "edRevPaper.csv"))
final_gcts<- read.csv(file=paste0(pipeRoot, "rawData.csv"), row.names = 1,check.names = F)
dim(final_gcts) ## 56638    40


## Remove sPE samples

ed.term<-ed.rev.cli[ed.rev.cli$groupB=="Control_NO",] ## 16 samples


# 1. Normalization ----------------------------------------------


final_gcts.term <- final_gcts[, colnames(final_gcts) %in% ed.term$sample]
dim(final_gcts.term) # 56638    16

group<- factor(ed.term$Preterm, levels = c("PRETERM","TERM"))

y.term<- DGEList(counts=final_gcts.term, group=group)

min_samples <- round(mean(table(ed.term$Preterm)))/2 

keep <- rowSums(cpm(y.term) > 1) >= min_samples 
table(keep)

y.term <- y.term[keep, keep.lib.sizes=FALSE]
y.term.norm <- calcNormFactors(y.term, method="TMM")

# 2. GLM models --------------------------------------------------------------

# 2.1 Design Matrix -------------------------------------------------------

design<- model.matrix(~0+group) ## Sin intercept
rownames(design)<- colnames(y.term.norm)

# 2.2 Estimating Dispersion -----------------------------------------------

## Estimating dispersion

y.term.norm<-estimateDisp(y.term.norm, design, robust=TRUE)


# 2.3 PCA plot term  ---------------------------------

## Data exploration: Multi-dimensional scaling (MDS) plot

points <- c(16,16)
colors <- rep(c("blue", "orange"), 2)

plotMDS(y.term.norm, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)


final_logcpm.term <- cpm(y.term.norm, prior.count=2, log=TRUE)
final_logcpm.term<- as.data.frame(final_logcpm.term)

final_fit_norm <- prcomp(t(na.omit(final_logcpm.term)))
final_pci_norm <- data.frame(final_fit_norm$x, 
                             color_group=ed.term[match(colnames(final_logcpm.term), ed.term$sample), "Preterm"],
                             paperName= ed.term[match(colnames(final_logcpm.term), ed.term$sample), "paperName"])

colores<- c("#2171b5","#6baed6") ## Blue light and blue dark

PCA12.norm<- ggplot(final_pci_norm, aes(x=PC1, y=PC2, color=color_group)) + 
  geom_point(size=3) +
  scale_x_continuous(name=paste("PC1 (", abs(round(summary(final_fit_norm)$importance[2, "PC1"]*100)), "%)", sep="")) +
  scale_y_continuous(name=paste("PC2 (", abs(round(summary(final_fit_norm)$importance[2, "PC2"]*100)), "%)", sep="")) +
  labs(title="PCA norm counts")+
  scale_color_manual(name = "Group", values = colores) +
  theme_light() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.title = element_text(hjust=0.5))+
  geom_text_repel(aes(x=PC1, y=PC2, color=color_group, label = paperName))


# 2.4 Fit model ---------------------------------------------------------------

fit<-glmQLFit(y.term.norm, design,robust=TRUE)
con<- makeContrasts( groupPRETERM-groupTERM, levels = design)

qlf<- glmQLFTest(fit, contrast=con)
topTags(qlf) 

summary(decideTests(qlf))
#         1*groupPRETERM -1*groupTERM
# Down                             0
# NotSig                       18476
# Up                               0

plotMD(qlf)


# 2.5 Treat model ---------------------------------------------------------

##################
### fc=1.2
##################

tr.FC12 <- glmTreat(fit, contrast=con, lfc=log2(1.2))

sPEctrl.treat.FC12<- topTags(tr.FC12, n=nrow(y.term.norm), p.value=1, adjust.method = "BH" )
sPEctrl.treat.FC12.deg<- sPEctrl.treat.FC12$table

summary(decideTests(tr.FC12))

dataPlot<- sPEctrl.treat.FC12$table %>%
  mutate(threshold=ifelse(FDR<0.05 & logFC>=log2(1.2), "UP",
                          ifelse(FDR<0.05 & logFC<log(1.2), "DOWN", "Non-Significant")))

## Volcano plot with all non-significant genes


VP.FC12<- ggplot(data=dataPlot, aes(x=logFC, y=-log10(FDR))) + 
  geom_point(aes(colour = threshold),alpha=0.8, size=2.5) + 
  geom_hline(yintercept = -log10(0.05), linetype="solid", colour="darkgrey") + 
  theme_bw() + theme(legend.position = "right") + 
  scale_colour_manual(values = c("UP"= "red", "DOWN"= "blue","Non-Significant"= "grey")) + 
  xlab("log2 fold change") + 
  ylab("-log10 FDR") +
  labs(title="TERM vs PRETERM: log2(FC=1.2) & FDR< 0.05 \n  ")


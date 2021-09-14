###################################################################
###################################################################
### Functions used at scripts



# Libraries ---------------------------------------------------------------

## DE analysis

suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(htmltools))
suppressMessages(library(openxlsx))
suppressMessages(library(DT))
suppressMessages(library(ggplot2))
suppressMessages(library(edgeR))
suppressMessages(library(viridis))
suppressMessages(library(xtable))
suppressMessages(library(ggrepel))
suppressMessages(library(pheatmap))
suppressMessages(library(Biobase))


# GO analysis

library(org.Hs.eg.db)
library(GO.db)



# Functions ---------------------------------------------------------------

get.df2plot <- function(df, ed) {
  df$gene <- rownames(df)
  df <- melt(df, id="gene")
  colnames(df) <- c("gene", "sample", "counts")
  df <- merge(df, ed, by = "sample")
  df$sample <- factor(df$sample, levels=unique(df$sample[order(df$group)]))
  return(df)
}


get.normalization <- function(ds_samples, ed, min_samples){  
  
  # INPUT
  # ds_samples: Raw counts, col:samples and row:genes. Mantener mismo orden en las col que ed$group 
  # ed: experimental design where group is at the column=group
  # min_samples: round(mean(table(ed$group)))/2
  
  # OUTPUT
  # y: Normalize Count
  
  y<- DGEList(ds_samples, group = ed$group)
  keep <- rowSums(cpm(y) > 1) >= min_samples ## CPM of 1 corresponds to a count of 6-7 and expressed at least min_samples
  y <- y[keep, keep.lib.sizes=FALSE]
  y <- calcNormFactors(y, method="TMM")
  return(y)
}


estimate.disp.pairw <- function(y){
  
  # INPUT
  # y: DGEList object, obtained previously. Norm counts.
  
  ## OUTPUT
  # y1: Add some components to previous object. Returns a list containing common.dispersion, trended.dispersion, 
  # tagwise.dispersion, span, prior.df and prior.n. 
  
  y1<-estimateDisp(y)
  cat("Dispersion", y1$common.dispersion, "\n")
  cat("Prior.n", y1$prior.n, "\n")
  return(y1)
}


pcaShape<- function(dataNormDE, ed, group, shape, mainTit, col){
    ## INPUT
    ## - dataNormDE: Normalize data applying cpm and only DE --> data.logcpm[mygenesDE,] 
    ## - ed: experimental design where group is at the column=group. Example --> nera_ed
    ## - group: "group". Si tienes mas clases
    ## - shape: "shape".Si tienes mas clases or "none"
    ##  - mainTit: Title
    ## - col: colores
    ## OUTPUT
    ## - Dataframe
    
    colores<- col
    #colores<- c("#1f78b4","#ff7f00") ## Azul Naranja-fuerte
    clasesCol<-levels(ed[,group])
    names(colores)<- clasesCol
    
    fit <- prcomp(t(na.omit(dataNormDE)), scale=F)
    
    if (shape=="none"){
      
      pci <- data.frame(fit$x, color_group=ed[rownames(fit$x),group])
      
      p_pca <- ggplot(pci) + geom_point(aes(x=PC1, y=PC2, color=color_group), size=3) 
      p_pca <- p_pca + scale_x_continuous(name=paste("PC1 (", abs(round(summary(fit)$importance[2, "PC1"]*100)), "%)", sep=""))
      p_pca <- p_pca + scale_y_continuous(name=paste("PC2 (", abs(round(summary(fit)$importance[2, "PC2"]*100)), "%)", sep=""))
      # p_pca <- p_pca + scale_color_viridis(discrete = TRUE)
      p_pca <- p_pca + scale_color_manual(name = "Group", values = colores) 
      p_pca <- p_pca + theme_light()
      p_pca <- p_pca + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      p_pca<- p_pca + geom_text_repel(aes(x=PC1, y=PC2, color=color_group, label = rownames(pci)))
      p_pca<- p_pca + labs(title=mainTit)
      
    } else {
      
      pci <- data.frame(fit$x, color_group=ed[rownames(fit$x),group],
                        shape_group=ed[rownames(fit$x),shape])
      p_pca <- ggplot(pci) + geom_point(aes(x=PC1, y=PC2, color=color_group, shape=shape_group), size=3) 
      p_pca <- p_pca + scale_x_continuous(name=paste("PC1 (", abs(round(summary(fit)$importance[2, "PC1"]*100)), "%)", sep=""))
      p_pca <- p_pca + scale_y_continuous(name=paste("PC2 (", abs(round(summary(fit)$importance[2, "PC2"]*100)), "%)", sep=""))
      # p_pca <- p_pca + scale_color_viridis(discrete = TRUE)
      p_pca <- p_pca + scale_color_manual(name = "Group", values = colores) 
      p_pca <- p_pca + theme_light()
      p_pca <- p_pca + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      p_pca<- p_pca + geom_text_repel(aes(x=PC1, y=PC2, color=color_group, label = rownames(pci)))
      p_pca<- p_pca + labs(title=mainTit)
      
    }
    
    
    
    return(list("dataFramePCA"=pci, "plotPCA"=p_pca))
    
  }  


makeHeatMap.3g<- function(dataNormDE, ed, group, mainTit, colores, sampleId) {
  
  ## INPUT
  ## - dataNormDE: Normalize data applying cpm and only DE --> data.logcpm[mygenesDE,] 
  ## - ed: experimental design where group is at the column=group. Example --> nera_ed
  ## - group: "group". Si tienes mas clases
  ## - sampleId: sampleId del ed
  ## - mainTit: main Title
  ## - colores: colors to use
  ## OUTPUT
  ## A list object type hclust
  
  col_group_heatmap <- ed[colnames(dataNormDE), group]
  sampl_heatmap <- ed[colnames(dataNormDE), sampleId]
  mat_colors<-  list(group=colores)
  names(mat_colors$group) <- levels(ed[, group])
  mat_col <- data.frame(group = col_group_heatmap)
  rownames(mat_col) <- sampl_heatmap
  
  pheatmap(
    mat               = dataNormDE,
    color             = viridis(200),
    border_color      = NA,
    show_colnames     = T,
    show_rownames     = F, ## Si hay pocos genes:T
    annotation_col    = mat_col,
    annotation_colors = mat_colors,
    drop_levels       = TRUE,
    fontsize          = 8,
    main              = mainTit)
  
}


makeHeatMap.dist<- function(dataNormDE, ed, group, mainTit, colores, sampleId, distance) {
  
  ## INPUT
  ## - dataNormDE: Normalize data applying cpm and only DE --> data.logcpm[mygenesDE,] 
  ## - ed: experimental design where group is at the column=group. Example --> nera_ed
  ## - group: "group". Si tienes mas clases
  ## - sampleId: sampleId del ed
  ## - mainTit: main Title
  ## - colores: colors to use
  ## - distance: It could be "euclidean","correlation","manhattan","canberra","minkowski",etc.
  ## OUTPUT
  ## A list object type hclust
  
  col_group_heatmap <- ed[colnames(dataNormDE), group]
  sampl_heatmap <- ed[colnames(dataNormDE), sampleId]
  mat_colors<-  list(group=colores)
  names(mat_colors$group) <- levels(ed[, group])
  mat_col <- data.frame(group = col_group_heatmap)
  rownames(mat_col) <- sampl_heatmap
  
  pheatmap(
    mat               = dataNormDE,
    color             = viridis(200),
    border_color      = NA,
    show_colnames     = T,
    clustering_distance_cols = distance,
    show_rownames     = F, ## Si hay pocos genes:T
    annotation_col    = mat_col,
    annotation_colors = mat_colors,
    drop_levels       = TRUE,
    fontsize          = 8,
    main              = mainTit)
  
}




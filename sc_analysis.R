library(BiocParallel)
library(data.table)
library(Seurat)
library(scTarNet)
library(sctransform)
library(igraph)
library(topGO)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AUCell)
library(reticulate)
library(doParallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)
library(progeny)
options(bitmapType='cairo')

counts <- #### path to count matrix
annotation <- fread('TabulaMuris/annotationations_facs.csv')
counts <- as.matrix(fread(counts))
rownames(counts) <- counts[,1]
counts <- counts[,-1]
mode(counts) <- 'numeric'

seurat_obj <- Seurat::CreateSeuratObject(counts,min.cells = 10)

seurat_obj <- FindVariableFeatures(object = seurat_obj,nfeatures = 4000)
seurat_obj <- NormalizeData(seurat_obj,normalization.method = "LogNormalize")
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(object = seurat_obj,ndims.print = 1:5,assay = 'RNA',nfeatures.print = 5,do.print=T,npcs = 50)
seurat_obj <- FindNeighbors(object = seurat_obj)
ElbowPlot(seurat_obj,ndims=50)
DimPlot(seurat_obj,reduction = 'pca')
seurat_obj <- FindClusters(object = seurat_obj, algorithm = 4,graph.name = 'RNA_snn',resolution = 0.3)
seurat_obj <- RunUMAP(object = seurat_obj,graph = 'RNA_nn',umap.method = 'umap-learn')

new_ident <- annotation[,3:4][which(annotation$cell %in% colnames(seurat_obj))]
new_ident <- rbind(new_ident,data.table(cell= colnames(seurat_obj)[(which(!(colnames(seurat_obj) %in% annotation$cell )))],cell_ontology_class = 'unknown'))
new_ident <- setNames(object = new_ident$cell_ontology_class,nm = new_ident$cell)#
new_ident <- as.factor(new_ident)
old_ident <- seurat_obj@active.ident
seurat_obj@active.ident <- new_ident

paga_path <- #path to paga output
paga_embeds <- fread(paga_path)
old.embeds <- seurat_obj@reductions$umap@cell.embeddings
paga_embeds <- as.matrix(paga_embeds)
rownames(paga_embeds) <- colnames(seurat_obj)
colnames(paga_embeds) <- c('UMAP_1','UMAP_2')

seurat_obj@reductions$umap@cell.embeddings <- paga_embeds

markers <- FindAllMarkers(object = seurat_obj,assay = 'RNA',test.use = 'MAST')
top20 <- markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 20, wt = avg_logFC)

png(file =paste0(output_directory,'/','heatmap_',dataset_name,'.png'),units = 'in',res = 300,type = 'cairo',width = 15,height =12)
DoHeatmap(seurat_obj,features = top20$gene, slot = 'data',draw.lines = T)+
  ggplot2::scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))
dev.off()



# sctar -------------------------------------------------------------------

TFs <- c("Erg","Fli1","Tal1","Lyl1","Lmo2","Runx1","Cbfb","Gata2","Gata1","Spi1","Ldb1","Cbfa2t3")
TFs %in% rownames(seurat_obj)
dCorEdges <- calculateTFstoTargets(as.matrix(seurat_obj@assays$RNA@data), TFs = TFs, n.cores=28, mt_correction="bon=0.05")
Interactions <- calculateConditionalCors(as.matrix(seurat_obj@assays$RNA@data), TFs, dCorEdges, n.cores=28, threshold.interaction=0.01, bidirectional=F, threshold.indirect=0.1, exclude.indirect=F)
Int <- combineInteractions(list(Interactions$Int))
Dep <- combineDependencies(list(Interactions$Dep))
graph_info <- plotInteractionsWithTargets(Dep, Int, recurr_threshold=1)
int_path <- #path to Int output
dep_path <- #path to Dep output
interactions <- fread(int_path)
dependencies <- fread(dedep_path)

zx_pos <- dependencies[Target %in% interactions$Target & 
                         (Gene %in% interactions$TF1 | Gene %in% interactions$TF2) &
                         direction ==1,-3,]
zx_neg <- dependencies[Target %in% interactions$Target & 
                         (Gene %in% interactions$TF1 | Gene %in% interactions$TF2) &
                         direction ==-1,-3,]

pos_corr <- split((zx_pos$Target),factor((zx_pos$Gene)))
neg_corr <- split((zx_neg$Target),factor((zx_neg$Gene)))
names(pos_corr) <- paste0(names(pos_corr),'_pos')
names(neg_corr) <- paste0(names(neg_corr),'_neg')


# heatmaps with target genes overlappings ------------------------------------------

heat_overlap <- function(vec1,vec2){
  targ_genes <- expand.grid(names(vec1),names(vec2),stringsAsFactors = F)
  targ_genes$overlap <- apply(targ_genes,MARGIN = 1,function(X){
    return(length(intersect(vec1[X[1]][[1]],vec2[X[2]][[1]])))
  })
  targ_genes <- targ_genes[order(targ_genes$overlap,decreasing = T),]
  colnames(targ_genes)[3] <- 'common genes'
  targ_genes[targ_genes$`common genes`>0,]
  t_mouse <- matrix(targ_genes[order(as.integer(rownames(targ_genes))),]$`common genes`,length(vec1),length(vec2))
  dimnames(t_mouse) <- list(names(table(targ_genes[order(as.integer(rownames(targ_genes))),]$Var1)),
                            names(table(targ_genes[order(as.integer(rownames(targ_genes))),]$Var2)))
  p <- pheatmap::pheatmap(t_mouse,color = rev(colorspace::sequential_hcl(n=200,"YlOrRd")),
                          border_color = 'black',clustering_method = "ward.D2",
                          display_numbers = T, angle_col = 315,
                          number_format = "%.f", number_color = "black",
                          treeheight_row = 0, treeheight_col = 0)
  return(p)
}

png(file =paste0(output_directory,'/','overlap++_',dataset_name,'.png'),units = 'in',res = 200,type = 'cairo',width = 6,height =4)
heat_overlap(pos_corr,pos_corr)
dev.off()
png(file =paste0(output_directory,'/','overlap+-_',dataset_name,'.png'),units = 'in',res = 200,type = 'cairo',width = 6,height =4)
heat_overlap(pos_corr,neg_corr)
dev.off()
png(file =paste0(output_directory,'/','overlap--_',dataset_name,'.png'),units = 'in',res = 200,type = 'cairo',width = 6,height =4)
heat_overlap(neg_corr,neg_corr)
dev.off()

# aucel -------------------------------------------------------------------



cells_rankings <- AUCell_buildRankings(as.matrix(seurat_obj@assays$RNA@data), nCores=14, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(pos_corr, cells_rankings)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=TRUE)
selectedThresholds <- getThresholdSelected(cells_assignment)
library(shiny); library(rbokeh)
aucellApp <- AUCell_createViewerApp(auc=cells_AUC,thresholds=selectedThresholds, exprMat=as.matrix(seurat_obj@assays$RNA@data))
options(shiny.host="0.0.0.0")
savedSelections <- runApp(aucellApp)

for (cell in names(cells_assignment)){
  savedSelections$thresholds[cell]
  cells_assignment[[cell]]$aucThr$selected[[1]] <- savedSelections$thresholds[cell][[1]]
  newSelectedCells <- names(which(getAUC(cells_AUC)[cell,]>savedSelections$thresholds[cell]))
  base <- rep(1,dim(cells_rankings)[1])
  base[which(rownames(cells_rankings) %in% pos_corr[[cell]])] <- 1:length(pos_corr[[cell]])+1
  model <- Rfast::kruskaltests(cells_rankings@assays@data$ranking,base)
  padj <- (p.adjust(model[,2],method = 'holm'))
  cells_assignment[[cell]]$assignment <- intersect(newSelectedCells,colnames(seurat_obj)[which(padj<0.05)])
  print(length(cells_assignment[[cell]]$assignment))
}



cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
library(NMF)
pdf(file =paste0(output_directory,'/auc_heat_',dataset_name,'.pdf'),width = 12,height = 9)
pheatmap::pheatmap(assignmentMat, scale="none",  legend=FALSE)
dev.off()

paga_embeds <- seurat_obj@reductions$umap@cell.embeddings
selectedThresholds <- getThresholdSelected(cells_assignment)
p <- list()
j <- 0
png(file =paste0(output_directory,'_positive','.png'),width = 20,height = 16,units = 'in',res = 300,type = 'cairo',bg = 'black')
for(geneSetName in names(selectedThresholds)){
  nBreaks <- 5 # Number of levels in the color palettes
  # Color palette for the cells that do not pass the threshold
  colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
  # Color palette for the cells that pass the threshold
  colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
  
  # Split cells according to their AUC value for the gene set
  passThreshold <- (getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]) & 
    names(getAUC(cells_AUC)[geneSetName,])%in%cellsAssigned[[geneSetName]]
  table(passThreshold)
  if(sum(passThreshold) > 0)
  {
    aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)
    
    # Assign cell color
    cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])), 
                   setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))
    
    # Plot
    cellsTsne <- paga_embeds[names(sort(getAUC(cells_AUC)[geneSetName,])),]
    cellsTsne <- data.frame(cellsTsne,
                            cell_color=rownames(cellsTsne),
                            thr=passThreshold[rownames(cellsTsne)],
                            cluster=seurat_obj@active.ident[rownames(cellsTsne)])
    clus_labels <- data.frame(x=numeric(),y=numeric(),name=character())
    for(i in (levels(cellsTsne$cluster))){
      centr_x <- setNames(mean(cellsTsne[cellsTsne$cluster==i,]$UMAP_1),i)
      centr_y <- setNames(mean(cellsTsne[cellsTsne$cluster==i,]$UMAP_2),i)
      clus_labels[i,] <- data.frame(x=centr_x,y=centr_y)
      
    }
    clus_labels$name <- rownames(clus_labels)
    g <- ggplot(cellsTsne,aes(x=UMAP_1,y=UMAP_2))+
      geom_point(aes(color=cell_color,shape=thr),show.legend = T)+
      scale_shape_discrete(name="Cells passed threshold",labels=c("No","Yes"))+
      scale_color_manual(values =(cellColor) ,guide="none")+guides(shape="legend")+
      ggrepel::geom_label_repel(data=clus_labels,aes(x=x,y=y,label=name),
                                segment.color = 'white',
                                segment.colour = 'black',
                                segment.size = 1,
                                box.padding = unit(0.95, "lines"),size=3,
                                alpha=0.85)+
      ggtitle(geneSetName)+
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),axis.text = element_blank(),
            panel.background = element_rect(fill='gray'),
            panel.grid = element_line(linetype = 'dashed',size = 0.1),
            text = element_text(colour = 'white'),
            plot.background = element_rect(fill = 'black'),
            title = element_text(color = '#ffd900',face = 'bold'),
            legend.title = element_text(color = "black"),legend.text = element_text(color='black'),
            legend.background = element_rect(fill = "gray"))
    j <- j+1
    p[[j]] <- g
    #plot(cellsTsne_new, main=geneSetName,sub=selectedThresholds[geneSetName],col=cellColor[rownames(cellsTsne_new)], pch=16,cex=1) 
  }
}
do.call(ggpubr::ggarrange,c(p,common.legend = T))
dev.off()

ggplot(NULL,aes(x=tissue,y=value))+geom_bar(data=coexpression_fin,aes(fill=factor(ntf)),
                                            stat="identity",position = position_fill(reverse = T))+
  scale_fill_viridis_d()+
  geom_errorbar(data = data.frame(value=wmeans_fin,tissue=unique(coexpression_fin$tissue)),aes(x=tissue,ymin=value,ymax=value),color="red")+
  labs(fill="number of coexpressed TFs",x=element_blank(),y="relative coexpression")+
  theme(axis.text.x = element_text(angle = 315,hjust = -0.021),
        axis.text.y = element_text(size = 12),axis.ticks.y = element_blank(),
        panel.grid = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_blank(),
        axis.line.y.left  = element_blank())+coord_flip()+DarkTheme()


wmeans_fin <- unlist(lapply(levels(coexpression_fin$tissue),function(X){
  (coexpression_fin[coexpression_fin$tissue==X,"value"]%*%(as.numeric(coexpression_fin[coexpression_fin$tissue==X,"ntf"])))
}))
wmeans_fin <- wmeans_fin/max(wmeans_fin)
wmeans <- wmeans/max(wmeans)
#Gene Ontology analysis -----------------------------------------------------------
go_script <- function(targets,genes){
  mm <- org.Mm.eg.db
  symbols <- targets
  symbols_all <- genes
  entrezs <- select(mm,
                    keys = symbols,
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL")
  entrezs_all <- select(mm,
                        keys = symbols_all,
                        columns = c("ENTREZID", "SYMBOL"),
                        keytype = "SYMBOL")
  geneNames <- entrezs$ENTREZID
  geneUniverse <- entrezs_all$ENTREZID
  geneList <- factor(as.integer(geneUniverse %in% geneNames))
  names(geneList) <- geneNames
  GOdata <- new("topGOdata",
                description = "GO",
                ontology = "BP",
                allGenes = geneList,
                annotation = annFUN.org,
                nodeSize = 7,
                mapping="org.Mm.eg.db")
  elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks",cutOff=0.05)
  allResKS <- GenTable(GOdata, `p-value` = elim.ks, orderBy = "p-value",  topNodes = 50)
  #printGraph(object = GOdata,firstSigNodes = 5,result = elim.ks,pdfSW=T)
  return(list(allResKS,GOdata))
}


# integration -------------------------------------------------------------

seurat.list <- lapply(X = list.files('filtered/',full.names = T,pattern = "-norm.csv$")[-c(2,12)], FUN = function(x) {
  tissue <- gsub(".*\\/(.*)\\-.*", "\\1",x)
  X <- as.matrix(fread(x))
  rownames(X) <- X[,1]
  X <- X[,-1]
  mode(X) <- 'numeric'
  new_ident <- annotation[,3:4][which(annotation$cell %in% colnames(X))]
  new_ident <- rbind(new_ident,data.table(cell= colnames(X)[(which(!(colnames(X) %in% annotation$cell )))],cell_ontology_class = 'unknown'))
  X <- X[,new_ident[grep("endothelial",new_ident$cell_ontology_class),cell]]
  new_ident <- setNames(object = new_ident$cell_ontology_class,nm = new_ident$cell)#
  new_ident <- as.factor(new_ident)
  X <- CreateSeuratObject(X)
  print(tissue)
  X <- FindVariableFeatures(X, nfeatures=6000,verbose = FALSE)
  #X <- ScaleData(X,verbose = FALSE)
  X@active.ident <- as.factor(setNames(object = rep(tissue,ncol(X)),nm = colnames(X)))
  return(X)
})

seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  return(x)
})

anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:40,anchor.features = 6000,
                                  verbose = F,k.filter = 40)

seurat_integrated <- IntegrateData(anchorset = anchors)
seurat_integrated <- NormalizeData(seurat_integrated,assay = 'integrated', verbose = FALSE)
seurat_integrated <- ScaleData(seurat_integrated,assay = 'integrated', verbose = FALSE)
seurat_integrated <- RunPCA(object = seurat_integrated, verbose = F)
PCAPlot(seurat_integrated,group.by="ident")
ElbowPlot(seurat_integrated,ndims = 50)
seurat_integrated$orig.ident <- seurat_integrated@active.ident
seurat_integrated <- FindNeighbors(object = seurat_integrated,"pca",1:45)
seurat_integrated <- FindClusters(seurat_integrated,graph.name = "integrated_snn",algorithm = 4,resolution = 0.65)
seurat_integrated <- RunTSNE(seurat_integrated,dims = 1:45)
seurat_integrated <- RunUMAP(seurat_integrated,dims = 1:45)

integ_markers <- FindAllMarkers(seurat_integrated)
top10_integ <- integ_markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_logFC)
png(file =paste0(output_directory,'heatmap_integrated.png'),units = 'in',res = 300,type = 'cairo',width = 30,height =23)
DoHeatmap(seurat_integrated,features = top10_integ$gene,slot = 'scale.data',draw.lines = F,
          group.colors = unname(Polychrome::light.colors()),raster = F,group.by = "disease")+
  ggplot2::scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))+
  DarkTheme()+theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())#+NoAxes()
dev.off()

seurat_integrated@active.ident <- ident

Runx1cells <- seurat_obj@active.ident[intersect(grep(pattern = "endo",seurat_obj@active.ident),which(seurat_obj@assays$RNA@data["Runx1",]>0))]

seurat_obj@active.ident[which(!names(seurat_obj@active.ident)%in%names(Runx1cells))]
xident <- c(setNames(as.character(seurat_obj@active.ident[which(!names(seurat_obj@active.ident)%in%names(Runx1cells))]),names(seurat_obj@active.ident[which(!names(seurat_obj@active.ident)%in%names(Runx1cells))])),
            (setNames(rep("x cells",340),
                      names(seurat_obj@active.ident[which(names(seurat_obj@active.ident)%in%names(Runx1cells))]))))
xident <- as.factor(xident)
seurat_obj@active.ident <- xident
xmarkers <- FindMarkers(seurat_obj,ident.1 = "endothelial cells",ident.2 = "x cells",test.use = "poisson")

prop <- pheatmap::pheatmap(prop.table(table(seurat_integrated$orig.ident,seurat_integrated@active.ident),margin = 2)*100,cluster_cols = F,color = rev(colorspace::sequential_hcl(n=200,"YlOrRd")),
                           border_color = 'black',clustering_method = "ward.D2",
                           display_numbers = T, angle_col = 315,
                           number_format = "%.f", number_color = "black",
                           treeheight_row = 0, treeheight_col = 0)



# Progeny -----------------------------------------------------------------

seurat_integrated@active.ident <- seurat_integrated$orig.ident
CellsClusters <- data.frame(Cell = names(Idents(seurat_integrated)),
                            CellType = as.character(Idents(seurat_integrated)),
                            stringsAsFactors = FALSE)



seurat_integrated <-  progeny(seurat_integrated, scale=FALSE, organism="Mouse", perm=20,
                              return_assay = TRUE)
seurat_integrated <- Seurat::ScaleData(seurat_integrated, assay = "progeny")

## We transform Progeny scores into a data frame to better handling the results
progeny_scores_df <-
  as.data.frame(t(GetAssayData(seurat_integrated, slot = "scale.data",
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell)

## We match Progeny scores with the cell clusters.
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>%
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0,
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength,
                      max(summarized_progeny_scores_df),
                      length.out=floor(paletteLength/2)))
progeny_hmap <- pheatmap::pheatmap((summarized_progeny_scores_df[,-1]),cluster_rows = T,cluster_cols = T,
                                   color=myColor, breaks = progenyBreaks,
                                   angle_col = 45,clustering_method = "ward.D2",
                                   treeheight_col = 0,  border_color = NA)

# mouse embryo GSE139389 --------------------------------------------------

barcodes <- fread('barcodes.txt')
matrices <- lapply(list.files(pattern = "GSM"),fread)
count_matrix <- Reduce(cbind,matrices)
count_matrix <- as.matrix(fread('data/mouse_embryo_gse139389/matrix.csv'))
rownames(count_matrix) <- count_matrix[,1]
count_matrix <- count_matrix[,-1]
mode(count_matrix) <- 'numeric'

seurat_obj <- CreateSeuratObject(count_matrix,min.cells = 10)
seurat_obj <- SCTransform(seurat_obj,variable.features.n = 4000)
seurat_obj <- RunPCA(object = seurat_obj,npcs = 50)
seurat_obj <- FindNeighbors(object = seurat_obj,assay = 'SCT')
seurat_obj <- FindClusters(object = seurat_obj,algorithm = 4,graph.name = 'SCT_snn',resolution = 0.3)
seurat_obj <- RunUMAP(seurat_obj,graph = 'SCT_nn',umap.method = 'umap-learn',n.neighbors = 15)
DimPlot(seurat_obj,reduction = 'umap',cols = unname(Polychrome::light.colors()))+DarkTheme()+NoAxes()#+ggtitle("EMBRYO")

paga_embeds <- fread('PAGA/embeddings/gse139389.csv')
old.embeds <- seurat_obj@reductions$umap@cell.embeddings
paga_embeds <- as.matrix(paga_embeds)
rownames(paga_embeds) <- colnames(seurat_obj)
colnames(paga_embeds) <- c('UMAP_1','UMAP_2')
seurat_obj@reductions$umap@cell.embeddings <- paga_embeds
DoHeatmap(seurat_obj,features =c(top20$gene), slot = 'data',draw.lines = F,
          group.colors = unname(Polychrome::light.colors()),raster = F)+
  ggplot2::scale_fill_gradientn(colors = rev(viridis::inferno(10)))+
  DarkTheme()+theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())

markers <- FindAllMarkers(seurat_obj,test.use = "MAST")
top20 <- markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 20, wt = avg_logFC)

#Progeny -------------------------------------------------------------------------


endo <- as.matrix(seurat_obj@assays$RNA@data)[which(rownames(seurat_obj) %in% TFs[-c(9,10)]),
                                              grep('x cells',seurat_obj@active.ident)]


mode(endo) <- "logical"
coexpr3[[13]] <- table(colSums(endo==T))/dim(endo)[2]


for(kk in 1:13){
  for (k in c(1:4)) {
    p[[kk]]$gtable$grobs[[k]]$gp <- gpar(col="white",fontsize=10)
    p[[kk]]$gtable$grobs[[1]]$vjust <- (1)
  }
}
png(file =paste0('network_tables.png'),units = 'in',res = 300,type = 'cairo',width = 13,height =12)
grid.draw(rectGrob(gp=gpar(fill="black", lwd=0)))
do.call(gridExtra::grid.arrange,c(list(p[[1]][[4]]),
                                  list(p[[2]][[4]]),
                                  list(p[[3]][[4]]),
                                  list(p[[4]][[4]]),
                                  list(p[[5]][[4]]),
                                  list(p[[6]][[4]]),
                                  list(p[[7]][[4]]),
                                  list(p[[8]][[4]]),
                                  list(p[[9]][[4]]),
                                  list(p[[10]][[4]]),
                                  list(p[[11]][[4]]),
                                  list(p[[12]][[4]]),
                                  list(p[[13]][[4]]),newpage=F))
dev.off()
gridExtra::grid.arrange()

# gse137116 ---------------------------------------------------------------

counts <- Matrix::readMM('data/gse137116/GSE137116_gene_by_cell_count_matrix.txt')
cellann <- fread('data/gse137116/GSE137116_cell_annotationation.csv')
length(cellann$V1)
head(cellann)
geneann <- fread('data/gse137116/GSE137116_gene_annotationation.csv')
head(geneann)
rownames(counts) <- geneann$symbol
colnames(counts) <- cellann$V1
seurat_obj <- CreateSeuratObject(counts,min.cells = 10,min.features = 1500,meta.data = cellann)
rownames(cellann) <- cellann$V1
seurat_obj <- AddMetaData(object = seurat_obj,metadata = cellann)
rm(counts)
qplot(seurat_obj$nFeature_SCT,binwidth=40)
seurat_subset <- SCTransform(seurat_subset,variable.features.n = 6000)
seurat_obj <- RunPCA(object = seurat_obj,npcs = 50,rev.pca = T)
PCAPlot(object = seurat_obj,group.by='orig.ident')
ElbowPlot(seurat_obj,ndims=50)
seurat_obj <- FindNeighbors(object = seurat_obj,assay = 'SCT',dims = 1:45,reduction = "pca")
seurat_obj <- FindClusters(object = seurat_obj,algorithm = 4,graph.name = 'SCT_snn',resolution = 0.3)
seurat_obj <- RunUMAP(seurat_obj,graph = 'SCT_nn',umap.method = 'umap-learn')
seurat_obj@reductions$umap@cell.embeddings
DimPlot(seurat_obj,reduction = 'umap',group.by = "ident",cols = unname(Polychrome::light.colors()))+DarkTheme()+NoAxes()#+ggtitle("EMBRYO")
table(seurat_obj$Combined_Dataset)
seurat_subset <- (Seurat::SubsetData(seurat_obj,cells=unlist(Seurat::CellsByIdentities(seurat_obj,"E10.5 E+HE+IAC"),use.names = F)))
seurat_obj@active.ident <- as.factor(seurat_obj$Combined_Dataset)
rm(seurat_obj)
table(seurat_subset$Cell_type_refined)
grep("Endo",seurat_subset$Cell_type_refined)

# gse143637 ---------------------------------------------------------------

counts <- as.matrix(fread("data/gse143637/norm_matrix.csv"))

rownames(counts) <- counts[,1]
counts <- counts[,-1]
mode(counts) <- 'numeric'
counts <- counts[,colnames(seurat_subset)]

seurat_obj <- CreateSeuratObject(counts,min.cells = 10,min.features = 100)
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj,normalization.method = "LogNormalize")
seurat_obj <- SCTransform(seurat_obj,variable.features.n = 6000)
seurat_obj <- FindVariableFeatures(seurat_obj,nfeatures = "4000",selection.method = "vst")
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(object = seurat_obj,npcs = 50,rev.pca = T)
PCAPlot(object = seurat_obj)
ElbowPlot(seurat_obj,ndims = 40)
seurat_obj <- FindNeighbors(object = seurat_obj,assay = 'RNA',reduction = "pca",dims = 1:31)
seurat_obj <- FindClusters(object = seurat_obj,algorithm = 4,graph.name = 'RNA_snn',resolution = 0.3)
seurat_obj <- RunUMAP(seurat_obj,graph = 'RNA_nn',umap.method = 'umap-learn')
DimPlot(seurat_obj,reduction = 'umap',group.by = "ident",cols = unname(Polychrome::light.colors()))+DarkTheme()+NoAxes()#+ggtitle("EMBRYO")
markers <- FindAllMarkers(seurat_obj,test.use = "MAST")
top20 <- markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_logFC)
FeaturePlot(seurat_obj,c("Cdh5"),reduction = 'umap',order = T,slot = 'data',cols = c("green", "blue", "red"))
levels(seurat_obj@active.ident) <- c("1","Embryo_endothelial_143637","3","4","5","6")
seurat_subset <- (Seurat::SubsetData(seurat_obj,cells=unlist(Seurat::CellsByIdentities(seurat_obj,"Embryo_endothelial_143637"),use.names = F)))
png(file =paste0('data/gse143637/umap_clustered','.png'),units = 'in',res = 300,type = 'cairo',width = 7,height =5)
DimPlot(seurat_obj,reduction = 'umap',group.by = "ident" ,cols = unname(Polychrome::light.colors()))+DarkTheme()+NoAxes()
dev.off()
png(file =paste0('data/gse143637/heatmap','.png'),units = 'in',res = 300,type = 'cairo',width = 18,height =23)
DoHeatmap(seurat_obj,features = "Cdh5", slot = 'data',draw.lines = F,
          group.colors = unname(Polychrome::light.colors()),raster = F,group.by = "ident")+
  ggplot2::scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))+
  DarkTheme()+theme(axis.line.x.bottom = element_blank(),axis.line.y.left = element_blank())#+NoAxes()
dev.off()



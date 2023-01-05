library(pheatmap)
library(Seurat)
library(data.table)
library(ggplot2)
library(Matrix)
library(dplyr)
library(data.table)
library(gplots)
library(pheatmap)
library(tidyr)
library(data.table)
library(RColorBrewer)
setwd("/data/Delta_data4/kysbbubbu/melanomaWM237")
# Read data
merge.data <- Read10X(data.dir = "./AGG_melanoma/outs/filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(merge.data, project = "merge", names.field = 2, names.delim = "-", min.cells = 5, min.features = 200)
rm(merge.data)
current.cluster.ids <- c(1,2,3,4)
new.cluster.ids <- c("WM239A","113/6-4L", "131/4-5B1", "131/4-5B2")
pbmc@meta.data$orig.ident  <- plyr::mapvalues(x = pbmc@meta.data$orig.ident , from = current.cluster.ids, to = new.cluster.ids)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


##QC criteria
ggplot(pbmc@meta.data, aes(x=nFeature_RNA, fill=orig.ident, color=orig.ident)) +geom_density(alpha=0.3)+  theme_classic()+scale_color_manual(values=brewer.pal(4,"Set1"))+scale_fill_manual(values=brewer.pal(4,"Set1"))+theme(text = element_text(size=12))+theme(axis.title.x=element_text(color="black", size=16, face="bold"),axis.title.y=element_text(color="black", size=16, face="bold"),axis.text.x=element_text(color="black", size=14, face="bold"),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 12, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14))+ylab("Denstiy")+xlab("Number of genes per cell")+scale_color_manual(values=brewer.pal(4,"Set1"))+scale_fill_manual(values=brewer.pal(4,"Set1"))+ geom_vline(xintercept=c(3000, 7000),linetype=3)+NoLegend()

ggplot(pbmc@meta.data, aes(x=percent.mt, fill=orig.ident, color=orig.ident)) +geom_density(alpha=0.3)+  theme_classic()+scale_color_manual(values=brewer.pal(4,"Set1"))+scale_fill_manual(values=brewer.pal(4,"Set1"))+theme(text = element_text(size=12))+theme(axis.title.x=element_text(color="black", size=16, face="bold"),axis.title.y=element_text(color="black", size=16, face="bold"),axis.text.x=element_text(color="black", size=14, face="bold"),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=14, face="bold"), strip.text = element_text(size = 12, face="bold"),legend.title=element_blank(),legend.text=element_text(size=14))+ylab("Denstiy")+ylab("Denstiy")+xlab("Mitochondrial gene (%)")+scale_color_manual(values=brewer.pal(4,"Set1"))+scale_fill_manual(values=brewer.pal(4,"Set1")) + geom_vline(xintercept=c(20),linetype=3)+NoLegend()

pbmc <- subset(x = pbmc, subset = nFeature_RNA > 3000 & nFeature_RNA < 7000 & percent.mt < 20)
pbmc=NormalizeData(object = pbmc, normalization.method = "LogNormalize",  scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc=ScaleData(object = pbmc, features = rownames(pbmc), vars.to.regress = c("nCount_RNA", "percent.mt"))
pbmc <- RunPCA(object = pbmc)
pbmc <- ProjectDim(object = pbmc)
ElbowPlot(object = pbmc)
pbmc <- FindNeighbors(object = pbmc, dims = 1:12)
pbmc <- RunUMAP(pbmc, reduction.use = "pca", dims = 1:12)
pbmc <- FindClusters(object = pbmc, resolution = 0.1)

#Rename clusters 
current.cluster.ids <- c(0,1,2)
new.cluster.ids <- c("Brain","Primary","Lung")
pbmc@meta.data$RNA_snn_res.0.1  <- plyr::mapvalues(x = pbmc@meta.data$RNA_snn_res.0.1 , from = current.cluster.ids, to = new.cluster.ids)
pbmc@active.ident  <- plyr::mapvalues(x = pbmc@active.ident  , from = current.cluster.ids, to = new.cluster.ids)
pbmc@active.ident=factor(pbmc@active.ident, levels=c("Primary","Lung","Brain"))


#UMAP
DimPlot(pbmc,cols =brewer.pal(3,"Set2"))+theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ NoLegend()

DimPlot(pbmc,cols =brewer.pal(4,"Set1"),group.by="orig.ident")+theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ NoLegend()


#PCA
DimPlot(pbmc,reduction="pca",cols =brewer.pal(3,"Set2"))+theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+xlab("PC1 (14.2% of variance explained)")+ylab("PC2 (13.4% of variance explained)")

pca_var = as.matrix(pbmc[["pca"]]@stdev)
eigValues = (pca_var)^2  ## EigenValues
varExplained = eigValues / sum(eigValues)
sum(varExplained)
FeaturePlot(object = pbmc, reduction="pca", features = c("Pseudotime"))
FeaturePlot(object = pbmc, features = c("Pseudotime"))
FeaturePlot(object = pbmc, features = c("PC_2"))

table(pbmc@active.ident,pbmc@meta.data$orig.ident)
pbmc@meta.data$seurat_clusterspbmc@active.ident

current.cluster.ids <- c("Primary","Lung","Brain1","Brain2")
new.cluster.ids <- c("WM239A","113/6-4L", "131/4-5B1", "131/4-5B2")
pbmc@meta.data$orig.ident  <- plyr::mapvalues(x = pbmc@meta.data$orig.ident , from = current.cluster.ids, to = new.cluster.ids)




Pseudotime

###Technical check 
###Technical check 
###Technical check 
VlnPlot(pbmc, features = c("nFeature_RNA"), pt.size=0)+scale_color_manual(values=brewer.pal(4,"Set1"))+scale_fill_manual(values=brewer.pal(3,"Set2"))+theme(axis.title.x=element_blank(),axis.title.y=element_text(color="black", size=16, face="bold"),axis.text.x=element_text(color="black", size=14, face="bold",angle=0,hjust=0.5),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=14), strip.text = element_text(size = 12, face="bold"),plot.title=element_blank())+NoLegend()+ylab("Number of expressed genes")
VlnPlot(pbmc, features = c("nCount_RNA"), pt.size=0)+scale_color_manual(values=brewer.pal(4,"Set1"))+scale_fill_manual(values=brewer.pal(3,"Set2"))+theme(axis.title.x=element_blank(),axis.title.y=element_text(color="black", size=16, face="bold"),axis.text.x=element_text(color="black", size=14, face="bold",angle=0,hjust=0.5),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=14), strip.text = element_text(size = 12, face="bold"),plot.title=element_blank())+NoLegend()+ylab("Number of UMIs")
VlnPlot(pbmc, features = c("percent.mito"), pt.size=0)+scale_color_manual(values=brewer.pal(4,"Set1"))+scale_fill_manual(values=brewer.pal(3,"Set2"))+theme(axis.title.x=element_blank(),axis.title.y=element_text(color="black", size=16, face="bold"),axis.text.x=element_text(color="black", size=14, face="bold",angle=0,hjust=0.5),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=14), strip.text = element_text(size = 12, face="bold"),plot.title=element_blank())+NoLegend()+ylab("Mitochondrial gene (%)")



FeaturePlot(object = pbmc, features = c("nFeature_RNA"))+ ggtitle("Number of expressed genes") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 20))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title = element_text(size=14,face="bold"),axis.text = element_blank(),axis.ticks=element_blank(),legend.title=element_text(size=14),legend.text=element_text(size=12))
FeaturePlot(object = pbmc, features = c("nCount_RNA"))+ ggtitle("Number of UMIs") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 20))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title = element_text(size=14,face="bold"),axis.text = element_blank(),axis.ticks=element_blank(),legend.title=element_text(size=14),legend.text=element_text(size=12))
FeaturePlot(object = pbmc, features = c("percent.mito"))+ ggtitle("Mitochondrial gene (%)") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 20))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title = element_text(size=14,face="bold"),axis.text = element_blank(),axis.ticks=element_blank(),legend.title=element_text(size=14),legend.text=element_text(size=12))



# convert raw data to  csv 
setwd("/home/kysbbubbu/melanomaWM237")
library(data.table)
data_to_write_out0 <- as.data.frame(as.matrix(pbmc$RNA@counts))
fwrite(x = data_to_write_out0, row.names = TRUE, file = "raw_counts.csv")

# convert relative expression data to  csv 
library(data.table)
data_to_write_out <- as.data.frame(as.matrix(pbmc@assays$RNA@data))
data_to_write_out = round(data_to_write_out, digits= 3)
write.csv(x = data_to_write_out, row.names = TRUE, file = "data.csv")

# convert sacled expression data to  csv 
data_to_write_out0 <- as.data.frame(as.matrix(pbmc$RNA@scale.data))
data_to_write_out0 = round(data_to_write_out0, digits=3)
fwrite(x = data_to_write_out0, row.names = TRUE, file = "pbmcseuratscaledata.csv")

#PCA
pca_matrix=as.data.frame(pbmc@reductions$pca@feature.loadings)
write.csv(x = pca_matrix, row.names = TRUE, file = "pca.csv")
pca_matrix=pca_matrix[order(pca_matrix$PC_2),]
head(pca_matrix)

# convert relative expression data to  csv 
data_to_write_out <- as.data.frame(pbmc@meta.data)
write.csv(x = data_to_write_out, row.names = TRUE, file = "metadata.csv")



pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


#Theme
fplot=theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 26),axis.title = element_text(size=16,face="bold"),axis.text = element_blank(),axis.ticks=element_blank(),legend.text=element_text(size=12))

#DEG expression check (Lung/Brain up-regulated vs Primary)
FeaturePlot(object = pbmc, features = c("IGFBP7"))+ fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))
FeaturePlot(object = pbmc, features = c("CD74"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))
FeaturePlot(object = pbmc, features = c("MT2A"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))
FeaturePlot(object = pbmc, features = c("SLC16A3"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))
FeaturePlot(object = pbmc, features = c("FN1"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))

#DEG expression check (Primary up-regulated vs Lung/Brain)
FeaturePlot(object = pbmc, features = c("CCND1"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))
FeaturePlot(object = pbmc, features = c("TMSB4X"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))
FeaturePlot(object = pbmc, features = c("SFRP1"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))
FeaturePlot(object = pbmc, features = c("APOE"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))
FeaturePlot(object = pbmc, features = c("PMEL"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))

#DEG expression check (Lung up-regulated vs Brain)
FeaturePlot(object = pbmc, features = c("COL9A3"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))
FeaturePlot(object = pbmc, features = c("PRDX1"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))
FeaturePlot(object = pbmc, features = c("ID3"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))
FeaturePlot(object = pbmc, features = c("VGF"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))

#DEG expression check (Brain up-regulated vs Lung)
FeaturePlot(object = pbmc, features = c("NUPR1"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))
FeaturePlot(object = pbmc, features = c("VCX3B"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))
FeaturePlot(object = pbmc, features = c("DCT"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))
FeaturePlot(object = pbmc, features = c("SDCBP"))+fplot+ scale_colour_gradientn(colours=rev(brewer.pal(11,"RdBu")[1:11]))

##DEG test
cluster.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.5,logfc.threshold = 0.5)
meta.markers <- FindMarkers(object = pbmc, ident.1=c("Lung","Brain"),ident.2=c("Primary"),only.pos =F, min.pct = 0.5,logfc.threshold = 0.5)
brain.markers <- FindMarkers(object = pbmc, ident.1=c("Brain"),ident.2=c("Lung"),only.pos =F, min.pct = 0.5,logfc.threshold = 0.5)
fwrite(as.data.frame(cluster.markers),sep="\t","clusterDEG.txt",row.names=TRUE)
fwrite(as.data.frame(meta.markers),sep="\t","metaDEG.txt",row.names=TRUE)
fwrite(as.data.frame(brain.markers),sep="\t","brainDEG.txt",row.names=TRUE)
##Meta/primary DEGS
#300/450 (n=10)
DotPlot(pbmc, features=rev(c("IGFBP7","CD74","CCND1","TMSB4X","SFRP1","MT2A","SLC16A3","FN1","APOE","PMEL")),dot.min=0)+theme(axis.title=element_blank(), axis.text.x.bottom = element_text(size = 14,face="bold"),axis.ticks.x = element_blank(), axis.text.y.left=element_text(size = 14,face="bold"))+coord_flip()+ scale_size(breaks = c(0, 20, 40, 60,80),range = c(0,6))+scale_colour_gradientn(colors=colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100),values = rescale(c(-2, -1, 0, 1, 2)),limits=c(-2,2))+ NoLegend()


##Brain/Lung DEGs
#300/360 (N=8)
DotPlot(pbmc, features=rev(c("COL9A3","PRDX1","ID3","VGF","NUPR1","VCX3B","DCT","SDCBP")),dot.min=0)+theme(axis.title=element_blank(), axis.text.x.bottom = element_text(size = 14,face="bold"),axis.ticks.x = element_blank(), axis.text.y.left=element_text(size = 14,face="bold"))+coord_flip()+ scale_size(breaks = c(0, 20, 40, 60,80),range = c(0,6))+scale_colour_gradientn(colors=colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100),values = rescale(c(-2, -1, 0, 1, 2)),limits=c(-2,2))+ NoLegend()



###Heatmap showing expression of DEG markers for three clusters 
###Heatmap showing expression of DEG markers for three clusters 
###Heatmap showing expression of DEG markers for three clusters 
pbmc4=subset(x = pbmc, downsample = 1000)
pbmc4 <- ScaleData(pbmc4,features=rownames(pbmc4))
pbmc4@meta.data$active.ident=pbmc4@active.ident
pbmc4@meta.data=pbmc4@meta.data[sample(1:nrow(pbmc4@meta.data)),]
pbmc4@meta.data=pbmc4@meta.data[order(pbmc4@meta.data$active.ident),]
my_sample_col=as.data.frame(pbmc4@meta.data$active.ident)
rownames(my_sample_col)=rownames(pbmc4@meta.data)
colnames(my_sample_col)=c("Cluster")
my_sample_col$Cluster=factor(my_sample_col$Cluster, levels=c("Primary","Lung","Brain"))
write_out= as.data.frame(as.matrix(pbmc4@assays$RNA@scale.data[,rownames(pbmc4@meta.data)]))
write_out=write_out[cluster.markers$gene,]
my_sample_row=as.data.frame(rownames(write_out))
colnames(my_sample_row)=c("DEGs")
rownames(my_sample_row)=my_sample_row$DEGs
my_sample_row$DEGs=cluster.markers$cluster
Var2 =brewer.pal(3,"Set2")
names(Var2) = c("Primary","Lung","Brain")
ann_colors = list(DEGs=Var2,Cluster=Var2)
pheatmap(write_out,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-2, 2, by = 0.04),annotation_row =my_sample_row,annotation_col =my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100), legend = F, annotation_legend = F)
rm(pbmc4,write_out)
table(cluster.markers$cluster)

###Heatmap showing correlation coefficients for three clusters 
###Heatmap showing correlation coefficients for three clusters 
###Heatmap showing correlation coefficients for three clusters 
pbmc4=subset(x = pbmc, downsample = 1000)
pbmc4 <- ScaleData(pbmc4)
?ScaleData
pbmc4@meta.data$active.ident=pbmc4@active.ident
pbmc4@meta.data=pbmc4@meta.data[sample(1:nrow(pbmc4@meta.data)),]
pbmc4@meta.data=pbmc4@meta.data[order(pbmc4@meta.data$active.ident),]
my_sample_col=as.data.frame(pbmc4@meta.data$active.ident)
rownames(my_sample_col)=rownames(pbmc4@meta.data)
colnames(my_sample_col)=c("Cluster")
Var2 =brewer.pal(3,"Set2")
names(Var2) = c("Primary","Lung","Brain")
ann_colors = list(Cluster=Var2)
my_sample_col$Cluster=factor(my_sample_col$Cluster, levels=c("Primary","Lung","Brain"))
write_out= as.data.frame(as.matrix(pbmc4@assays$RNA@scale.data[,rownames(pbmc4@meta.data)]))
write_out=write_out[pbmc@assays$RNA@var.features,]
write_out2=round(cor(write_out),digits=2)
pheatmap(write_out2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-0.3, 0.3, by = 0.006),annotation_col =my_sample_col,annotation_row =my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100), legend = F, annotation_legend = F)

#primary correlation average 0.064
dat.n3=write_out2[1:1000,1:1000]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)
#Lung correlation average 0.038
dat.n3=write_out2[1001:2000,1001:2000]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)
#Brain correlation average 0.028
dat.n3=write_out2[2001:3000,2001:3000]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)
##Others
mean(write_out2[1001:2000,1:1000])
mean(write_out2[2001:3000,1:1000])
mean(write_out2[2001:3000,1001:2000])
rm(pbmc4,write_out,write_out2)


###Heatmap showing expression of DEGs between primary and lung/brain DEG 
###Heatmap showing expression of DEGs between primary and lung/brain DEG
###Heatmap showing expression of DEGs between primary and lung/brain DEG
pbmc4=subset(x = pbmc, downsample = 1000)
pbmc4 <- ScaleData(pbmc4,features=rownames(pbmc4))
pbmc4@meta.data$active.ident=pbmc4@active.ident
pbmc4@meta.data=pbmc4@meta.data[sample(1:nrow(pbmc4@meta.data)),]
pbmc4@meta.data=pbmc4@meta.data[order(pbmc4@meta.data$active.ident),]
my_sample_col=as.data.frame(pbmc4@meta.data$active.ident)
colnames(my_sample_col)=c("Cluster")
rownames(my_sample_col)=rownames(pbmc4@meta.data)
Var2 =brewer.pal(3,"Set2")
names(Var2) = c("Primary","Lung","Brain")
ann_colors = list(Cluster=Var2)
my_sample_col$Cluster=factor(my_sample_col$Cluster, levels=c("Primary","Lung","Brain"))
write_out= as.data.frame(as.matrix(pbmc4@assays$RNA@scale.data[,rownames(pbmc4@meta.data)]))
meta.markers2=meta.markers[order(meta.markers$avg_log2FC, decreasing =F),]
write_out2=write_out[rownames(meta.markers2),]
table(meta.markers$avg_log2FC>0)
DEG=as.data.frame(colSums(write_out2[1:37,])/37-colSums(write_out2[38:69,])/32)
colnames(DEG)="DEG"
write_out3=write_out2[,order(DEG$DEG, decreasing=T)] 
##550/300
pheatmap(write_out3,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-2, 2, by = 0.04),annotation_col =my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100), legend = F, annotation_legend = F)
rm(pbmc4,write_out,write_out2,write_out3,meta.markers2)
dim(meta.markers)

###Heatmap showing expression of DEGs between lung and brain
###Heatmap showing expression of DEGs between lung and brain
###Heatmap showing expression of DEGs between lung and brain
pbmc4=subset(x = pbmc, downsample = 1000)
pbmc4 <- ScaleData(pbmc4,features=rownames(pbmc4))
pbmc4@meta.data$active.ident=pbmc4@active.ident
pbmc4@meta.data=pbmc4@meta.data[sample(1:nrow(pbmc4@meta.data)),]
pbmc4@meta.data=pbmc4@meta.data[order(pbmc4@meta.data$active.ident),]
my_sample_col=as.data.frame(pbmc4@meta.data$active.ident)
colnames(my_sample_col)=c("Cluster")
rownames(my_sample_col)=rownames(pbmc4@meta.data)
Var2 =brewer.pal(3,"Set2")
names(Var2) = c("Primary","Lung","Brain")
ann_colors = list(Cluster=Var2)
my_sample_col$Cluster=factor(my_sample_col$Cluster, levels=c("Primary","Lung","Brain"))
write_out= as.data.frame(as.matrix(pbmc4@assays$RNA@scale.data[,rownames(pbmc4@meta.data)]))
meta.markers2=brain.markers[order(brain.markers$avg_log2FC, decreasing =F),]
write_out2=write_out[rownames(meta.markers2),]
table(brain.markers$avg_log2FC>0)
DEG=as.data.frame(colSums(write_out2[1:15,])/15 - colSums(write_out2[16:36,])/21)
colnames(DEG)="DEG"
write_out3=write_out2[,order(DEG$DEG, decreasing=T)] 
##550/300
pheatmap(write_out3,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-2, 2, by = 0.04),annotation_col =my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100), legend = F, annotation_legend = F)
rm(pbmc4,write_out,write_out2,write_out3,meta.markers2)
dim(brain.markers)



##PC2 density plot
pbmc@meta.data$PC1=pbmc@reductions$pca@cell.embeddings[,1]
pbmc@meta.data$PC2=pbmc@reductions$pca@cell.embeddings[,2]
pbmc@meta.data$RNA_snn_res.0.1=factor(pbmc@meta.data$RNA_snn_res.0.1,levels=c("Primary","Lung","Brain"))
ggplot(pbmc@meta.data, aes(x=PC1, fill=RNA_snn_res.0.1, color=RNA_snn_res.0.1)) +geom_density(alpha=0.3)+  theme_classic()+theme(text = element_text(size=12))+theme(axis.title.x=element_text(color="black", size=12, face="bold"),axis.title.y=element_text(color="black", size=12, face="bold"),axis.text.x=element_text(color="black", size=12, face="bold"),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"), strip.text = element_text(size = 12, face="bold"))+ylab("Denstiy")+xlab("PC1")+scale_color_manual(values=brewer.pal(3,"Set2"))+scale_fill_manual(values=brewer.pal(3,"Set2")) + theme(legend.position = "none")
ggplot(pbmc@meta.data, aes(x=PC2, fill=RNA_snn_res.0.1, color=RNA_snn_res.0.1)) +geom_density(alpha=0.3)+  theme_classic()+theme(text = element_text(size=12))+theme(axis.title.x=element_text(color="black", size=14, face="bold"),axis.title.y=element_text(color="black", size=14, face="bold"),axis.text.x=element_text(color="black", size=12, face="bold"),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"), strip.text = element_text(size = 12, face="bold"))+ylab("Denstiy")+xlab("PC2")+scale_color_manual(values=brewer.pal(3,"Set2"))+scale_fill_manual(values=brewer.pal(3,"Set2")) + theme(legend.position = "none")

merge.markers <- FindMarkers(object = pbmc, ident.1=c("Primary"),ident.2=c("Lung","Brain"),only.pos = FALSE, group.by="RNA_snn_res.0.1")
write.csv(x = merge.markers, row.names = TRUE, file = "DEG.markers.csv")

rownames(merge.markers)
###DEG / correlation plot 
###DEG / correlation plot 
###DEG / correlation plot 
my_sample_col=as.data.frame(cbind(as.data.frame(pbmc@meta.data$orig.ident),as.data.frame(pbmc@meta.data$RNA_snn_res.0.1)))
colnames(my_sample_col)=c("Type","Cluster")
rownames(my_sample_col)=rownames(pbmc@meta.data)
Var1 =brewer.pal(4,"Set1")
names(Var1) = c("WM239A","113/6-4L", "131/4-5B1", "131/4-5B2")
Var2 =brewer.pal(3,"Set2")
names(Var2) = c("Primary","Lung","Brain")
ann_colors = list(Type = Var1, Cluster=Var2)
my_sample_col$Cluster=factor(my_sample_col$Cluster, levels=c("Primary","Lung","Brain"))
my_sample_col2=my_sample_col[order(my_sample_col$Cluster),]
head(my_sample_col2)
write_out= as.data.frame(as.matrix(pbmc@assays$RNA@scale.data))
merge.markers=merge.markers[order(merge.markers$avg_logFC, decreasing =TRUE),]
write_out2=write_out[rownames(merge.markers),order(my_sample_col$Cluster)]
rm(write_out)
dat.n=round(write_out2,digits=2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-0.5, 0.5, by = 0.01),annotation_col =my_sample_col2,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))
dim(dat.n)

DEG=as.data.frame(colSums(dat.n[1:90,])/90 - colSums(dat.n[91:159,])/69)
colnames(DEG)="DEG"

dat.n2=dat.n[,order(DEG$DEG, decreasing=TRUE)] 
##1050/400
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-0.5, 0.5, by = 0.01),annotation_col =my_sample_col2,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))
dim(dat.n)

##Correlation plot 
write_out= as.data.frame(as.matrix(pbmc@assays$RNA@scale.data))
write_out2=write_out[pbmc@assays$RNA@var.features,]
rm(write_out)
write_out3=write_out2[,order(my_sample_col$Cluster)]
cor_dat.n=round(cor(write_out3),digits=2)

pheatmap(cor_dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-0.3, 0.3, by = 0.006),annotation_col =my_sample_col2,annotation_row =my_sample_col2,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))
table(pbmc@active.ident)
#primary correlation average 0.065
dat.n3=cor_dat.n[1:2167,1:2167]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)
table(pbmc@meta.data$RNA_snn_res.0.1)
#Lung correlation average 0.038
dat.n3=cor_dat.n[2168:4002,2168:4002]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)

#Brain correlation average 0.028
dat.n3=cor_dat.n[4003:8173,4003:8173]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)

##Others
mean(cor_dat.n[2168:4002,1:2167])
mean(cor_dat.n[4003:8173,1:2167])
mean(cor_dat.n[2168:4002,4003:8173])






dat.n3=cor_dat.n[2583:4156,2583:4156]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)








##PC2 density plot
pbmc@meta.data$PC1=pbmc@reductions$pca@cell.embeddings[,1]
pbmc@meta.data$PC2=pbmc@reductions$pca@cell.embeddings[,2]

pbmc@meta.data$RNA_snn_res.0.1=factor(pbmc@meta.data$RNA_snn_res.0.1,levels=c("Primary","Lung","Brain"))
ggplot(pbmc@meta.data, aes(x=PC1, fill=RNA_snn_res.0.1, color=RNA_snn_res.0.1)) +geom_density(alpha=0.3)+  theme_classic()+theme(text = element_text(size=12))+theme(axis.title.x=element_text(color="black", size=12, face="bold"),axis.title.y=element_text(color="black", size=12, face="bold"),axis.text.x=element_text(color="black", size=12, face="bold"),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"), strip.text = element_text(size = 12, face="bold"))+ylab("Denstiy")+xlab("PC1")+scale_color_manual(values=brewer.pal(3,"Set2"))+scale_fill_manual(values=brewer.pal(3,"Set2")) + theme(legend.position = "none")
ggplot(pbmc@meta.data, aes(x=PC2, fill=RNA_snn_res.0.1, color=RNA_snn_res.0.1)) +geom_density(alpha=0.3)+  theme_classic()+theme(text = element_text(size=12))+theme(axis.title.x=element_text(color="black", size=14, face="bold"),axis.title.y=element_text(color="black", size=14, face="bold"),axis.text.x=element_text(color="black", size=12, face="bold"),plot.background=element_blank(),axis.line = element_line(colour="black"),axis.text.y = element_text(color="black", size=12, face="bold"), strip.text = element_text(size = 12, face="bold"))+ylab("Denstiy")+xlab("PC2")+scale_color_manual(values=brewer.pal(3,"Set2"))+scale_fill_manual(values=brewer.pal(3,"Set2")) + theme(legend.position = "none")



merge.markers2 <- FindMarkers(object = pbmc, ident.1=c("Lung"),ident.2=c("Brain"),only.pos = FALSE, group.by="RNA_snn_res.0.1")
write.csv(x = merge.markers2, row.names = TRUE, file = "DEG.markers2.csv")




###DEG / correlation plot 
###DEG / correlation plot 
###DEG / correlation plot 
my_sample_col=as.data.frame(cbind(as.data.frame(pbmc@meta.data$orig.ident),as.data.frame(pbmc@meta.data$RNA_snn_res.0.1)))
colnames(my_sample_col)=c("Type","Cluster")
rownames(my_sample_col)=rownames(pbmc@meta.data)
Var1 =brewer.pal(4,"Set1")
names(Var1) = c("WM239A","113/6-4L", "131/4-5B1", "131/4-5B2")
Var2 =brewer.pal(3,"Set2")
names(Var2) = c("Primary","Lung","Brain")
ann_colors = list(Type = Var1, Cluster=Var2)
my_sample_col$Cluster=factor(my_sample_col$Cluster, levels=c("Primary","Lung","Brain"))
my_sample_col2=my_sample_col[order(my_sample_col$Cluster),]
head(my_sample_col2)
write_out= as.data.frame(as.matrix(pbmc@assays$RNA@scale.data))
merge.markers2=merge.markers2[order(merge.markers2$avg_logFC, decreasing =TRUE),]
write_out2=write_out[rownames(merge.markers2),order(my_sample_col$Cluster)]
rm(write_out)
dat.n=round(write_out2,digits=2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-0.5, 0.5, by = 0.01),annotation_col =my_sample_col2,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))
rownames(merge.markers2)

DEG=as.data.frame(colSums(dat.n[1:40,])/40 - colSums(dat.n[41:81,])/41)
colnames(DEG)="DEG"

dat.n2=dat.n[,order(DEG$DEG, decreasing=TRUE)] 
##1050/400
pheatmap(dat.n2,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-0.5, 0.5, by = 0.01),annotation_col =my_sample_col2,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))
dim(dat.n)




###DEG Visualization by Dotplot
###DEG Visualization by Dotplot
###DEG Visualization by Dotplot

##Primary up-regulated
markers.to.plot=c("FABP7", "MMP1", "BCAS3", "MFSD12", "IGFBP7", "MGST1", "S100B", "CAPG", "PRNP", "FMN2","COTL1", "CD74", "SLC16A3", "NQO1", "FN1", "MT2A","ANXA1")


##Lung & Brain up-
markers.to.plot=c("TTC39A", "IGFBP2", "PMEL", "TMSB4X", "CCND1", "VCX3B", "APOD", "BST2", "APOC1", "APOE", "PYROXD2", "PAGE5")
markers.to.plot=c("COL9A3", "S100A6", "PRDX1", "ID3", "AC147651.4", "VGF", "SDCBP", "DCT", "NUPR1")


##Selected between primary and meta 
markers.to.plot1=c("CD74", "IGFBP7","CCND1","TMSB4X", "PMEL","APOE","SLC16A3","FN1","MT2A")

pbmc@active.ident=factor(pbmc@active.ident,levels=c("Primary","Lung","Brain"))

##Preferential lung/brain 
markers.to.plot2=c("S100A6", "PRDX1",  "SDCBP", "DCT", "NUPR1")

##800/300RotatedAxis()
?DotPlot
DotPlot(pbmc, features = markers.to.plot1, cols = c("grey", "red"),dot.min=0)+theme(axis.title=element_blank(), axis.text.x.bottom = element_text(size=12,face="bold"), axis.text.y.left=element_text(size = 12,face="bold")) +coord_flip()
DotPlot(pbmc, features = markers.to.plot2, cols = c("grey", "red"),dot.min=0)+theme(axis.title=element_blank(), axis.text.x.bottom = element_text(size=12,face="bold"), axis.text.y.left=element_text(size = 12,face="bold")) +coord_flip()


s.gene <- c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8")
g2m.gene <- c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA")
pbmc <- CellCycleScoring(pbmc, s.features = s.gene, g2m.features = g2m.gene, set.ident = FALSE)

MITF.score=list(c("MITF", "TYR", "PMEL", "PLP1", "GPR143", "MLANA", "STX7", "IRF4", "ERBB3", "CDH1", "GPNMB", "IGSF11", "SLC24A5", "SLC45A2", "RAP2B", "ASAH1", "MYO10", "GRN", "DOCK10", "ACSL3", "SORT1", "QPCT", "S100B", "MYC", "LZTS1", "GYG2", "SDCBP", "LOXL4", "ETV5", "C1orf85", "HMCN1", "OSTM1", "ALDH7A1", "FOSB", "RAB38", "ELOVL2", "MLPH", "PLK2", "CHL1", "RDH11", "LINC00473", "RELL1", "C21orf91", "SCAMP3", "SGK3", "ABCB5", "SLC7A5", "SIRPA", "WDR91", "PIGS", "CYP27A1", "TM7SF3", "PTPRZ1", "CNDP2", "CTSK", "BNC2", "TOB1", "CELF2", "ROPN1", "TMEM98", "CTSA", "LIMA1", "CD99", "IGSF8", "FDFT1", "CPNE3", "SLC35B4", "EIF3E", "TNFRSF14", "VAT1", "HPS5", "CDK2", "CAPN3", "SUSD5", "ADSL", "PIGY", "PON2", "SLC19A1", "KLF6", "MAGED1", "ERGIC3", "PIR", "SLC25A5", "JUN", "ARPC1B", "SLC19A2", "AKR7A2", "HPGD", "TBC1D7", "TFAP2A", "PTPLAD1", "SNCA", "GNPTAB", "DNAJA4", "APOE", "MTMR2", "ATP6V1B2", "C16orf62", "EXOSC4", "STAM"))
AXL.score=list(c("ANGPTL4", "FSTL3", "GPC1", "TMSB10", "SH3BGRL3", "PLAUR", "NGFR", "SEC14L2", "FOSL1", "SERPINE1", "IGFBP3", "TNFRSF12A", "GBE1", "AXL", "PHLDA2", "MAP1B", "GEM", "SLC22A4", "TYMP", "TREM1", "RIN1", "S100A4", "COL6A2", "FAM46A", "CITED1", "S100A10", "UCN2", "SPHK1", "TRIML2", "S100A6", "TMEM45A", "CDKN1A", "UBE2C", "ERO1L", "SLC16A6", "CHI3L1", "FN1", "S100A16", "CRIP1", "SLC25A37", "LCN2", "ENO2", "PFKFB4", "SLC16A3", "DBNDD2", "LOXL2", "CFB", "CADM1", "LTBP3", "CD109", "AIM2", "TCN1", "STRA6", "C9orf89", "DDR1", "TBC1D8", "METTL7B", "GADD45A", "UPP1", "SPATA13", "GLRX", "PPFIBP1", "PMAIP1", "COL6A1", "JMJD6", "CIB1", "HPCAL1", "MT2A", "ZCCHC6", "IL8", "TRIM47", "SESN2", "PVRL2", "DRAP1", "MTHFD2", "SDC4", "NNMT", "PPL", "TIMP1", "RHOC", "GNB2", "PDXK", "CTNNA1", "CD52", "SLC2A1", "BACH1", "ARHGEF2", "UBE2J1", "CD82", "ZYX", "P4HA2", "PEA15", "GLRX2", "HAPLN3", "RAB36", "SOD2", "ESYT2", "IL18BP", "FGFRL1", "PLEC"))

s.g2m.gene = list(c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8","HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA"))
pbmc <- AddModuleScore(object = pbmc, features = MITF.score, name = "MITF.score")
pbmc <- AddModuleScore(object = pbmc, features = AXL.score, name = "AXL.score")
pbmc <- AddModuleScore(object = pbmc, features = s.g2m.gene, name = "s.g2m.score")


PC2_high=list(c("BST2", "PAGE5", "APOE", "PYROXD2", "PMEL", "LGALS3", "IGFBP2", "APOC1", "CCND1", "NBL1", "CASP1", "EMP3", "APOD", "VCX3B", "TTC39A", "BCL2A1", "SFRP1", "S100A11", "HLA-C", "HLA-A", "APOC2", "CPN1", "GYPC", "VIM", "GSTO1", "MSRB2", "SIGIRR", "ST3GAL6", "VCX", "PPAP2C", "GTSF1", "NUPR1", "MLANA", "GPR137B", "UGCG", "LINC00221", "DLL3", "ANXA6", "CPVL", "MCOLN3", "SEPT6", "HLA-B", "SPSB1", "AC145110.1", "NELL1", "FAM167B", "ANXA3", "QPRT", "MAGEA1", "PYCARD", "CDH3", "GPR143", "MYEOV", "RP11-401O9.3", "GLTSCR2", "SNHG18", "NME4", "PXDN", "QPCT", "TM4SF19.1", "PLP2", "STK32A", "RGS20", "TUBA3C", "TIMM50", "CABLES1", "HES6", "FBLN1", "RAB20", "RAMP2", "DCT", "PSMB8", "BHLHE41", "BCAN", "BAIAP2L2", "TMEM220", "RLBP1", "SWAP70", "ARHGAP29", "FAM26F", "HTRA1", "DMKN", "CT55", "TMSB4X", "TNFSF12", "GSG1L", "NPAS2", "CNIH3", "RPL37", "ZNF503", "SLC27A5", "JUN", "GRN", "RP11-408E5.4", "MPPED2", "ZNF703", "C10orf11", "PRSS33", "LEPREL1", "AC104655.3"))
PC2_low=list(c("AGMO", "RND3", "CHCHD3", "G6PD", "KLF9", "CBR1", "ANXA2", "TMEM106C", "PKM", "ITGA4", "EVI2A", "FN1", "TNC", "RP11-51M18.1", "ID1", "HES1", "SH2B3", "MMP8", "TJP1", "ANKS1A", "PFN2", "LY96", "LXN", "SLC20A1", "SNTB1", "SRPX", "NAA38", "CTD-3049M7.1", "HSP90AA1", "PMP2", "HLA-DPA1", "PFKFB4", "HSPE1", "HLA-DMA", "ATP6V0E2", "FABP7", "DHRS3", "LUM", "TAF7", "EXTL1", "DTYMK", "RDH5", "PLP1", "MCAM", "PPP2R4", "RASGRF1", "TXNRD1", "CADM4", "SLC30A10", "ATP2B1", "ACSS3", "CHCHD2", "TGFA", "CCT6A", "PDIA4", "C3orf67", "RP11-574H6.1", "CTSD", "RP11-19N8.4", "CDH2", "C1QTNF3", "TIMP3", "SPP1", "PDLIM2", "MFI2", "POPDC3", "KCNN4", "HLA-DRA", "NQO1", "MGST1", "MTUS1", "LSM5", "CTHRC1", "S100A6", "CSRP1", "CAPG", "SH3BGRL3", "SVIP", "NGFRAP1", "COPS8", "TCTEX1D2", "MMP1", "SLC16A3", "MFSD12", "ZNF706", "CD74", "CACNA2D4", "OAF", "HMG20B", "BCHE", "MPZ", "IGFBP7", "S100B", "ANXA1", "FMN2", "MT2A", "PRNP", "COTL1", "AC147651.4", "BCAS3"))
PC2=c("BST2", "PAGE5", "APOE", "PYROXD2", "PMEL", "LGALS3", "IGFBP2", "APOC1", "CCND1", "NBL1", "CASP1", "EMP3", "APOD", "VCX3B", "TTC39A", "BCL2A1", "SFRP1", "S100A11", "HLA-C", "HLA-A", "APOC2", "CPN1", "GYPC", "VIM", "GSTO1", "MSRB2", "SIGIRR", "ST3GAL6", "VCX", "PPAP2C", "GTSF1", "NUPR1", "MLANA", "GPR137B", "UGCG", "LINC00221", "DLL3", "ANXA6", "CPVL", "MCOLN3", "SEPT6", "HLA-B", "SPSB1", "AC145110.1", "NELL1", "FAM167B", "ANXA3", "QPRT", "MAGEA1", "PYCARD", "CDH3", "GPR143", "MYEOV", "RP11-401O9.3", "GLTSCR2", "SNHG18", "NME4", "PXDN", "QPCT", "TM4SF19.1", "PLP2", "STK32A", "RGS20", "TUBA3C", "TIMM50", "CABLES1", "HES6", "FBLN1", "RAB20", "RAMP2", "DCT", "PSMB8", "BHLHE41", "BCAN", "BAIAP2L2", "TMEM220", "RLBP1", "SWAP70", "ARHGAP29", "FAM26F", "HTRA1", "DMKN", "CT55", "TMSB4X", "TNFSF12", "GSG1L", "NPAS2", "CNIH3", "RPL37", "ZNF503", "SLC27A5", "JUN", "GRN", "RP11-408E5.4", "MPPED2", "ZNF703", "C10orf11", "PRSS33", "LEPREL1", "AC104655.3","AGMO", "RND3", "CHCHD3", "G6PD", "KLF9", "CBR1", "ANXA2", "TMEM106C", "PKM", "ITGA4", "EVI2A", "FN1", "TNC", "RP11-51M18.1", "ID1", "HES1", "SH2B3", "MMP8", "TJP1", "ANKS1A", "PFN2", "LY96", "LXN", "SLC20A1", "SNTB1", "SRPX", "NAA38", "CTD-3049M7.1", "HSP90AA1", "PMP2", "HLA-DPA1", "PFKFB4", "HSPE1", "HLA-DMA", "ATP6V0E2", "FABP7", "DHRS3", "LUM", "TAF7", "EXTL1", "DTYMK", "RDH5", "PLP1", "MCAM", "PPP2R4", "RASGRF1", "TXNRD1", "CADM4", "SLC30A10", "ATP2B1", "ACSS3", "CHCHD2", "TGFA", "CCT6A", "PDIA4", "C3orf67", "RP11-574H6.1", "CTSD", "RP11-19N8.4", "CDH2", "C1QTNF3", "TIMP3", "SPP1", "PDLIM2", "MFI2", "POPDC3", "KCNN4", "HLA-DRA", "NQO1", "MGST1", "MTUS1", "LSM5", "CTHRC1", "S100A6", "CSRP1", "CAPG", "SH3BGRL3", "SVIP", "NGFRAP1", "COPS8", "TCTEX1D2", "MMP1", "SLC16A3", "MFSD12", "ZNF706", "CD74", "CACNA2D4", "OAF", "HMG20B", "BCHE", "MPZ", "IGFBP7", "S100B", "ANXA1", "FMN2", "MT2A", "PRNP", "COTL1", "AC147651.4", "BCAS3")

PC1_low=c("CXXC5", "LY6E", "BTG1", "EGR3", "HSPB1", "SGK1", "BRI3", "C17orf89", "AC074389.9", "SNTB1", "BNIP3L", "LINC01003", "L1CAM", "RP11-51M18.1", "MMP1", "CCDC85B", "SLC16A3", "TIMP3", "HLA-DPA1", "MTUS1", "G6PD", "MYC", "LARP6", "RP3-395M20.12", "TP53TG1", "KIAA0930", "IER3", "PKM", "MBP", "MFI2", "ZFAS1", "CYGB", "RDH5", "PHLDA3", "MPZ", "C6orf48", "CACNA2D4", "HLA-DRA", "UBC", "CITED1", "EXTL1", "PCBP4", "RAB13", "CHCHD2", "BCHE", "CD27-AS1", "GDF15", "COPZ2", "KLF9", "FTL", "SPP1", "PLAT", "WIPI1", "PRSS23", "VIM", "VAMP5", "RAB40B", "IGFBP7", "RPL37", "RP11-19N8.4", "SQSTM1", "CD74", "TPST1", "METRN", "UBB", "BCAS3", "ALDOA", "OAF", "CD44", "ZMAT3", "PDLIM2", "GS1-124K5.4", "MFSD12", "FN1", "RPS12", "TIMP1", "RPL38", "SH3BGRL3", "CTSD", "RND3", "PNRC1", "MXD4", "IGFBP5", "S100A6", "FMN2", "RP4-718J7.4", "S100B", "RPL23A", "TKT", "GPX1", "GPR56", "PHLDA1", "PLP1", "RPL37A", "NEAT1", "RPL31", "RPL29", "RPL32", "RPL7A", "RPL30")
pbmc <- AddModuleScore(object = pbmc, features = PC2_high, name = "PC2_high")
pbmc <- AddModuleScore(object = pbmc, features = PC2_low, name = "PC2_low")

pbmc@meta.data$PC2.score=pbmc@meta.data$PC2_high1-pbmc@meta.data$PC2_low1

###PC1 heatmap
###PC1 heatmap
###PC1 heatmap
my_sample_col=as.data.frame(cbind(as.data.frame(pbmc@meta.data$orig.ident),as.data.frame(pbmc@meta.data$RNA_snn_res.0.1)))
colnames(my_sample_col)=c("Type","Cluster")
rownames(my_sample_col)=rownames(pbmc@meta.data)
Var1 =brewer.pal(4,"Set1")
names(Var1) = c("WM239A","113/6-4L", "131/4-5B1", "131/4-5B2")
Var2 =brewer.pal(3,"Set2")
names(Var2) = c("Primary","Lung","Brain")
ann_colors = list(Type = Var1, Cluster=Var2)
my_sample_col$Cluster=factor(my_sample_col$Cluster, levels=c("Primary","Lung","Brain"))
write_out= as.data.frame(as.matrix(pbmc@assays$RNA@scale.data))
write_out2=write_out[,order(pbmc@meta.data$PC1)]
rm(write_out)
write_out3=as.data.frame(rbind(write_out2[PC1,],write_out2[PC1_low,]))

#1050/400
pheatmap(write_out3,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-0.5, 0.5, by = 0.01),annotation_col =my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))




###PC2 heatmap
###PC2 heatmap
###PC2 heatmap
my_sample_col=as.data.frame(cbind(as.data.frame(pbmc@meta.data$orig.ident),as.data.frame(pbmc@meta.data$RNA_snn_res.0.1)))
colnames(my_sample_col)=c("Type","Cluster")
rownames(my_sample_col)=rownames(pbmc@meta.data)
Var1 =brewer.pal(4,"Set1")
names(Var1) = c("WM239A","113/6-4L", "131/4-5B1", "131/4-5B2")
Var2 =brewer.pal(3,"Set2")
names(Var2) = c("Primary","Lung","Brain")
ann_colors = list(Type = Var1, Cluster=Var2)
my_sample_col$Cluster=factor(my_sample_col$Cluster, levels=c("Primary","Lung","Brain"))
key.markers=PC2
write_out= as.data.frame(as.matrix(pbmc@assays$RNA@scale.data))
write_out2=write_out[key.markers,order(pbmc@meta.data$PC2)]
rm(write_out)
dat.n=round(write_out2,digits=2)
dim(dat.n)
#1050/400
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-0.5, 0.5, by = 0.01),annotation_col =my_sample_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))


##Cell-to-cell Correlation plot 
write_out= as.data.frame(as.matrix(pbmc@assays$RNA@scale.data))
write_out2=write_out[pbmc@assays$RNA@var.features,]
rm(write_out)
write_out3=write_out2[,order(my_sample_col$Cluster)]
dat.n=round(cor(write_out3),digits=2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-0.3, 0.3, by = 0.006),annotation_col =my_sample_col,annotation_row =my_sample_col2,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))


##Gene-to-gene Correlation plot 
write_out= as.data.frame(as.matrix(pbmc@assays$RNA@scale.data))
write_out2=rbind(write_out[PC1,],write_out[PC2_high1,],write_out[PC2_low1,])
rm(write_out)
dat.n=round(cor(t(write_out2)),digits=2)

my_gene_col=as.data.frame(colnames(dat.n))
colnames(my_gene_col)="Genes"
rownames(my_gene_col)=my_gene_col$Genes
my_gene_col$Genes=rep(c("PC1_high","PC2_high","PC2_low"),c(100,100,100))

Var2 = brewer.pal(3,"Accent")
names(Var2) = c("PC1_high","PC2_high","PC2_low")
ann_colors = list(Genes=Var2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-0.5, 0.5, by = 0.01),annotation_col =my_gene_col,annotation_row =my_gene_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),annotation_color=ann_colors,legend=FALSE, annotation_legend	=FALSE)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-0.5, 0.5, by = 0.01),annotation_col =my_gene_col,annotation_row =my_gene_col,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors,color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),annotation_color=ann_colors)

#PC1 high correlation average 0.367
dat.n3=dat.n[1:100,1:100]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)
table(pbmc@meta.data$RNA_snn_res.0.1)
#PC2 high correlation average 0.155
dat.n3=dat.n[101:200,101:200]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)

#PC2 low correlation average 0.180
dat.n3=dat.n[201:300,201:300]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)

#-0.010
dat.n3=dat.n[1:100,101:200]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)

#-0.005
dat.n3=dat.n[1:100,201:300]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)

#-0.151
dat.n3=dat.n[101:200,201:300]
dat.n3[dat.n3==1] <- NA
mean(dat.n3,na.rm=TRUE)



###Expression of pathway 
###Expression of pathway 
###Expression of pathway 

FeaturePlot(object = pbmc, reduction="pca", features = c("S.Score"))+ ggtitle("Average S phase related gene expression") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+xlab("PC1 (14.2% of variance explained)")+ylab("PC2 (13.4% of variance explained)")
FeaturePlot(object = pbmc, reduction="pca", features = c("G2M.Score"))+ ggtitle("Average G2/M phase related gene expression") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+xlab("PC1 (14.2% of variance explained)")+ylab("PC2 (13.4% of variance explained)")


FeaturePlot(object = pbmc, features = c("S.Score"))+ ggtitle("Average S phase related gene expression") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

FeaturePlot(object = pbmc, features = c("G2M.Score"))+ ggtitle("Average G2/M phase related gene expression") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

FeaturePlot(object = pbmc, features = c("MITF.score1"))+ ggtitle("Average MITF pathway related gene expression") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

FeaturePlot(object = pbmc, features = c("AXL.score1"))+ ggtitle("Average AXL pathway related gene expression") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

FeaturePlot(object = pbmc, features = c("PC2.score"))+ ggtitle("Average PC2 expression") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

FeaturePlot(object = pbmc, features = c("PC2"))+ ggtitle("Average PC2 high (100 genes) expression") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))


FeaturePlot(object = pbmc, features = c("PC2_high1"))+ ggtitle("Average PC2 high (100 genes) expression") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

FeaturePlot(object = pbmc, features = c("PC2_low1"))+ ggtitle("Average PC2 low (100 genes) expression") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

FeaturePlot(object = pbmc, features = c("Pseudotime"))+ ggtitle("Monocle Pseudotime") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))
pbmc@meta.data$s.g2m.score1=pbmc@meta.data$S.Score+pbmc@meta.data$G2M.Score
pbmc@meta.data$s.g2m.score1=pbmc@meta.data$s.g2m.score1/2
FeaturePlot(object = pbmc, features = c("s.g2m.score1"))+ ggtitle("Average S or G2/M phase related gene expression") + theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+ scale_colour_gradientn(colours=brewer.pal(9,"OrRd")[1:9])+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))


##Correlation check - FeatureScatter
##Correlation check - FeatureScatter
##Correlation check - FeatureScatter
FeatureScatter(object = pbmc, feature1 = "PC1", feature2 = "S.Score")+ theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+  theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

FeatureScatter(object = pbmc, feature1 = "PC1", feature2 = "G2M.Score")+ theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+  theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))


rect=data.frame(xmin=-Inf, xmax=0, ymin=-Inf, ymax=0)

rect2 <- data.frame(xmin=0, xmax=Inf, ymin=-Inf, ymax=Inf)
rect3 <- data.frame(xmin=-Inf, xmax=0, ymin=0, ymax=Inf)
?FeatureScatter
FeatureScatter(object = pbmc, feature1 = "G2M.Score", feature2 = "S.Score",pt.size=0.7, cells	= rownames(pbmc@meta.data[pbmc@meta.data$RNA_snn_res.0.1=="Brain",]))+theme(plot.title = element_blank(), axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.1, inherit.aes = FALSE)+geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE)+geom_rect(data=rect3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE)+xlab("Average expression of G2/M phase genes")+ylab("Average expression of S phase genes")+NoLegend()+  scale_color_manual(values="#8DA0CB")

FeatureScatter(object = pbmc, feature1 = "G2M.Score", feature2 = "S.Score",pt.size=0.7, cells	= rownames(pbmc@meta.data[pbmc@meta.data$RNA_snn_res.0.1=="Lung",]))+theme(plot.title = element_blank(), axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.1, inherit.aes = FALSE)+geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE)+geom_rect(data=rect3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE)+xlab("Average expression of G2/M phase genes")+ylab("Average expression of S phase genes")+NoLegend()+  scale_color_manual(values="#FC8D62")

FeatureScatter(object = pbmc, feature1 = "G2M.Score", feature2 = "S.Score",pt.size=0.7, cells	= rownames(pbmc@meta.data[pbmc@meta.data$RNA_snn_res.0.1=="Primary",]))+theme(plot.title = element_blank(), axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.1, inherit.aes = FALSE)+geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE)+geom_rect(data=rect3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE)+xlab("Average expression of G2/M phase genes")+ylab("Average expression of S phase genes")+NoLegend()+  scale_color_manual(values="#66C2A5")


FeatureScatter(object = pbmc, feature1 = "G2M.Score", feature2 = "S.Score",pt.size=0.7, cells	= rownames(pbmc@meta.data[pbmc@meta.data$RNA_snn_res.0.1=="Brain",]))+theme(plot.title = element_blank(), axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.1, inherit.aes = FALSE)+geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE)+geom_rect(data=rect3, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="red", alpha=0.1, inherit.aes = FALSE)+xlab("Average expression of G2/M phase genes")+ylab("Average expression of S phase genes")+NoLegend()



?FeatureScatter

FeatureScatter(object = pbmc, feature1 = "AXL.score1", feature2 = "MITF.score1", group.by="orig.ident")

FeatureScatter(object = pbmc, feature1 = "PC2_high1", feature2 = "Pseudotime")+ theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+  theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

FeatureScatter(object = pbmc,feature1 = "PC_2", feature2 = "MITF.score1",cols=brewer.pal(3,"Set2"))+ theme(plot.title = element_blank())+  theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+xlab("PC2")+ylab("Average expression of MITF pathway")+NoLegend()
cor.test(pbmc@meta.data$PC2, pbmc@meta.data$MITF.score1)

FeatureScatter(object = pbmc,feature1 = "PC_2", feature2 = "AXL.score1",cols=brewer.pal(3,"Set2"))+ theme(plot.title = element_blank())+  theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+xlab("PC2")+ylab("Average expression of AXL pathway")+NoLegend()
cor.test(pbmc@meta.data$PC2, pbmc@meta.data$AXL.score1)
colnames(pbmc@meta.data)
cor.test(pbmc@meta.data[pbmc@meta.data$RNA_snn_res.0.1 %in% c("Lung","Brain"),]$PC2.score, pbmc@meta.data[pbmc@meta.data$RNA_snn_res.0.1 %in% c("Lung","Brain"),]$MITF.score1)
FeatureScatter(object = pbmc,cells=rownames(pbmc@meta.data[pbmc@meta.data$RNA_snn_res.0.1 %in% c("Lung","Brain"),]), feature1 = "PC2.score", feature2 = "MITF.score1")+ theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+  theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))




FeatureScatter(object = pbmc,feature1 = "PC2.score", feature2 = "MITF.score1")+ theme(plot.title = element_blank())+  theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+ylab("PC2 signature")+d
cor.test(pbmc@meta.data$PC2.score, pbmc@meta.data$MITF.score1)


FeatureScatter(object = pbmc,cells=rownames(pbmc@meta.data[pbmc@meta.data$RNA_snn_res.0.1 %in% c("Primary"),]), feature1 = "PC2.score", feature2 = "MITF.score1")+ theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+  theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))




?FeatureScatter
FeatureScatter(object = pbmc, feature1 = "PC2.score", feature2 = "MITF.score1")+ theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+  theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

FeatureScatter(object = pbmc, feature1 = "PC2", feature2 = "Pseudotime",cols =brewer.pal(3,"Set2") )+ theme(plot.title = element_blank(),axis.title = element_text(size=16,face="bold"), axis.text = element_text(size=14))+NoLegend()+xlab("PC2")
cor.test(pbmc@meta.data$PC2,pbmc@meta.data$Pseudotime)



FeatureScatter(object = pbmc, feature1 = "PC_1", feature2 = "S.Score")+ theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+  theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+scale_color_manual(values=brewer.pal(3,"Set2"))+NoLegend()+xlab("PC1")+ylab("Average expression of S phase genes")+ggtitle("")


FeatureScatter(object = pbmc, feature1 = "PC_1", feature2 = "G2M.Score")+ theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+  theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+scale_color_manual(values=brewer.pal(3,"Set2"))+xlab("PC1")+ylab("Average expression of G2/M phase genes")+ggtitle("")
cor.test(pbmc@meta.data$PC1, pbmc@meta.data$G2M.Score)
cor.test(pbmc@meta.data$PC1, pbmc@meta.data$S.Score)
?cor.test


FeatureScatter(object = pbmc, feature1 = "PC_1", feature2 = "G2M.Score", group.by="orig.ident")+ theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+  theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))




FeatureScatter(object = pbmc, feature1 = "PC_2", feature2 = "PC_3")

FeatureScatter(object = pbmc, cells=rownames(pbmc@meta.data[pbmc@active.ident == "Primary",]), feature1 = "PC2_high1", feature2 = "MITF.score1")

?FeatureScatter
FeatureScatter(object = pbmc, feature1 = "PC2", feature2 = "Pseudotime")

FeatureScatter(object = pbmc, feature1 = "PC2_high1", feature2 = "AXL.score1")

FeatureScatter(object = pbmc, feature1 = "PC2_low1", feature2 = "AXL.score1")
FeatureScatter(object = pbmc, feature1 = "PC2_low1", feature2 = "MITF.score1")

gp=FeatureScatter(object = pbmc, feature1 = "PC2_high1", feature2 = "Pseudotime")

gp+theme(plot.title = element_blank(), axis.title.x =  element_text(size=14,face="bold"), axis.title.y =  element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')+xlab("PC2 signature average expression")+ylab("Pseudotime in Monocle")

cor.test(pbmc@meta.data$PC2_high1,pbmc@meta.data$Pseudotime)

?FeatureScatter
FeatureScatter(object = pbmc, feature1 = "PC_2", feature2 = "s.g2m.score1")
FeatureScatter(object = pbmc, feature1 = "PC_1", feature2 = "s.score1")
FeatureScatter(object = pbmc, feature1 = "PC_1", feature2 = "g2m.score1")
FeatureScatter(object = pbmc, feature1 = "AXL.score1", feature2 = "MITF.score1", group.by="orig.ident")

FeatureScatter(object = pbmc, feature1 = "PC_2", feature2 = "MITF.score1", group.by="orig.ident")
FeatureScatter(object = pbmc, feature1 = "PC_2", feature2 = "AXL.score1", group.by="orig.ident")
FeatureScatter(object = pbmc, cells = rownames(pbmc@meta.data[pbmc@meta.data$RNA_snn_res.0.1=="Primary",]),feature1 = "S.Score", feature2 = "G2M.Score")
FeatureScatter(object = pbmc, cells = rownames(pbmc@meta.data[pbmc@meta.data$RNA_snn_res.0.1=="Lung",]),feature1 = "S.Score", feature2 = "G2M.Score")
FeatureScatter(object = pbmc, cells = rownames(pbmc@meta.data[pbmc@meta.data$RNA_snn_res.0.1=="Brain",]),feature1 = "S.Score", feature2 = "G2M.Score")
FeatureScatter(object = pbmc, feature1 = "G2M.Score", feature2 = "KDM5B")

####Pseudotime expression 
####Pseudotime expression 
####Pseudotime expression 
#Up-regulated
FeatureScatter(object = pbmc, feature1 = "Pseudotime", feature2 = "CCND1",pt.size=0.5,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc, feature1 = "Pseudotime", feature2 = "TMSB4X",pt.size=0.5,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')

FeatureScatter(object = pbmc, feature1 = "Pseudotime", feature2 = "PMEL",pt.size=0.5,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc, feature1 = "Pseudotime", feature2 = "APOE",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')

FeatureScatter(object = pbmc, feature1 = "Pseudotime", feature2 = "DCT",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')

#Down-regulated
FeatureScatter(object = pbmc, feature1 = "Pseudotime", feature2 = "CD74",pt.size=0.5,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc, feature1 = "Pseudotime", feature2 = "IGFBP7",pt.size=0.5,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc, feature1 = "Pseudotime", feature2 = "SLC16A3",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc, feature1 = "Pseudotime", feature2 = "MT2A",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc, feature1 = "Pseudotime", feature2 = "FN1",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')






###cell cycle active cell subsampling 
###cell cycle active cell subsampling 
###cell cycle active cell subsampling 
pbmc@active.ident
pbmc <- CellCycleScoring(pbmc, s.features = s.gene, g2m.features = g2m.gene, set.ident = TRUE)
pbmc2 <- subset(x = pbmc, idents=c("S","G2M"))
pbmc2 = FindVariableFeatures(object = pbmc2, selection.method = "vst", nfeatures = 2000)
pbmc2 <- ScaleData(object = pbmc2, features=rownames(pbmc2),vars.to.regress = c("nCount_RNA", "percent.mito"),block.size	= 5000)


pbmc2 <- RunPCA(object = pbmc2, features = VariableFeatures(object = pbmc2), verbose = FALSE)
VizDimLoadings(object = pbmc2, dims = 1:2)
DimPlot(object = pbmc2)
pbmc2 <- ProjectDim(object = pbmc2)
pbmc2 <- JackStraw(object = pbmc2, num.replicate = 100)
pbmc2 <- ScoreJackStraw(object = pbmc2, dims = 1:20)
##ElbowPlot(object = pbmc2)
pbmc2 <- FindNeighbors(object = pbmc2, dims = 1:12)
pbmc2 <- RunUMAP(pbmc2, reduction.use = "pca", dims = 1:12)
pbmc2 <- FindClusters(object = pbmc2, resolution = 0.1)
#UMAP

table(pbmc@active.ident,pbmc@meta.data$orig.ident)
pbmc@meta.data$seurat_clusterspbmc@active.ident
DimPlot(pbmc2)
DimPlot(pbmc2,group.by="orig.ident")
DimPlot(pbmc2,reduction="pca")
DimPlot(pbmc2,reduction="pca",group.by="orig.ident")
DimPlot(pbmc2,group.by="orig.ident",reduction="pca")




current.cluster.ids <- c(0,1,2)
new.cluster.ids <- c("Brain","Primary","Lung")
pbmc2@meta.data$RNA_snn_res.0.1  <- plyr::mapvalues(x = pbmc2@meta.data$RNA_snn_res.0.1 , from = current.cluster.ids, to = new.cluster.ids)
pbmc2@active.ident  <- plyr::mapvalues(x = pbmc2@active.ident  , from = current.cluster.ids, to = new.cluster.ids)
pbmc2@active.ident=factor(pbmc2@active.ident, levels=c("Primary","Lung","Brain"))


current.cluster.ids <- c(0,1)
new.cluster.ids <- c("Brain","Lung")
pbmc2@meta.data$RNA_snn_res.0.1  <- plyr::mapvalues(x = pbmc2@meta.data$RNA_snn_res.0.1 , from = current.cluster.ids, to = new.cluster.ids)
current.cluster.ids <- c("Primary","Lung","Brain1","Brain2")
new.cluster.ids <- c("WM239A","113/6-4L", "131/4-5B1", "131/4-5B2")
pbmc2@meta.data$orig.ident  <- plyr::mapvalues(x = pbmc2@meta.data$orig.ident , from = current.cluster.ids, to = new.cluster.ids)
#Confirm the object

merge.markers2 <- FindMarkers(object = pbmc2, ident.1=c("Lung"),ident.2=c("Brain"),only.pos = FALSE, group.by="RNA_snn_res.0.1")
write.csv(x = merge.markers2, row.names = TRUE, file = "DEG.markers2.csv")

my_sample_col=as.data.frame(cbind(as.data.frame(pbmc2@meta.data$orig.ident),as.data.frame(pbmc2@meta.data$RNA_snn_res.0.1)))
colnames(my_sample_col)=c("Type","Cluster")
colnames(my_sample_col)
rownames(my_sample_col)=rownames(pbmc2@meta.data)
Var1 =brewer.pal(4,"Set1")
names(Var1) = c("WM239A","113/6-4L", "131/4-5B1", "131/4-5B2")
Var2 =brewer.pal(2,"Set1")
names(Var2) = c("Lung","Brain")
ann_colors= list(Type = Var1, Cluster=Var2)
my_sample_col$Cluster=factor(my_sample_col$Cluster, levels=c("Lung","Brain"))
my_sample_col2=my_sample_col[order(my_sample_col$Cluster),]
head(my_sample_col2)
write_out= as.data.frame(as.matrix(pbmc2@assays$RNA@scale.data))


merge.markers2=merge.markers2[order(merge.markers2$avg_logFC, decreasing =TRUE),]

setwd("/home/kysbbubbu/melanomaWM237")
write.csv(x = pbmc2@meta.data, row.names = TRUE, file = "pbmc2_meta.data.csv")


data=pbmc2@meta.data[order(pbmc2@meta.data$Pseudotime2),]
data$Pseudotime2_bin=rep(c("Bin1","Bin2","Bin3","Bin4","Bin5"), c(1021,1021,1021,1021,1021))
ggplot(data, aes(Pseudotime2_bin, fill = seurat_clusters)) +  geom_bar(position = "fill") +  scale_y_continuous(labels = scales::percent)+theme_classic()+theme(plot.title = element_blank(), axis.title.x =element_blank(), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12,face="bold"), axis.text.y = element_text(size=12,hjust=0.5,face="bold"),legend.title=element_blank(),legend.text=element_text(size=12))+ylab("Proportion")



dim(write_out)
dim(my_sample_col)
write_out2=write_out[rownames(merge.markers2),order(my_sample_col$Cluster)]
dat.n=round(write_out2,digits=2)
pheatmap(dat.n,cluster_cols=FALSE, show_rownames=FALSE,show_colnames=FALSE,cluster_rows=FALSE,breaks = seq(-0.5, 0.5, by = 0.01),annotation_col =my_sample_col2,annotation_names_row=FALSE,annotation_names_col=FALSE, annotation_colors=ann_colors)
my_sample_col2$Cluster



FeatureScatter(object = pbmc2, feature1 = "PC2", feature2 = "Pseudotime")+ theme(plot.title = element_text(family = "arial", face = "bold", hjust = 0.5, size = 16))+  theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text.x = element_text(size=12), axis.text.y = element_text(size=12),legend.title=element_text(size=14),legend.text=element_text(size=12))

FeatureScatter(object = pbmc2, feature1 = "Pseudotime2", feature2 = "PC2",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_text(face="bold", size=14), axis.title.y =  element_text(face="bold", size=14),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')+xlab("Pseudotime")+ylab("PC2")
cor.test(pbmc2@meta.data$Pseudotime2,pbmc2@meta.data$PC2)


FeatureScatter(object = pbmc2, feature1 = "Pseudotime2", feature2 = "CCND1",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc2, feature1 = "Pseudotime2", feature2 = "TMSB4X",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc2, feature1 = "Pseudotime2", feature2 = "CD74",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc2, feature1 = "Pseudotime2", feature2 = "IGFBP7",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')


FeatureScatter(object = pbmc2, feature1 = "Pseudotime2", feature2 = "PMEL",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc2, feature1 = "Pseudotime2", feature2 = "APOE",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc2, feature1 = "Pseudotime2", feature2 = "DCT",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')

FeatureScatter(object = pbmc2, feature1 = "Pseudotime2", feature2 = "SLC16A3",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc2, feature1 = "Pseudotime2", feature2 = "MT2A",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')
FeatureScatter(object = pbmc2, feature1 = "Pseudotime2", feature2 = "FN1",pt.size=0.7,cols=brewer.pal(3,"Set2"))+theme(plot.title = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),legend.title=element_text(size=14),legend.text=element_text(size=12))+NoLegend()+ geom_smooth(method = 'loess', aes(colour = NA), colour = 'black')


## Monocle
## Monocle
## Monocle
gp=plot_cell_trajectory(HSMM, show_branch_points = F,cell_size = 0.2, color_by = "seurat_clusters")+ theme(axis.title.x = element_text(size=14,face="bold"), axis.title.y = element_text(size=14,face="bold"),axis.text = element_blank(),axis.ticks=element_blank(),legend.title=element_blank(),legend.text=element_text(size=12),legend.position = "none")+scale_colour_manual(name="orig.ident", values= brewer.pal(3,"Set2"))
ldat <- read.loom.matrices("merge.loom")
emat <- ldat$spliced
nmat <- ldat$unspliced
emb <- as.data.frame(cbind(gp$data$data_dim_1, gp$data$data_dim_2))
rownames(emb)=rownames(pbmc@meta.data)
cell.dist <- as.dist(1-armaCor(t(emb)))
MyData <- read.csv(file="emat_modify.csv", header=FALSE, sep=",", stringsAsFactors = FALSE)
colnames(emat) = MyData$V1
colnames(nmat) = MyData$V1
emat.new = emat[,colnames(emat) %in% rownames(emb)]
nmat.new = nmat[,colnames(nmat) %in% rownames(emb)]
emat.new2 = emat.new[rownames(emat.new) %in% rownames(pbmc@assays$RNA@data),]
nmat.new2 = nmat.new[rownames(emat.new) %in% rownames(pbmc@assays$RNA@data),]
rvel.cd_monocle <- gene.relative.velocity.estimates(emat.new2,nmat.new2, n.cores=24,verbose=TRUE,cell.dist=cell.dist)

colors <- as.list(ggplot_build(gp)$data[[2]]$colour)
names(colors) <- rownames(emb)
par(oma=c(2,2,2,2))
par(mar=c(2,2,2,2))
par(mfrow=c(1,1))
emb2 = cbind(ggplot_build(gp)$data[[2]][,"x"],ggplot_build(gp)$data[[2]][,"y"])
rownames(emb2)=names(colors)
show.velocity.on.embedding.cor(emb2,rvel.cd_monocle2,n=300,scale='sqrt',cell.colors=ac(colors,alpha=0.5),cex=0.8,arrow.scale=5,show.grid.flow=T,min.grid.cell.mass=0.5,grid.n=40,arrow.lwd=1,do.par=T,cell.border.alpha = 0.1,n.cores=24)




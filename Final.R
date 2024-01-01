##Human Ovarian Cancer, 11 mm Capture Area (FFPE) 
##https://www.10xgenomics.com/resources/datasets/human-ovarian-cancer-11-mm-capture-area-ffpe-2-standard

library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleR)
library(celldex)
library(SingleCellExperiment)


#dir.create("C:/Users/Dear User/Desktop/Single/Result", recursive = TRUE)
#setwd("D:/Drive (Y)/AAA.SingleCell workshop(FINAL)/GSM6329647_Combo/")

#countmat1=read.table("C:/Users/Dear User/Desktop/Single/step16.txt",header = TRUE)

countmat_COMBO=Read10X("C:/Users/Dear User/Desktop/Single/GSM6329644_DMSO/")
dim(countmat_COMBO)  #Col1 : Genes , Col2 : Cells   .. Genes : 36601    Cell : 9757

##############################################
        #Step1: Create Seurat Object
##############################################
srobj_COMBO  = CreateSeuratObject(countmat_COMBO, min.cells = 3,min.features = 200,project = "GSM6329647_Combo" ) #mincells : at least 3 genes per one cell , minfeature : at least in each cell 200 Geneg exist .Try to reduce and increase mincell and min features 
dim(srobj_COMBO)  # 25699  9702


##########################################################################
         #Step2: Filtering (By number of features and MT percent)
##########################################################################
srobj_COMBO[["MTpercent"]] =PercentageFeatureSet(srobj_COMBO,pattern = "^MT-" )

srobj_COMBO$MTpercent


##VlnPlot

VlnPlot(srobj_COMBO  , features = c("nFeature_RNA" ,"MTpercent" ,"nCount_RNA"))
#pdf(file = "d:/   /.pdf", title = "vlnplot", width = 8)
#png("d:/   /.png", width = 1900, height = 1500, res = 200)

#### FeatureScatter (each one should be Run after it's VlnPlot)############################################################

Plot1= FeatureScatter(srobj_COMBO, feature1 ="nCount_RNA", feature2 = "MTpercent")
Plot2=  FeatureScatter(srobj_COMBO, feature1 ="nCount_RNA", feature2 = "nFeature_RNA")
Plot1+Plot2
dev.off()

srobj_COMBO  = subset(srobj_COMBO,subset  = nFeature_RNA > 200 & nFeature_RNA<9000 & MTpercent <15)

dim(srobj_COMBO)  # Genes 25699   Cell 8977

##########################################################################
                  #Step3: Normalization
##########################################################################
srobj= NormalizeData(srobj_COMBO,normalization.method="LogNormalize", scale.factor = 10000)

srobj[["RNA"]]@data

##########################################################################
             #Step 4: Highest variable features
##########################################################################
#install.packages("ggplot2")
#require(ggplot2)
srobj = FindVariableFeatures(srobj_COMBO  , selection.method = "vst" , nfeatures = 2000)
topVariableFeatures=head(VariableFeatures(srobj), 10)
Plot1=VariableFeaturePlot(srobj)
Plot2=LabelPoints(Plot1, topVariableFeatures, size=3,xnudge = 0,ynudge =.5)
Plot2+Plot1

#########################################################################################################
           #Step 5: centering and Scaling the data (should: mean=0, variance=1)(row:genes, columns:Cells)
#########################################################################################################
srobj=ScaleData(srobj, features = rownames(srobj) )
srobj[["RNA"]]@scale.data
dim(srobj[["RNA"]]@scale.data)  # Gen 25699   cell :8977
mean(srobj[["RNA"]]@scale.data[,1])
var(srobj[["RNA"]]@scale.data[,1])



##########################################################################
               #Step 6: Dimension Reduction (PCA)
##########################################################################
dim(srobj[["RNA"]]@scale.data) # 23172   Points = 8262  ---- pc Example = 4
srobj=RunPCA(srobj,features = VariableFeatures(srobj)  )
print(srobj[["pca"]], dims= 1:10, nfeatures = 5)
VizDimLoadings(srobj , dims= 1:2 , reduction= "pca")
DimPlot(srobj , reduction= "pca")
DimHeatmap(srobj, dims = 1 , cells= 300)
DimHeatmap(srobj, dims = 1:10 , cells= 300)
ElbowPlot(srobj)

##########################################################################
               #Step 7: Determine number of Dimensions
##########################################################################
srobj = JackStraw(srobj,num.replicate  = 10 )
srobj=ScoreJackStraw(srobj, dims = 1:20)
JackStrawPlot(srobj, dims = 1:20)
ElbowPlot(srobj)
ElbowPlot(srobj,ndims = 50)

##########################################################################
#Step 8: Clustering (dims:1 to number of selected PCA that we found them good after jackStraw and Elbow plot)(resolution: choosing the best cluster visualization until 3000 cell)
##########################################################################


#(Lenghth(Levels): How many clusters?)(Idents:each cells in which cluster)
srobj=FindNeighbors(srobj,dims = 1:20)
srobj=FindClusters(srobj,resolution = 1 ) #(resolution defult:(0.4:1,2 up to 3K cells))
Idents(srobj)
levels(Idents(srobj))
length(levels(Idents(srobj)))
head(Idents(srobj),5)

##########################################################################
#Step 9: Dimension Reduction (UMAP - Non-Linear Dimention Reduction)
##########################################################################
srobj=RunUMAP(srobj,dims = 1:10) #100
DimPlot(srobj,label = T,reduction = "umap")

srobj=RunTSNE(srobj,dims = 1:50) #100
DimPlot(srobj,label = T,reduction = "tsne")

##########################################################################
#Step 6: FindAllMarkers
##########################################################################

markers=FindAllMarkers(srobj,min.pct=0.25, logfc.threshold= 1 )#  0.25  # Delayed too much
dim(markers)
View(markers)
nrow(markers)

summary(markers$avg_log2FC)
min(markers$avg_log2FC)
min(markers$avg_log2FC[which(markers$avg_log2FC>0)])
markers[which(markers$avg_log2FC>=1),]
dim(markers[which(markers$avg_log2FC ==1),])
View(markers[which(markers$avg_log2FC>=1),])          #290
###################################################################################

      #Manual Markers 

###################################################################################

markers=markers[which(markers$avg_log2FC>=0),]
cl=markers[which(markers$cluster==1),]

DoHeatmap(srobj, features= c("HIST1H4C", "NEAT1", "SOD2", "MYBL2", "IFITM2"))
DoHeatmap(srobj, features= cl$gene)

VlnPlot(srobj, features= c("HIST1H4C", "NEAT1", "SOD2", "MYBL2", "IFITM2"))
RidgePlot(srobj, features= c("HIST1H4C", "NEAT1", "SOD2", "MYBL2", "IFITM2"))
RidgePlot(srobj, features= cl$gene)

FeaturePlot(srobj, features= cl$gene)
FeaturePlot(srobj, features= c("HIST1H4C", "NEAT1", "SOD2", "MYBL2", "IFITM2"))

head(srobj@meta.data)

###################################################################################

#Manual Celltype 

###################################################################################
cellclusterIDs=0:5
celltypes=c("Actrocytes","Tcell","Epithelialcell","Neuron","Ips","B-Cell")
names(x=celltypes)=levels(x=srobj)
srobj2=RenameIdents(object = srobj,celltypes)
DimPlot(srobj2,label = TRUE)
x=GetAssayData(subset(srobj2,idents = "Neuron"),"counts")
write.csv(x,file = "e:/Neuron.csv")

############################ ENDDDDDDDDDDDDDDDDDDDDDDDDDDDDD





Ref=celldex::HumanPrimaryCellAtlasData(ensembl = F)
predLabs=SingleR(test=as.SingleCellExperiment(srobj),ref= Ref,labels =Ref$label.main  )
predLabs$scores
length(predLabs$labels)
View(as.data.frame(table(predLabs$labels)))

#x=as.data.frame(table(predLabs$labels))
#write.csv(x,file = "d:/x.csv")

View(predLabs$scores)

srobj$Celltypes=predLabs$labels
DimPlot(srobj,reduction = "umap",label = T)
DimPlot(srobj,reduction = "umap",label = T,group.by = "Celltypes")



### test witch gen belongd to wich celltype ? 


which(predLabs$labels == "Endothelial_cells") # this Cells 
srobj$seurat_clusters

table(srobj$seurat_clusters[which(predLabs$labels == "Neuron")])


table(srobj$seurat_clusters[which(predLabs$labels == "BM & Prog.")]   )          # 5 
table(srobj$seurat_clusters[which(predLabs$labels == "Smooth_muscle_cells")] )#   5
table(srobj$seurat_clusters[which(predLabs$labels == "Fibroblasts")]       ) # 5
table(srobj$seurat_clusters[which(predLabs$labels == "Keratinocytes")]      )  #  5
table(srobj$seurat_clusters[which(predLabs$labels == "Endothelial_cells")] ) #  5 ,  6
table(srobj$seurat_clusters[which(predLabs$labels == "Neuroepithelial_cell")] ) # 5 , 4 , 8
table(srobj$seurat_clusters[which(predLabs$labels == "Neurons")]  )            # 0 , 4 , 5 , 6 , 7 , 8 
table(srobj$seurat_clusters[which(predLabs$labels == "MSC")]     )             # 5 , 4 , 1 
table(srobj$seurat_clusters[which(predLabs$labels == "Embryonic_stem_cells")] )# 124 - 
table(srobj$seurat_clusters[which(predLabs$labels == "iPS_cells")])
table(srobj$seurat_clusters[which(predLabs$labels == "Epithelial_cells")])

ClustersID = 0:8
Celltypes = c("Epithelial_cells","iPS_cells","iPS_cells","Epithelial_cells","Embryonic_stem_cells","MSC","Endothelial_cells","Epithelial_cells","Neurons")

names(Celltypes)=levels(srobj)
srobj=RenameIdents(obj=srobj, Celltypes)

DimPlot(srobj,reduction = "umap",label = T)
DimPlot(srobj,reduction = "umap",label = T,group.by = "Celltypes")

GetAssayData(subset(srobj,idents="Epithelial_cells"),"counts")

x=GetAssayData(subset(srobj,idents = "Epithelial_cells"),"counts")
write.csv(x,file = "e:/Epithelial_cells202422.csv")




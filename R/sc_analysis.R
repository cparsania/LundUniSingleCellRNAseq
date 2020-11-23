library(Seurat)
??Seurat

## load data 
d14_ventral_cdna.data <- Seurat::Read10X(data.dir = "../data/")


## initialize seurat object with at least 200 genes in each cell

d14_ventral_cdna.seurat <- CreateSeuratObject(counts = d14_ventral_cdna.data, 
                           project = "Day14_Ventral_cDNA", 
                           min.features = 200)

## identify  mito genes (start with "MT-")
d14_ventral_cdna.seurat[["percent.mt"]] <- PercentageFeatureSet(d14_ventral_cdna.seurat, pattern = "^MT-")
head(d14_ventral_cdna.seurat@meta.data, 5)

## plot feature counts
VlnPlot(d14_ventral_cdna.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

## meta data average 
meta_data_avg <- d14_ventral_cdna.seurat@meta.data %>% 
        tibble::as_tibble() %>% 
        dplyr::summarise_if(is.numeric , mean)

## average feature 
avg_nFeature_RNA <- meta_data_avg %>% 
        dplyr::pull(nFeature_RNA)

## remove cells expressing more than 10% mitochondrial genes and 
## cells expressing genes more than 2.5 times of average number of genes expressed in dataset

d14_ventral_cdna.seurat <- subset(d14_ventral_cdna.seurat, 
                                  subset = nFeature_RNA < (nFeature_RNA * 2.5) & percent.mt < 10)

## normalization 
d14_ventral_cdna.seurat <- NormalizeData(d14_ventral_cdna.seurat, 
                                         normalization.method = "LogNormalize", 
                                         scale.factor = 10000)

## Find highly variable features

d14_ventral_cdna.seurat <- FindVariableFeatures(d14_ventral_cdna.seurat, 
                                                selection.method = "vst", 
                                                nfeatures = 2000)

top10 <- head(VariableFeatures(d14_ventral_cdna.seurat), 10)

plot1 <- VariableFeaturePlot(d14_ventral_cdna.seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## scale variable features
d14_ventral_cdna.seurat <- ScaleData(d14_ventral_cdna.seurat, 
                                     features = VariableFeatures(d14_ventral_cdna.seurat))


## dimension reduction 
## run PCA using variable features 
d14_ventral_cdna.seurat <- RunPCA(d14_ventral_cdna.seurat, 
                                  features = VariableFeatures(d14_ventral_cdna.seurat))

DimPlot(d14_ventral_cdna.seurat, reduction = "pca")



## find optimum PCs.
DimHeatmap(d14_ventral_cdna.seurat, dims = 1:20, cells = 500, balanced = TRUE)
elb_plt <- ElbowPlot(d14_ventral_cdna.seurat,ndims = 50) 
elb_plt + geom_vline(xintercept = 14)


## test number of cluster with different sets of PCs
num_of_pcs <- seq(10,20, by =2 )

umap_plots <- purrr::map(num_of_pcs , function(x){
        d14_ventral_cdna.seurat <- FindNeighbors(d14_ventral_cdna.seurat, dims = 1:x)
        d14_ventral_cdna.seurat <- FindClusters(d14_ventral_cdna.seurat, resolution = 0.6)
        d14_ventral_cdna.seurat <- RunUMAP(d14_ventral_cdna.seurat, dims = 1:x)
        dim_plot <- DimPlot(d14_ventral_cdna.seurat, reduction = "umap",label = T)
        dim_plot + ggtitle(glue::glue("No. of PCs : {x}"))
})

ggpubr::ggarrange(umap_plots[[1]] + theme(legend.position = "none"),
                  umap_plots[[2]] + theme(legend.position = "none"),
                  umap_plots[[3]]+ theme(legend.position = "none"),
                  umap_plots[[4]]+ theme(legend.position = "none"),
                  umap_plots[[5]]+ theme(legend.position = "none"),
                  umap_plots[[6]]+ theme(legend.position = "none"))


## Elbow plot and UMAP plot with different number of PCs agree to have optimum number of PCs to 16. 

num_of_pc = 16

## cell clustering based on previously identified PCs. 

d14_ventral_cdna.seurat <- FindNeighbors(d14_ventral_cdna.seurat, 
                                         dims = 1:num_of_pc,
                                         reduction = "pca")

d14_ventral_cdna.seurat <- FindClusters(d14_ventral_cdna.seurat, resolution =0.6)

## how to decide resolution in the function FindCluster? 

head(Idents(d14_ventral_cdna.seurat), 5)


## non linear dimensional reduction UMAP

d14_ventral_cdna.seurat <- RunUMAP(d14_ventral_cdna.seurat , dims = 1:num_of_pc)
dim_plot <- DimPlot(d14_ventral_cdna.seurat, reduction = "umap",label = T)
dim_plot + ggtitle(glue::glue("No. of PCs : {num_of_pc}"))



## Find top markers of each cluster 

top.markers <- FindAllMarkers(d14_ventral_cdna.seurat, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

## no. of  markers by cluster  
top.markers %>% rownames_to_column() %>% as_tibble() %>% group_by(cluster) %>% summarise(n = n())

##plot top marker of each cluster 

top.markers.1 <- top.markers %>% 
        tibble::rownames_to_column() %>% 
        tibble::as_tibble() %>% 
        dplyr::group_by(cluster) %>% 
        dplyr::arrange(-avg_logFC) %>% 
        slice(1)

fp <- FeaturePlot(d14_ventral_cdna.seurat , features = top.markers.1$gene )


## cell cycle regression 

s.genes <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8")
g2m.genes <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")

d14_ventral_cdna.seurat <- CellCycleScoring(d14_ventral_cdna.seurat, 
                                            s.features = s.genes, 
                                            g2m.features = g2m.genes, 
                                            set.ident = FALSE)


## scale data  without cell cycle regression 
d14_ventral_cdna.seurat <- ScaleData(d14_ventral_cdna.seurat)

## PCA with cell cycle features 
d14_ventral_cdna.seurat <- RunPCA(d14_ventral_cdna.seurat, features = c(s.genes, g2m.genes))
DimPlot(d14_ventral_cdna.seurat , reduction = "pca")

## PCA with all features 
d14_ventral_cdna.seurat <- RunPCA(d14_ventral_cdna.seurat)
DimPlot(d14_ventral_cdna.seurat , reduction = "pca")


## scale data  with cell cycle regression 
d14_ventral_cdna.seurat <- ScaleData(d14_ventral_cdna.seurat, vars.to.regress = c("S.Score", "G2M.Score"))

## PCA with cell cycle features 
d14_ventral_cdna.seurat <- RunPCA(d14_ventral_cdna.seurat, features = c(s.genes, g2m.genes))
DimPlot(d14_ventral_cdna.seurat , reduction = "pca")

## PCA with cell cycle features 
d14_ventral_cdna.seurat <- RunPCA(d14_ventral_cdna.seurat)
DimPlot(d14_ventral_cdna.seurat , reduction = "pca")


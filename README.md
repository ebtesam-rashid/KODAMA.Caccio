# KODAMA.tkcaccia
Enhanced dimensionality reduction for High Throughput omics data

## Overview 

KODAMA is an unsupervised and semi-supervised learning algorithm that performs feature extraction from noisy and high-dimensional data. It facilitates identification of patterns representing underlying groups on all samples in a data set. 

This is a version in developing of KODAMA with landmarks to adapt the algorithm to the analysis of data set with more than 10,000 entries. The wrapper for the C++ implementation of Barnes-Hut t-Distributed Stochastic Neighbor Embedding has been integrated to convert the KODAMA's dissimilarity matrix in a low dimensional space.

KODAMA was built on accuracy maximization algorithms described in detail in the following publication:
[Zinga, M. M., Abdel-Shafy, E., Melak, T., Vignoli, A., Piazza, S., Zerbini, L. F., ... & Cacciatore, S. (2022). KODAMA exploratory analysis in metabolic phenotyping. Frontiers in Molecular Biosciences, 9.] (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9887019/)
[Cacciatore, S., Tenori, L., Luchinat, C., Bennett, P. R., & MacIntyre, D. A. (2017). KODAMA: an R package for knowledge discovery and data mining. Bioinformatics, 33(4), 621-623.] (https://academic.oup.com/bioinformatics/article/33/4/621/2667156?login=false)
[Cacciatore, S., Luchinat, C., & Tenori, L. (2014). Knowledge discovery by accuracy maximization. Proceedings of the National Academy of Sciences, 111(14), 5117-5122.] (https://www.pnas.org/doi/abs/10.1073/pnas.1220873111)


 

## Installation

```



## Metabolomic data

The data belong to a cohort of 22 healthy donors (11 male and 11 female) where each provided about 40 urine samples over the time course of approximately 2 months, for a total of 873 samples. Each sample was analysed by Nuclear Magnetic Resonance Spectroscopy. Each spectrum was divided in 450 spectral bins.

```
data(MetRef)
u=MetRef$data;
u=u[,-which(colSums(u)==0)]
u=normalization(u)$newXtrain
u=scaling(u)$newXtrain
class=as.numeric(as.factor(MetRef$gender))
plot(cc,pch=21,bg=class,xlab="First Component",ylab="Second Component")

class=as.numeric(as.factor(MetRef$donor))
kk=KODAMA(u,landmarks = 100,perplexity = 10,f.par = 50)
tt=Rtsne(u,perplexity = 10)

par(mfrow=c(1,2))
plot(kk$scores,pch=21,bg=rainbow(22)[class], main="KODAMA-tSNE")
plot(tt$Y,pch=21,bg=rainbow(22)[class], main="tSNE, xlab="First Component", ylab="Second Component")
```
![This is an image](https://github.com/tkcaccia/Documents/blob/main/MetRef.png)


## Single-cell data

The data set from Tasic et al. encompasses 23,822 cells from adult mouse cortex, split by the authors into 133 clusters with strong hierarchical organisation. A standard preprocessing pipeline consisting of sequencing depth normalisation, feature selection, log-transformation, and reducing the dimensionality to 50 PCs was applied as described by Kobak & Berens in [The art of using t-SNE for single-cell transcriptomics](https://www.nature.com/articles/s41467-019-13056-x).

Download the data from [here](http://celltypes.brain-map.org/rnaseq) and unpack. Direct links: [VISp](http://celltypes.brain-map.org/api/v2/well_known_file_download/694413985), [ALM](http://celltypes.brain-map.org/api/v2/well_known_file_download/694413179).
To get the information about cluster colors and labels (sample_heatmap_plot_data.csv), open the interactive [data browser](http://celltypes.brain-map.org/rnaseq/mouse/v1-alm), go to "Sample Heatmaps", click "Build Plot!" and then "Download data as CSV".

```
ta=read.csv("tasic-sample_heatmap_plot_data.txt")
rownames(ta)=ta[,1]
VIS=read.csv("mouse_VISp_gene_expression_matrices_2018-06-14/mouse_VISp_2018-06-14_exon-matrix.csv")
ALM=read.csv("mouse_ALM_gene_expression_matrices_2018-06-14/mouse_ALM_2018-06-14_exon-matrix.csv")

data=t(cbind(ALM,VIS))
colnames(data)=as.character(data[1,])

data=data[-1,]

ii=intersect(rownames(data),rownames(ta))
data=data[ii,]

data=data[,colSums(data)!=0]

near.zero.counts=colMeans(data<32)
temp=data
temp[temp<=32]=NA
temp=log2(temp)
m=colMeans(temp,na.rm = TRUE)

y=exp(-1.5*(m-6.56))+0.02


data=data[,which(near.zero.counts>y)]

su=rowSums(data)
data=((data/su)*10^6)*median(su)
data=log2(data+1)

pca=prcomp(data)$x[,1:50]

kk=KODAMA(pca,landmarks = 1000)
plot(kk$scores,pch=21,bg=ta[,"cluster_color"])

```
![This is an image](https://github.com/tkcaccia/Documents/blob/main/Tasic.png)



# GeoMX data
The GeoMx Digital Spatial Profiler (DSP) is a platform for capturing spatially resolved high-plex gene (or protein) expression data from tissue [Merritt et al., 2020](https://pubmed.ncbi.nlm.nih.gov/32393914/). In particular, formalin-fixed paraffin-embedded (FFPE) or fresh-frozen (FF) tissue sections are stained with barcoded in-situ hybridization probes that bind to endogenous mRNA transcripts. 
GeoMx kidney dataset has been created with the human whole transcriptome atlas (WTA) assay. The dataset includes 4 diabetic kidney disease (DKD) and 3 healthy kidney tissue samples. Regions of interest (ROI) were spatially profiled to focus on two different kidney structures: tubules or glomeruli. One glomerular ROI contains the entirety of a single glomerulus. Each tubular ROI contains multiple tubules that were segmented into distal (PanCK+) and proximal (PanCK-) tubule areas of illumination (AOI). The preprocessing workflow is described [here](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html).
An imputing procedure was added to the original [R script](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.R).


```
data=t(log2(assayDataElement(target_demoData , elt = "q_norm")))
data[is.infinite(data)]=NA
data=impute.knn(data)$data

#data=prcomp(data)$x[,1:100]


# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42
# run UMAP
umap_out <-
  umap(data,  
       config = custom_umap)
pData(target_demoData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
plot1=ggplot(pData(target_demoData),
       aes(x = UMAP1, y = UMAP2, color = region, shape = class)) +
  geom_point(size = 3) +
  theme_bw()

# run tSNE
set.seed(42) # set the seed for tSNE as well
tsne_out <-
  Rtsne(data,
        perplexity = ncol(target_demoData)*.15)
pData(target_demoData)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
plot2=ggplot(pData(target_demoData),
       aes(x = tSNE1, y = tSNE2, color = segment, shape = class)) +
  geom_point(size = 3) +
  theme_bw()


kk=KODAMA(data)
pData(target_demoData)[, c("KODAMA1", "KODAMA2")] <- kk$scores
plot3=ggplot(pData(target_demoData),
       aes(x = KODAMA1, y = KODAMA2, color = segment, shape = class)) +
  geom_point(size = 3) +
  theme_bw()



require("gridExtra")
grid.arrange(plot1, plot2, plot3, ncol=3)
```

![This is an image](https://github.com/tkcaccia/Documents/blob/main/GeoMX.png)

To reduce the computational time, the first 50 principal components can be used as input of KODAMA

```


data.pca=prcomp(data)$x[,1:50]

kkpca=KODAMA(data.pca)
pData(target_demoData)[, c("KODAMA1", "KODAMA2")] <- kkpca$scores
ggplot(pData(target_demoData),
             aes(x = KODAMA1, y = KODAMA2, color = segment, shape = class)) +
  geom_point(size = 3) +
  theme_bw()
```
![This is an image](https://github.com/tkcaccia/Documents/blob/main/GeoMX2.png)





# KODAMA for Spatial Transcriptomics

To show the how KODAMA could deal with noisy datasets compared to other dimensionality reduction widely used methods such as Uniform Manifold Approximation and Projection (UMAP) and t-Distributed Stochastic Neighbour Embedding (t-SNE), the following simulated data will be used.

## Example 1: Simulated data set 

KODAMA, tSNE, UMAP are applied to simulated data set of two dimention with different degrees of noisy (from 0 to 20). The following script to compare the effect of different dimentionallity reduction algorithms on a simulated data set of 8 noisy dimensions.

### Tutorial

1. The data is simulated with vertix function from KODAMA package with 2 dimensions and 8 noisy dimensions.
2. The generated data is scaled

```
ma=scale(ma)
```

3. Apply MDS, tSNE, KODAMA
```
res_MDS=cmdscale(dist(ma))
res_tSNE=Rtsne(ma)$Y
res_UMAP = umap(ma)$layout
```
4. Apply KODAMA
```
kk=KODAMA.matrix(ma)
res_KODAMA_MDS=KODAMA.visualization(kk,method = "MDS")
res_KODAMA_tSNE=KODAMA.visualization(kk,method = "t-SNE")
res_KODAMA_UMAP=KODAMA.visualization(kk,method = "UMAP")
```

5. Plot the results

```
par(mfrow = c(2,3))
labels <- rep(c("#FF0000","#0000FF","#008000","#FFFF00"),each= 50)
plot(res_MDS,pch=21,bg=labels,main="MDS")
plot(res_tSNE,pch=21,bg=labels,main="tSNE")
plot(res_UMAP,pch=21,bg=labels,main="UMAP")
plot(res_KODAMA_MDS,pch=21,bg=labels,main="KODAMA_MDS",ylim=range(res_KODAMA_MDS[,1]))
plot(res_KODAMA_tSNE,pch=21,bg=labels,main="KODAMA_tSNE")
plot(res_KODAMA_UMAP,pch=21,bg=labels,main="KODAMA_UMAP")
```
<p>
  <p align="center">
    <img src="https://github.com/ebtesam-rashid/KODAMA.Caccio/blob/main/Figures/final%20one%20simulated.png" alt="hello-light" height="500" width="700" />
  </p>
</p>

## Example 2: GEOMx data
The GeoMx Digital Spatial Profiler (DSP) is a platform for capturing spatially resolved high-plex gene (or protein) expression data from tissue [Merritt et al., 2020](https://pubmed.ncbi.nlm.nih.gov/32393914/). In particular, formalin-fixed paraffin-embedded (FFPE) or fresh-frozen (FF) tissue sections are stained with barcoded in-situ hybridization probes that bind to endogenous mRNA transcripts. 
GeoMx kidney dataset has been created with the human whole transcriptome atlas (WTA) assay. The dataset includes 4 diabetic kidney disease (DKD) and 3 healthy kidney tissue samples. Regions of interest (ROI) were spatially profiled to focus on two different kidney structures: tubules or glomeruli. One glomerular ROI contains the entirety of a single glomerulus. Each tubular ROI contains multiple tubules that were segmented into distal (PanCK+) and proximal (PanCK-) tubule areas of illumination (AOI). The preprocessing workflow is described [here](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html).
An imputing procedure was added to the original [R script](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.R).

### Tutorial

#### Install required packages
```
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("impute")
install.packages("KODAMA")
library(impute)
library(KODAMA)
```
#### Ulpoad data
```
data=t(log2(assayDataElement(target_demoData , elt = "q_norm")))
data[is.infinite(data)]=NA
data=impute.knn(data)$data
```
#### Run PCA
```
data=prcomp(data)$x[,1:100]
```
#### Run MDA
```
MDS_out=cmdscale(dist(data))
pData(target_demoData)[, c("MDS1", "MDS2")] <- MDS_out[, c(1,2)]
```
#### run tSNE
```
set.seed(42) # set the seed for tSNE as well
tsne_out <- Rtsne(data, perplexity = ncol(target_demoData)*.15)
pData(target_demoData)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
```
#### run UMAP
```
umap_out <- umap(data, config = custom_umap)
pData(target_demoData)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
```
# run KODAMA
```
kk=KODAMA.matrix(data)
res= KODAMA.visualization(kk)
res1= KODAMA.visualization(kk,method = "MDS")
res2= KODAMA.visualization(kk,method = "t-SNE")
res3= KODAMA.visualization(kk,method = "UMAP")
pData(target_demoData)[, c("KODAMA1.MDS", "KODAMA2.MDS")] <- res1
pData(target_demoData)[, c("KODAMA1.tSNE", "KODAMA2.tSNE")] <- res2
pData(target_demoData)[, c("KODAMA1.UMAP", "KODAMA2.UMAP")] <- res3
```
#### MDS vs KODAMA.MDS
```
require("gridExtra")
plot1=ggplot(pData(target_demoData), aes(x = MDS1, y = MDS2,color = segment, shape = class)) + geom_point(size = 3) + theme_bw()
plot2=ggplot(pData(target_demoData), aes(x = KODAMA1.MDS, y = KODAMA2.MDS, color = segment, shape = class)) + geom_point(size = 3) + theme_bw()
grid.arrange(plot1, plot4, ncol=2)
```
<p>
  <p align="center">
    <img src="https://github.com/ebtesam-rashid/KODAMA.Caccio/blob/main/Figures/MDS%20geomx.png" alt="hello-light" height="500" width="700" />
  </p>
</p>

#### tSNA vs KODAMA.tSNE
```
plot3=ggplot(pData(target_demoData), aes(x = tSNE1, y = tSNE2, color = segment, shape = class)) + geom_point(size = 3) + theme_bw()
plot4=ggplot(pData(target_demoData), aes(x = KODAMA1.tSNE, y = KODAMA2.tSNE, color = segment, shape = class)) + geom_point(size = 3) + theme_bw()
grid.arrange(plot2, plot5, ncol=2)
```
<p>
  <p align="center">
    <img src="https://github.com/ebtesam-rashid/KODAMA.Caccio/blob/main/Figures/tsne%20geomx.png" alt="hello-light" height="500" width="700" />
  </p>
</p>

#### UMAP vs KODAM.UMAP
```
plot5=ggplot(pData(target_demoData), aes(x = UMAP1, y = UMAP2, color = segment, shape = class)) + geom_point(size = 3) + theme_bw()
plot6=ggplot(pData(target_demoData), aes(x = KODAMA1.UMAP, y = KODAMA2.UMAP, color = segment, shape = class)) + geom_point(size = 3) + theme_bw()
grid.arrange(plot3, plot6, ncol=2)
```
<p>
  <p align="center">
    <img src="https://github.com/ebtesam-rashid/KODAMA.Caccio/blob/main/Figures/umap%20geomx.png" alt="hello-light" height="500" width="700" />
  </p>
</p>


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





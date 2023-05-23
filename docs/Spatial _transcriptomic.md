## KODAMA for Spatial Transcriptomics

To show the how KODAMA could deal with noisy datasets compared to other dimensionality reduction widely used methods such as Uniform Manifold Approximation and Projection (UMAP) and t-Distributed Stochastic Neighbour Embedding (t-SNE), the following simulated data will be used.

```
library("KODAMA")
library("cluster")
library("vertex")

dimensions=2
size_cluster=50
cluster_number=2^dimensions
noisy_dimension <- 9
v=matrix(rep(vertex(c(0,10),dims=dimensions),each=size_cluster),ncol=dimensions)
ma=v+rnorm(length(v),sd = 0.2)
if(noisy_dimension>0){
  ma=cbind(ma,matrix(rnorm(cluster_number*size_cluster*noisy_dimension),ncol=noisy_dimension))}
ma=scale(ma)

res_MDS=cmdscale(dist(ma))
res_tSNE=Rtsne(ma,perplexity = 20)$Y
custom.settings = umap.defaults
custom.settings$n_neighbors=20
res_UMAP = umap(ma, config = custom.settings)$layout

kk=KODAMA.matrix(ma)
res_KODAMA_MDS=KODAMA.visualization(kk,method = "MDS")
res_KODAMA_tSNE=KODAMA.visualization(kk,method = "t-SNE")
res_KODAMA_UMAP=KODAMA.visualization(kk,method = "UMAP")

par(mfrow = c(2,3), oma = c(2,2,0,0) + 0.1, mar = c(1,1,1,1) + 1)

labels <- rep(c("#FF0000","#0000FF","#008000","#FFFF00"),each= 50)
plot(res_MDS,pch=21,bg=labels,main="MDS")
plot(res_tSNE,pch=21,bg=labels,main="tSNE")
plot(res_UMAP,pch=21,bg=labels,main="UMAP")
plot(res_KODAMA_MDS,pch=21,bg=labels,main="KODAMA_MDS",ylim=range(res_KODAMA_MDS[,1]))
plot(res_KODAMA_tSNE,pch=21,bg=labels,main="KODAMA_tSNE")
plot(res_KODAMA_UMAP,pch=21,bg=labels,main="KODAMA_UMAP")
title(xlab = "Fisrt dimention", ylab = "Second dimention", outer = TRUE, line = 0.01)

```
![This is an image](https://github.com/ebtesam-rashid/KODAMA.Caccio/blob/main/Figures/simulated%20data%201.png)


## Examples

## KODAMA for Spatial Transcriptomics

To show the how KODAMA could deal with noisy datasets compared to other dimensionality reduction widely used methods such as Uniform Manifold Approximation and Projection (UMAP) and t-Distributed Stochastic Neighbour Embedding (t-SNE), the following simulated data will be used.

```
library("KODAMA")
library("cluster")
library("vertex")

dimensions=2
size_cluster=50
cluster_number=2^dimensions
noisy_dimension <- 5
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

custom.settings = Rtsne.defaults
custom.settings$perplexity=20
res_KODAMA_tSNE=KODAMA.visualization(kk,method = "t-SNE",config = custom.settings)

custom.settings = umap.defaults
custom.settings$n_neighbors=20
res_KODAMA_UMAP=KODAMA.visualization(kk,method = "UMAP",config = custom.settings)

par(mar = c(3,3,3,3))
par(mfrow = c(2,3))

plot1 <-plot(res_MDS,pch=21,bg=rainbow(4),main="MDS")
plot2 <-plot(res_tSNE,pch=21,bg=rainbow(4),main="tSNE")
plot3 <-plot(res_UMAP,pch=21,bg=rainbow(4),main="UMAP")

plot4 <-plot(res_KODAMA_MDS,pch=21,bg=rainbow(4),main="KODAMA_MDS",ylim=range(res_KODAMA_MDS[,1]))
plot5 <-plot(res_KODAMA_tSNE,pch=21,bg=rainbow(4),main="KODAMA_tSNE")
plot6 <-plot(res_KODAMA_UMAP,pch=21,bg=rainbow(4),main="KODAMA_UMAP")

```
![This is an image](https://github.com/ebtesam-rashid/KODAMA.Caccio/blob/main/Figures/simulated.data%205.png)


## Examples
1. Simulated data
2. GEOMx data

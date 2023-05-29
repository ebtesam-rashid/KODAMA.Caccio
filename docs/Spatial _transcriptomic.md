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
    <img src="https://github.com/ebtesam-rashid/KODAMA.Caccio/blob/main/Figures/one%20simulated.png" alt="hello-light" height="700" width="800" />
  </p>
</p>



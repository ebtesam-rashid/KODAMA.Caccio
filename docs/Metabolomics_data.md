## Metabolomic data

The data belong to a cohort of 22 healthy donors (11 male and 11 female) where each provided about 40 urine samples over the time course of approximately 2 months, for a total of 873 samples. Each sample was analysed by Nuclear Magnetic Resonance Spectroscopy. Each spectrum was divided in 450 spectral bins.

```

library(KODAMA)
data(MetRef)
u=MetRef$data
u=u[,-which(colSums(u)==0)]
u=normalization(u)$newXtrain
u=scaling(u)$newXtrain
class=as.numeric(as.factor(MetRef$donor))
res_MDS=cmdscale(dist(u))
res_tSNE=Rtsne(u,perplexity = 20)$Y
custom.settings = umap.defaults
custom.settings$n_neighbors=20
res_UMAP = umap(u, config = custom.settings)$layout

kk=KODAMA.matrix(u,f.par = 50)

res_KODAMA_MDS=KODAMA.visualization(kk,method = "MDS")
res_KODAMA_tSNE=KODAMA.visualization(kk,method = "t-SNE")
res_KODAMA_UMAP=KODAMA.visualization(kk,method = "UMAP",config = custom.settings)

par(mfrow = c(2,3), oma = c(2,2,0,0) + 0.1, mar = c(1,1,1,1) + 1)

plot(res_MDS,pch=21,bg=rainbow(22)[class],main="MDS")
plot(res_tSNE,pch=21,bg=rainbow(22)[class],main="tSNE")
plot(res_UMAP,pch=21,bg=rainbow(22)[class],main="UMAP")
plot(res_KODAMA_MDS,pch=21,bg=rainbow(22)[class],main="KODAMA_MDS",ylim=range(res_KODAMA_MDS[,1]))
plot(res_KODAMA_tSNE,pch=21,bg=rainbow(22)[class],main="KODAMA_tSNE")
plot(res_KODAMA_UMAP,pch=21,bg=rainbow(22)[class],main="KODAMA_UMAP")
title(xlab = "Fisrt dimention", ylab = "Second dimention", outer = TRUE, line = 0.01)
```
<p>
  <p align="center">
    <img src="https://github.com/ebtesam-rashid/KODAMA.Caccio/blob/main/Figures/metabolomics.png" alt="hello-light" height="700" width="800" />
  </p>
</p>




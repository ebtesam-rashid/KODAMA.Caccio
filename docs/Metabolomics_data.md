## Metabolomic data

The data belong to a cohort of 22 healthy donors (11 male and 11 female) where each provided about 40 urine samples over the time course of approximately 2 months, for a total of 873 samples. Each sample was analysed by Nuclear Magnetic Resonance Spectroscopy. Each spectrum was divided in 450 spectral bins.

```

library(KODAMA)
data(MetRef)
u=MetRef$data
u=u[,-which(colSums(u)==0)]
u=normalization(u)$newXtrain
u=scaling(u)$newXtrain
class=as.numeric(as.factor(MetRef$gender))
class2=as.numeric(as.factor(MetRef$donor))

tt=Rtsne(u,perplexity = 10)

kk=KODAMA.matrix(u)
custom.settings = Rtsne.defaults
custom.settings$perplexity=20
res_KODAMA_tSNE=KODAMA.visualization(kk,method = "t-SNE",config = custom.settings)

par(mfrow=c(1,2))
plot(tt$Y,pch=21,bg=rainbow(22)[class2], main="tSNE", ylim= range(tt$Y[,1]), xlim = range(tt$Y[,2]))
plot(res_KODAMA_tSNE,pch=21,bg=rainbow(22)[class2], main="KODAMA-tSNE", ylim= range(tt$Y[,1]), xlim = range(tt$Y[,2]))

```

![This is an image](https://github.com/ebtesam-rashid/KODAMA.Caccio/blob/main/Figures/KODAMA-tsne.png))

# all KODAMA

```

library(KODAMA)
data(MetRef)
u=MetRef$data
u=u[,-which(colSums(u)==0)]
u=normalization(u)$newXtrain
u=scaling(u)$newXtrain
class=as.numeric(as.factor(MetRef$gender))
class2=as.numeric(as.factor(MetRef$donor))

res_MDS=cmdscale(dist(u))
res_tSNE=Rtsne(u,perplexity = 20)$Y
custom.settings = umap.defaults
custom.settings$n_neighbors=20
res_UMAP = umap(u, config = custom.settings)$layout

kk=KODAMA.matrix(u)

custom.settings = Rtsne.defaults
custom.settings$perplexity=20
res_KODAMA_tSNE=KODAMA.visualization(kk1,method = "t-SNE",config = custom.settings)
    
custom.settings = umap.defaults
custom.settings$n_neighbors=20
res_KODAMA_UMAP=KODAMA.visualization(kk1,method = "UMAP",config = custom.settings)

par(mfrow=c(2,3))

par(mfrow=c(2,3))

plot1 <-plot(res_MDS,pch=21,bg=rainbow(22)[class2],main="MDS", xlab= "Fisrt dimention", ylab= "Second dimention")
plot2 <-plot(res_tSNE,pch=21,bg=rainbow(22)[class2],main="tSNE")
plot3 <-plot(res_UMAP,pch=21,bg=rainbow(22)[class2],main="UMAP")

plot4 <-plot(res_KODAMA_MDS,pch=21,bg=rainbow(22)[class2],main="KODAMA_MDS",ylim=range(res_KODAMA_MDS[,1]))
plot5 <-plot(res_KODAMA_tSNE,pch=21,bg=rainbow(22)[class2],main="KODAMA_tSNE")
plot6 <-plot(res_KODAMA_UMAP,pch=21,bg=rainbow(22)[class2],main="KODAMA_UMAP")

```

![This is an image](https://github.com/ebtesam-rashid/KODAMA.Caccio/blob/main/Figures/all%20METAB%20KODAMA.png)



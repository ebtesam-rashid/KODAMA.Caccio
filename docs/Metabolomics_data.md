library(KODAMA)
data(MetRef)
u=MetRef$data
u=u[,-which(colSums(u)==0)]
u=normalization(u)$newXtrain
u=scaling(u)$newXtrain
class=as.numeric(as.factor(MetRef$gender))
plot(u,pch=21,bg=class,xlab="First Component",ylab="Second Component")

class=as.numeric(as.factor(MetRef$donor))
kk=KODAMA.matrix(u)
tt=Rtsne(u,perplexity = 10)

par(mfrow=c(1,2))
plot(kk$scores,pch=21,bg=rainbow(22)[class], main="KODAMA-tSNE")
plot(tt$Y,pch=21,bg=rainbow(22)[class], main="tSNE, xlab="First Component", ylab="Second Component")

res_MDS=cmdscale(dist(u))
res_tSNE=Rtsne(u,perplexity = 20)$Y
custom.settings = umap.defaults
custom.settings$n_neighbors=20
res_UMAP = umap(u, config = custom.settings)$layout

kk=KODAMA.matrix(u)

res_KODAMA_MDS=KODAMA.visualization(kk,method = "MDS")
    
custom.settings = Rtsne.defaults
custom.settings$perplexity=20
res_KODAMA_tSNE=KODAMA.visualization(kk,method = "t-SNE",config = custom.settings)
    
custom.settings = umap.defaults
custom.settings$n_neighbors=20
res_KODAMA_UMAP=KODAMA.visualization(kk,method = "UMAP",config = custom.settings)

plot1 <-plot(res_MDS,pch=21,bg=rainbow(300),main="MDS")
plot2 <-plot(res_tSNE,pch=21,bg=rainbow(300),main="tSNE")
plot3 <-plot(res_UMAP,pch=21,bg=rainbow(300),main="UMAP")
    
plot4 <-plot(res_KODAMA_MDS,pch=21,bg=rainbow(300),main="KODAMA_MDS",ylim=range(res_KODAMA_MDS[,1]))
plot5 <-plot(res_KODAMA_tSNE,pch=21,bg=rainbow(300),main="KODAMA_tSNE")
plot6 <-plot(res_KODAMA_UMAP,pch=21,bg=rainbow(300),main="KODAMA_UMAP")

![This is an image](https://github.com/tkcaccia/Documents/blob/main/MetRef.png)

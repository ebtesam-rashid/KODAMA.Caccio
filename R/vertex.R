
vertex = function(vertix=c(0,10),dims=2){
  out=as.matrix(vertix)
  if(dims>1){
    for(i in 2:dims){
      nr=nrow(out)
      out=cbind(out,NA)
      out=rbind(out,out)
      
      out[1:nr,i]=vertix[1]
      out[1:nr+nr,i]=vertix[2]
    }
  }
  out=as.numeric(out)
  out
}

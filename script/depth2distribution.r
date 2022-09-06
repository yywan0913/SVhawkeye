GenomeDepthDistribution <- function(depth,output){
  library(data.table)

  df = fread(depth,header=F,sep="\t")
  df = as.data.frame(df)
  len = nrow(df)
  wlen = nchar(len)-3
  if(wlen>0){
      windows = as.numeric(paste0(1,paste0(rep(0,wlen),collapse="")))
  }else{
      windows = 1
  }
  chrom = df[1,1]
  start = df[1,2]
  end = df[nrow(df),2]

  df[,2] = ceiling(df[,2]/windows)
  df = aggregate(V3~V1+V2,data=df,mean)
  df[,2] = df[,2]*windows

  pdf(output,width=12,height=6)
  par(cex.axis=1.2,cex.lab=1.5,cex.main=1.8)
  plot(1,xlim=c(max(start-1,0),end+1),ylim=c(0,max(df[,3])),type="n",main=sprintf("%s:%s-%s",chrom,start,end),xlab="",ylab="Depth",bty="l",las=1,xaxt="n")
  rect(df[,2]-0.5,0,df[,2]+0.5,df[,3],col="grey70",border="grey70")
  a1 = axis(1,tick=FALSE,labels=FALSE)
  if(a1[length(a1)]>1000000) a1x = paste0(a1/1000000,"Mb") else 
      if(a1[length(a1)]>1000) a1x = paste0(a1/1000000,"Kb") else
	   a1x = a1
  axis(1,a1,a1x)
  dev.off()
}



args = commandArgs(TRUE)
depthfile = args[1]
outpdf = args[2]
GenomeDepthDistribution(depthfile,outpdf)

GenomeDepthDistribution <- function(depthf,i=1){
  library(data.table)

  df = fread(depthf,header=F,sep="\t")
  df = as.data.frame(df)
  len = nrow(df)
  windows = df[1,3]-df[1,2]
  chrom = df[1,1]
  start = df[1,2]
  end = df[nrow(df),3]
  df[,4] = df[,4]/(df[,3]-df[,2])
  #df[,2] = ceiling(df[,2]/windows)
  #df = aggregate(V3~V1+V2,data=df,mean)
  #df[,2] = df[,2]*windows
  main=sprintf("%s:%s-%s",chrom,start,end)
  xlab = ''
  ylab = gsub(sprintf(".%s_.*",chrom),'',basename(depthf))
  if(i==1) main= main else main = ""
  par(cex.axis=1.2,cex.lab=1.2,cex.main=1.5)
  plot(1,xlim=c(max(start-1,0),end+1),ylim=c(0,max(df[,ncol(df)])),type="n",main=main,xlab=xlab,ylab=ylab,bty="l",las=1,xaxt="n")
  rect(df[,2]-0.5*windows,0,df[,3]+0.5*windows,df[,4],col="grey80",border="grey90")
  a1 = axis(1,tick=FALSE,labels=FALSE)
  if(a1[length(a1)]>1000000) a1x = paste0(a1/1000000,"Mb") else 
      if(a1[length(a1)]>1000) a1x = paste0(a1/1000000,"Kb") else
	   a1x = a1
  axis(1,a1,a1x)
}



args = commandArgs(TRUE)
depthfile = args[1]
output = args[2] # out.pdf  out.png
options(digits = 10)
files = strsplit(depthfile,",")[[1]]
n = length(files)
if(substr(output,nchar(output)-3,nchar(output))==".pdf") pdf(output,width=12,height=6)
if(substr(output,nchar(output)-3,nchar(output))==".png") png(output,width=6000,height=3000,type="cairo",res=400)
par(mfrow=c(n,1))
for(i in 1:n){
  filei = files[i]
  par(mar=c(3.5,5,3,2))
  GenomeDepthDistribution(filei,i=i)
}
dev.off()



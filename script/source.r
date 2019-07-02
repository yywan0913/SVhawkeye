drawigv = function(data,MinInsert=10,MinINS=40,Insert_col = rgb(118/255,24/255,220/255),MinDEL=40,sample=NULL,TRA=FALSE,splitreads=NULL,left=TRUE,nsamples=NULL){
	data0 = data
	if (length(unique(data$Reads))>=99) Depth = "depth:>99" else Depth = paste0("depth:",length(unique(data$Reads)))
	#if (length(unique(data$Reads))>=99) sample = paste(sample,"depth:>99",sep="\n") else sample = paste(sample,paste0("depth:",length(unique(data$Reads))),sep="\n")
	if(nsamples<=2&nrow(data)>40) data=head(data,n=40) 
	if(nsamples>=3&nrow(data)>25) data=head(data,n=25)
	normcol = rgb(202/255,202/255,202/255) #"#9C9C9C"
	svcolor = "#9C9C9C" #"#BCD2EE"  #######
	peak_h =2
	d = 0.1
	refstart = data$RefStart
	refend = data$RefEnd
	strand = data$Strand
	colorindex = data$Color
	Reads = data$Reads
	type = data$Type
	if(TRA) split_color = rainbow(length(splitreads),s=0.4,v=0.4)  
	if(max(colorindex)>0) Reads_color = rainbow(max(colorindex),s=0.6,v=0.6)
	#Reads_col
	Reads_col = c()
	for(i in 1:length(colorindex)){
		j = colorindex[i]
		if(TRA){
			if(!is.na(match(Reads[i],splitreads))){
				Reads_col = c(Reads_col,split_color[match(Reads[i],splitreads)])
			}else{
				if(j==0) Reads_col=c(Reads_col,normcol) else{
					if(j== -1) Reads_col=c(Reads_col,svcolor) else Reads_col=c(Reads_col,Reads_color[j])
				}
			}
		}else{
			if(j==0) Reads_col=c(Reads_col,normcol) else{
				if(j== -1) Reads_col=c(Reads_col,svcolor) else Reads_col=c(Reads_col,Reads_color[j])
			}
		}

	}
	## region
	Data = data[,-c(1:match('Reads',colnames(data)))]
	region = as.numeric(colnames(Data))   ## Absolute value
	
	if(is.null(sample)) ylab = "" else ylab = sample
	par(cex.lab=2)
	ylim1 = ifelse(nrow(data)<8,8,nrow(data))
	plot(1,type="n",xlab="",ylab=ylab,axes=F,xlim=c(min(region)-1,max(region)+1),ylim=c(-ylim1-1,peak_h),xaxs="i",yaxs="i")
	par(xpd=T)
	yrow = 0
	for(i in 1:nrow(Data)){    # zong   i relatiuve value
		## Y coordinate
		DELsplit = FALSE
		INSsplit = FALSE
		if(i==1){yrow=yrow+1}else{
			if(Reads[i]==Reads[i-1] & grepl('DEL--',type[i]) & grepl('DEL--',type[i-1])){
				typei = gsub('DEL--','',gsub('@@.*','',gsub('.*(DEL--.*)','\\1',type[i])))
				typei_ = gsub('DEL--','',gsub('@@.*','',gsub('.*(DEL--.*)','\\1',type[i-1])))
				if(length(intersect(strsplit(typei,",")[[1]],strsplit(typei_,",")[[1]] ))==0){
					yrow = yrow + 1
				}else{
					yrow = yrow
					DELsplit = TRUE
				}
			}else if(Reads[i]==Reads[i-1] & grepl('INS--',type[i]) & grepl('INS--',type[i-1])){
				typei = gsub('INS--','',gsub('@@.*','',gsub('.*(INS--.*)','\\1',type[i])))
				typei_ = gsub('INS--','',gsub('@@.*','',gsub('.*(INS--.*)','\\1',type[i-1])))
				if(length(intersect(strsplit(typei,",")[[1]],strsplit(typei_,",")[[1]] ))==0){
					yrow = yrow + 1
				}else{
					yrow = yrow
					INSsplit = TRUE
				}
				
			}else{
				yrow = yrow + 1
			}
		}
		readbaselist = unlist(Data[i,])
		## bianjie
		string = readbaselist
		string[string==""]="+"
		strings = paste(string,collapse="",sep="+")
		strings = gsub('NA','+',strings)
		sleft = nchar(strings)-nchar(gsub('^\\++','',strings))
		sright = length(region) - (nchar(strings)-nchar(gsub('\\++$','',strings)))
		## strand arrows
		triangle = ifelse(TRA,1/80,1/150)
		arrow = ifelse(TRA,1/25,1/40)
		if(strand[i]=="+") {
			if(data$RefEnd[i]<=max(region)){
				polygon(x=c(region[sright],region[sright]+triangle*(par('usr')[2]-par('usr')[1]),region[sright]),y=c(-yrow-1+d,-yrow-0.5,-yrow-d),col=Reads_col[i],border=Reads_col[i])
			}else{
				arrows(par('usr')[2], -yrow-0.5, par('usr')[2]+arrow*(par('usr')[2]-par('usr')[1]), -yrow-0.5,col="grey",code=2,lwd=2,length=0.1)
			}
		}else{
			if(data$RefStart[i]>=min(region)){
				polygon(x=c(region[sleft+1],region[sleft+1]-triangle*(par('usr')[2]-par('usr')[1]),region[sleft+1]),y=c(-yrow-1+d,-yrow-0.5,-yrow-d),col=Reads_col[i],border=Reads_col[i])
			}else{
				arrows(par('usr')[1],-yrow-0.5,par('usr')[1]-arrow*(par('usr')[2]-par('usr')[1]),-yrow-0.5,col="grey",code=2,lwd=2,length=0.1)
			}
		}
		## base
		if(ncol(Data)>=20000){
			if(grepl('DEL_',type[i])){
				typei = gsub('@@.*','',gsub('.*(DEL_.*)','\\1',type[i]))
				pos_len = strsplit(gsub('DEL_','',typei),",")[[1]]
				pos = as.numeric(sapply(strsplit(pos_len,':'),'[',1))  ## centerpos 
				len = as.numeric(sapply(strsplit(pos_len,':'),'[',2))
				for(p in 1:length(pos)){
					pleft = pos[p]-ceiling(1/2*len[p])
					pright = pos[p]+trunc(1/2*len[p])
					if(p==1){
						pleftp = sleft+1
						rect(region[pleftp]-0.5,-yrow-1+d,region[pleft]+0.5,-yrow-d,col=Reads_col[i],border=Reads_col[i])
						segments(region[pleft]-0.5,-yrow-0.5,region[pright]+0.5,-yrow-0.5,col=svcolor,lwd=2)
						pleftp = pright
						
					}else{
						rect(region[pleftp]-0.5,-yrow-1+d,region[pleft]+0.5,-yrow-d,col=Reads_col[i],border=Reads_col[i])
						segments(region[pleft]-0.5,-yrow-0.5,region[pright]+0.5,-yrow-0.5,col=svcolor,lwd=2)
						pleftp = pright
					}
					if(p==length(pos)) {
						rect(region[pleftp]-0.5,-yrow-1+d,region[sright]+0.5,-yrow-d,col=Reads_col[i],border=Reads_col[i])
					}
					
				}
			}else{
				rect(region[sleft+1],-yrow-1+d,region[sright],-yrow-d,col=Reads_col[i],border=Reads_col[i])
			}
			
		}else{
			for(j in 1:ncol(Data)){   ## heng  X coordinate
				if((readbaselist[j]==""|readbaselist[j]=="+")&(j<=sleft|j>sright)) next
				if((readbaselist[j]==""|readbaselist[j]=="+")&(j>sleft&j<=sright)) {
					segments(region[j]-0.5,-yrow-0.5,region[j]+0.5,-yrow-0.5,lwd=2,col=normcol);next
				}
				if(readbaselist[j]=="-") {
					segments(region[j]-0.5,-yrow-0.5,region[j]+0.5,-yrow-0.5,lwd=2,col=svcolor) 
				}else{ 
					rect(region[j]-0.5,-yrow-1+d,region[j]+0.5,-yrow-d,col=Reads_col[i],border=Reads_col[i])
				}
				if(nchar(readbaselist[j])>=MinInsert&&nchar(readbaselist[j])<MinINS) {
					rect(region[j]-0.5,-yrow-1+d,region[j]+0.5,-yrow-d,col=Insert_col,border=Insert_col)
				}
			}
		}
		## ins or del labels
		if(grepl('DEL_',type[i])){
			typei = gsub('@@.*','',gsub('.*(DEL_.*)','\\1',type[i]))
			pos_len = strsplit(gsub('DEL_','',typei),",")[[1]]
			pos = as.numeric(sapply(strsplit(pos_len,':'),'[',1))
			len = as.numeric(sapply(strsplit(pos_len,':'),'[',2))
			text(region[pos],-yrow-0.5,len,cex=1,col=rgb(118/255,24/255,220/255),font=2)
		}
		if(grepl('INS_',type[i])){
			typei = gsub('@@.*','',gsub('.*(INS_.*)','\\1',type[i]))
			pos_len = strsplit(gsub('INS_','',typei),",")[[1]]
			pos = as.numeric(sapply(strsplit(pos_len,':'),'[',1))
			len = as.numeric(sapply(strsplit(pos_len,':'),'[',2))
			rect(region[pos]-0.5,-yrow-1+d,region[pos]+0.5+strwidth(len),-yrow-d,col=Insert_col,border=Insert_col)
			text(region[pos],-yrow-0.5,len,cex=1,adj=0,col="white",font=2)
		}
		## DEL split mapping 
		if(DELsplit){
			if(!grepl('INV',type[i]) & !grepl('INV',type[i-1])){
			#if(strand[i]=="+"){
				delx1 = region[match(refend[i-1],region)]
				delx2 = region[match(refstart[i],region)]
			#}else{
			#	delx1 = region[match(refend[i],region)]
			#	delx2 = region[match(refstart[i-1],region)]
			#}
			segments(delx1,-yrow-0.5,delx2,-yrow-0.5,col=Reads_col[i],lty=2,lwd=2)
			text((delx1+delx2)/2,-yrow-0.5,abs(refstart[i]-refend[i-1]),cex=1,col=rgb(118/255,24/255,220/255),font=2)
			}
		}
		## INS split mapping
		if(INSsplit){
			typei = gsub('@@.*','',gsub('.*(INS--.*)','\\1',type[i]))
			len = strsplit(gsub('INS--','',typei),",")[[1]]
			delx1 = region[match(refend[i-1],region)]
			delx2 = region[match(refstart[i],region)]
			pos = (delx1+delx2)/2
			rect(pos-0.5,-yrow-1+d,pos+0.5+strwidth(len),-yrow-d,col=Insert_col,border=Insert_col)
			text(pos,-yrow-0.5,len,cex=1,adj=0,col="white",font=2)
		}
	}
	## kurtosis depth or coverage
	y_peakbottom = 0
	peakdata = data0[,-c(1:match('Reads',colnames(data0)))]
	peakdata = ifelse(is.na(peakdata)|peakdata=="-"|peakdata==""|peakdata=="+",0,1)
	peak = apply(peakdata,2,sum)
	peak = peak/max(peak)*peak_h
	rect(region-0.5,y_peakbottom,region+0.5,y_peakbottom+peak,col="grey",border="#563624")
	#axis(1)
	par(xpd=T)
	if(left){
		#segments(par('usr')[1]-1/3*(par('usr')[2]-par('usr')[1]),par('usr')[4],par('usr')[1],par('usr')[4],lwd=2,lty=2)
		text(par('usr')[1],y_peakbottom+peak_h/2,Depth,adj=1,cex=1.5,font=2)
	}else{
		text(par('usr')[2],y_peakbottom+peak_h/2,Depth,adj=0,cex=1.5,font=2)
	}
}



Readcytoband = function(file,chrom){
        options(stringsAsFactors=F)
        cyto = read.table(file,header=F,sep="\t")
        colnames(cyto)<-c("chr","start","end","name","type")
        cyto$Color = "white"
        cyto[cyto$type=="gneg",]$Color<-rgb(255,255,255, maxColorValue=255)
        cyto[cyto$type=="gpos25",]$Color<-rgb(200,200,200, maxColorValue=255)
        cyto[cyto$type=="gpos50",]$Color<-rgb(150,150,150, maxColorValue=255)
        cyto[cyto$type=="gpos75",]$Color<-rgb(130,130,130, maxColorValue=255)
        cyto[cyto$type=="gpos100",]$Color<-rgb(100,100,100, maxColorValue=255)
        cyto[cyto$type=="acen",]$Color<-"red" # red centromere 
        cyto[cyto$type=="stalk",]$Color<-rgb(100,127,164, maxColorValue=255) # repeat regions 
        cyto[cyto$type=="gvar",]$Color<-rgb(220,220,220, maxColorValue=255) # indented region
        cyto$chr = gsub('chr','',cyto$chr)
        cyto = cyto[cyto$chr==chrom,]

        return(cyto)
}

drawcytoband = function(cytoband,chrom,start,end,xlim,main,genome,left=TRUE){
	cyto = Readcytoband(file=cytoband,chrom=chrom)
	bei = (max(cyto$end)-min(cyto$start))/(end-start)
	cyto$start = cyto$start/bei +start
	cyto$end = cyto$end/bei +start
	cstart = start/bei + start
	cend = end/bei +start
	
	main=paste0(main,"(",prettyNum(abs(end-start),","),"bp)")
	par(cex.axis=1.5,cex.lab=2,lab=c(8,5,0))
	plot(3:6,xlim=xlim,type="n",xlab=main,ylab="",xaxs="i",yaxs="i",axes=F)
	abline(h=par('usr')[3])
	text((par('usr')[1]+par('usr')[2])/2,5.5,paste0("chromosome:",chrom),cex=1.5,font=2)
	axis(1)
	par(xpd=T)
	if(nrow(cyto)>0){
		xacen = grep('acen',cyto$type)
		cytoi1 = cyto[1:xacen[1],]
		cytoi2 = cyto[xacen[2]:nrow(cyto),]
		cytoacenF = cyto[cyto$type!="acen",]
		cytoacen = cyto[cyto$type=="acen",]
		cytoacenF$Color[1]=NA
		cytoacenF$Color[nrow(cytoacenF)] = NA
		rect(cytoacenF$start,4,cytoacenF$end,5,col=cytoacenF$Color,border="NA",lwd=2)

		roundrect(min(cytoi1$start),4,max(cytoi1$end),5,strand=1,col_1=NA,border_1=NA,col_2='red',border_2='red',lwd=2,col_middle=NA,border_middle="black",circlelength=cytoacen$end[1]-cytoacen$start[1])
                roundrect(min(cytoi2$start),4,max(cytoi2$end),5,strand=1,col_1='red',border_1='red',col_2=NA,border_2=NA,lwd=2,col_middle=NA,border_middle="black",circlelength=cytoacen$end[2]-cytoacen$start[2])
		par(xpd=T)
        	for(i in 1:nrow(cyto)){
                	if(i%%4==3 & strwidth(cyto$name[i],cex=1) > cyto$end[i]-cyto$start[i]) {
                	        text((cyto$start[i]+cyto$end[i])/2,4-0.5,cyto$name[i],cex=1)
                	}
                	if(strwidth(cyto$name[i],cex=1) <= cyto$end[i]-cyto$start[i]) text((cyto$start[i]+cyto$end[i])/2,4.5,cyto$name[i],cex=1)
        	}
        rect(cstart,4,cend,5,col="darkred",border="darkred")
	par(xpd=T)
	if(left) text(par('usr')[1],4.5,genome,adj=1,cex=2)
}


}


GetregionIntersect = function(start1,end1,start2,end2){
        TF = !( end1<start2 | start1>end2)
        return(TF)
}

getannot = function(scriptdir,genome,chrom,start,end){
	annot = list()
	databasedir = paste0(dirname(scriptdir),"/database")
	if(genome=="hg19"){
        	segdup = gettextf("gzip -dc %s",paste0(databasedir,"/hg19/hg19_genomicSuperDups.bed.gz"))
        	rmsk = gettextf("gzip -dc %s",paste0(databasedir,"/hg19/hg19_rmsk.bed.gz"))
        	cytoband = paste0(databasedir,"/hg19/cytoBand.txt")
        	refgene = gettextf("gzip -dc %s",paste0(databasedir,"/hg19/hg19.bed.gtf.gz"))
	}else{
        	segdup = gettextf("gzip -dc %s",paste0(databasedir,"/hg38/hg38_genomicSuperDups.bed.gz"))
        	rmsk = gettextf("gzip -dc %s",paste0(databasedir,"/hg38/hg38_rmsk.bed.gz"))
        	cytoband = paste0(databasedir,"/hg38/cytoBand.txt")
        	refgene = gettextf("gzip -dc %s",paste0(databasedir,"/hg38/hg38.bed.gtf.gz"))
	}
	Segdup = fread(segdup,header=F,sep="\t")
	Segdup = as.data.frame(Segdup)
	Segdup[,1] = gsub('chr','',Segdup[,1])
	Segdup = Segdup[Segdup[,1]==chrom,]
	Segdup = Segdup[GetregionIntersect(Segdup[,2],Segdup[,3],start,end),,drop=F]

	Rmsk = fread(rmsk,sep="\t",header=F)
	Rmsk = as.data.frame(Rmsk)
	Rmsk[,1] = gsub('chr','',Rmsk[,1])
	Rmsk = Rmsk[Rmsk[,1]==chrom,]
	Rmsk = Rmsk[GetregionIntersect(Rmsk[,2],Rmsk[,3],start,end),,drop=F]

	Refgene = fread(refgene,sep="\t",header=F)
	Refgene = as.data.frame(Refgene)
	Refgene[,1] = gsub('chr','',Refgene[,1])
	Refgene = Refgene[Refgene[,1]==chrom,]
	Refgene = Refgene[GetregionIntersect(Refgene[,2],Refgene[,3],start,end),,drop=F]

	annot[['Refgene']] = Refgene
	annot[['Rmsk']] = Rmsk
	annot[['Segdup']] = Segdup
	annot[['cytoband']] = cytoband
	
	return (annot)
}


overlap = function(start1,end1,start2,end2){
        start = apply(data.frame(start1,start2),1,max)
        end = apply(data.frame(end1,end2),1,min)
        return(data.frame(start=start,end=end))
}


recttext = function(rmskintersect,textcex=0.6,up,down,start,end){
	par(xpd=F)
	rect(rmskintersect[,2],down,rmskintersect[,3],up,col="blue",border="blue")
        par(xpd=T)
        textpos = apply(overlap(rmskintersect[,2],rmskintersect[,3],start,end),1,sum)/2
	Text = sapply(strsplit(as.character(rmskintersect[,4]),"-"),"[",1)
	textcex = textcex
        textwidth = strwidth(Text,cex=textcex)
	nextn = 0
        for(k in 1:length(textpos)){
               textkregion = c(textpos[k]-textwidth[k]/2,textpos[k]+textwidth[k]/2)
               if(k==1) {
                     text(textpos[k],down-0.5, Text[k],cex=textcex,font=2)
                     textregion = textkregion
               }else{
                     if(GetregionIntersect(textregion[1],textregion[2],textkregion[1],textkregion[2])){
			nextn = nextn +1
			if(nextn==1){
				text(textpos[k],down-1, Text[k],cex=textcex,font=2)
				textregion = textregion
				text2region = textregion
			}else{
				if(!GetregionIntersect(text2region[1],text2region[2],textkregion[1],textkregion[2])){
                        		text(textpos[k],down-1, Text[k],cex=textcex,font=2)
                        		textregion = textregion
					text2region = textregion
				}
			}
                     }else{
                        text(textpos[k],down-0.5, Text[k],cex=textcex,font=2)
                        textregion = c(min(textregion[1],textkregion[1]),max(textregion[2],textkregion[2]))
                    }
               }
        }
}



drawannot = function(xlim,Annot,start,end,TRA=FALSE){
	plot(-8:10,xlim=xlim,type="n",xlab="",ylab="",xaxs="i",yaxs="i",axes=F)
	## refgene
	refgeneintersect = Annot[['Refgene']]
	if(nrow(refgeneintersect)>0){
        	par(xpd=F)
        	rect(refgeneintersect[,2],8,refgeneintersect[,3],9,col="darkblue",border="darkblue")
        	par(xpd=T)
		textpos = apply(overlap(refgeneintersect[,2],refgeneintersect[,3],start,end),1,sum)/2
        	#text(textpos,7, refgeneintersect[,4],cex=1.2,font=2)
		Text = refgeneintersect[,4]
		textcex = 0.8
		textwidth = strwidth(Text,cex=textcex)
		for(k in 1:length(textpos)){
			textkregion = c(textpos[k]-textwidth[k]/2,textpos[k]+textwidth[k]/2)
			if(k==1) {
				text(textpos[k],7.5, Text[k],cex=1.2,font=2)
				textregion = textkregion
			}else{
				if(GetregionIntersect(textregion[1],textregion[2],textkregion[1],textkregion[2])){
					text(textpos[k],7, Text[k],cex=1.2,font=2)
					textregion = textregion
				}else{
					text(textpos[k],7.5, Text[k],cex=1.2,font=2)
					textregion = c(min(textregion[1],textkregion[1]),max(textregion[2],textkregion[2]))
				}
			}
		}
	}
	par(xpd=T)
	if(!TRA) text(par('usr')[1],8.5,"RefGene",adj=1,cex=1.5)

	## --- rmsk
	rmskintersect = Annot[['Rmsk']]
	if(nrow(rmskintersect)>0){
		RMSK = Text = sapply(strsplit(as.character(rmskintersect[,4]),"-"),"[",2)
		if(any(grepl('LINE',RMSK))) {
			LINE = rmskintersect[RMSK=="LINE",,drop=F]
			recttext(LINE,textcex=0.8,6,5,start,end)
		}
		if(any(grepl('SINE',RMSK))) {
			SINE = rmskintersect[RMSK=="SINE",,drop=F]
			recttext(SINE,textcex=0.8,3,2,start,end)
		}
		if(any(grepl('Simple_repeat',RMSK))){
			Simple = rmskintersect[RMSK=="Simple_repeat",,drop=F]
			recttext(Simple,textcex=0.8,0,-1,start,end)
		}
		Other = rmskintersect[RMSK!="LINE"&RMSK!="SINE"&RMSK!="Simple_repeat",,drop=F]
		if(nrow(Other)>0){
			recttext(Other,textcex=0.8,-3,-4,start,end)
		}

	}
	par(xpd=T)
	if(!TRA)text(par('usr')[1],5.5,"LINE",adj=1,cex=1.3)
	if(!TRA)text(par('usr')[1],2.5,"SINE",adj=1,cex=1.3)
	if(!TRA)text(par('usr')[1],-0.5,"Simple_repeat",adj=1,cex=1.3)
	if(!TRA)text(par('usr')[1],-3.5,"Other",adj=1,cex=1.3)
	#if(!TRA)text(par('usr')[1],5.5,"Rmsk",adj=1,cex=2)

	##--- segdup
	Segdupintersect = Annot[['Segdup']]
	if(nrow(Segdupintersect)>0){
        	par(xpd=F)
        	rect(Segdupintersect[,2],-7,Segdupintersect[,3],-6,col="darkgreen",border="darkgreen")
	}
	par(xpd=T)
	if(!TRA)text(par('usr')[1],-6.5,"Segdup",adj=1,cex=1.3)
}


roundrect = function(x1,y1,x2,y2,strand=2,col_1=NULL,border_1=NA,col_2=NULL,border_2=NA,angle=160,circlelength=NULL,col_middle=NA,border_middle='black',lwd=1,...){
	## col_1  top/left
	## col_2  bottom/right
	usr = par('usr')
	bei = (usr[4]-usr[3]) / (usr[2] - usr[1])
	
	if(strand==2){
		xa = 1/2*abs((x2-x1))  ## length for x
		ya = xa*bei       ## length for y
	}
	if(strand==1){
		ya = 1/2*abs((y2-y1))
		xa = ya/bei
	}
	stopifnot(strand!=1|strand!=2)

	if(is.null(circlelength)){
		angle = ifelse(angle>180,180,angle)
		thetaAngle = 1/2*angle
		theta = 1/2*angle/180*pi
	}else{
		if(strand==2){
			theta = asin((2*ya*circlelength)/(ya^2+circlelength^2))
		}else{
			theta = asin((2*xa*circlelength)/(xa^2+circlelength^2))
		}
		thetaAngle = theta/pi*180
		angle = 2*thetaAngle
	}
	## bottom / left
	Ry = ya/sin(theta)  ## radius for y  length
	Rx = xa/sin(theta)  ## radius for x  length
	if(strand==2){
		r1_x = 1/2*(x1+x2)  ## Center of mind for x  bottom  coordinate
		r1_y = 	y1 + Ry     ## Center of mind for y  bottom
		radian1 = (180+(90-thetaAngle) ) : (180 + (90-thetaAngle)+angle) 
		radian1 = radian1/180 * pi  ## circle theta
		x1_radian1 = r1_x+Rx*cos(radian1) 
		y1_radian1 = r1_y+Ry*sin(radian1)
		x11 = x1  ##
		y11 = r1_y-Ry*cos(theta)   ## y1+circlelength
	}else{
		r1_x = x1 + Rx
		r1_y = 1/2 * (y1+y2)
		radian1 = (180-thetaAngle) : (180+thetaAngle)
		radian1 = radian1/180 * pi
		x1_radian1 = r1_x+Rx*cos(radian1)
		y1_radian1 = r1_y+Ry*sin(radian1)
		x11 = r1_x - Rx*cos(theta)
		y11 = y1
		
	}

	##  top / right
	if(strand==2){
		r2_x = r1_x
		r2_y = y2-Ry
		radian2 = (90-thetaAngle):(90-thetaAngle+angle)
		radian2 = radian2/180 * pi
		x2_radian2 = r2_x + Rx*cos(radian2)
		y2_radian2 = r2_y+Ry*sin(radian2)

		x22 = x2
		y22 = r2_y + Ry*cos(theta)
	}else{
		r2_x = x2 - Rx
		r2_y = r1_y
		radian2 = -thetaAngle:thetaAngle
		radian2 = radian2/180 * pi
		x2_radian2 = r2_x + Rx*cos(radian2)
		y2_radian2 = r2_y+Ry*sin(radian2)
		x22 = r2_x + Rx*cos(theta)
		y22 = y2
	
	}

	##
	if(strand==2){
		all_x = c(x11,x1_radian1,x22,x22,x2_radian2,x11,x11)
		all_y = c(y11,y1_radian1,y11,y22,y2_radian2,y22,y11)

		all_xbottom = c(x11,x1_radian1,x22,x11)
		all_ybottom = c(y11,y1_radian1,y11,y11)
		all_xtop = c(x22,x2_radian2,x11,x22)
		all_ytop = c(y22,y2_radian2,y22,y22)
		## zong
		#polygon(x=c(x22,x2_radian2,x11,x11,x1_radian1,x22,x22),y=c(y22,y2_radian2,y22,y11,y1_radian1,y11,y22),col=NA,border=border_middle,lwd=lwd,...)
		polygon(x=all_x,y=all_y,col=col_middle,border=border_middle,lwd=lwd,...)
		## top
		polygon(x=all_xtop,y=all_ytop,col=col_1,border=border_1,lwd=lwd,...)
		## bottom
		polygon(x=all_xbottom,y=all_ybottom,col=col_2,border=border_2,lwd=lwd,...)
		## cencer
		#rect(x11,y11,x22,y22,col=col_middle,border=NA,lwd=lwd,...)
	}else{
		all_x = c(x11,x1_radian1,x11,x22,x2_radian2,x22,x11)
		all_y = c(y22,y1_radian1,y11,y11,y2_radian2,y22,y22)
		all_xleft = c(x11,x1_radian1,x11,x11)
		all_yleft = c(y22,y1_radian1,y11,y22)
		all_xright = c(x22,x2_radian2,x22,x22)
		all_yright = c(y11,y2_radian2,y22,y11)
		polygon(x=all_x,y=all_y,col=col_middle,border=border_middle,lwd=lwd,...)
		polygon(x=all_xleft,y=all_yleft,col=col_1,border=border_1,lwd=lwd,...)
		polygon(x=all_xright,y=all_yright,col=col_2,border=border_2,lwd=lwd,...)
	}

}



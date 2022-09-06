set.seed(123456)
mycolors <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#7E6148","#B09C85","#631879","#5F559B","#A20056","#FFDC91","#EE4C97","#00468B","#925E9F","#FDAF91","#ADB6B6","#374E55","#DF8F44","#6A6599","#80796B","#EFC000","#CD534C","#7AA6DC","#8F7700","#3B3B3B","#4A6990","#6699FF","#CC33FF","#99991E","#FF00CC","#FFCCCC","#99CCFF","#CCFFFF","#9900CC","#996600","#CCCCCC","#79CC3D","#CCCC99","#AEC7E8","#FFBB78","#98DF8A","#FF9896","#C5B0D5","#C49C94","#F7B6D2","#C7C7C7","#DBDB8D","#9EDAE5","#393B79","#637939","#8C6D31","#843C39","#7B4173","#5254A3","#8CA252","#BD9E39","#AD494A","#A55194","#6B6ECF","#B5CF6B","#E7BA52","#D6616B","#CE6DBD","#9C9EDE","#CEDB9C","#E7CB94","#E7969C","#DE9ED6","#3182BD","#E6550D","#31A354","#756BB1","#636363","#6BAED6","#FD8D3C","#74C476","#9E9AC8","#969696","#9ECAE1","#FDAE6B","#A1D99B","#BCBDDC","#BDBDBD","#C6DBEF","#FDD0A2","#C7E9C0","#DADAEB","#D9D9D9","#5050FF","#CE3D32","#749B58","#F0E685","#466983","#BA6338","#5DB1DD","#802268","#6BD76B","#D595A7","#924822","#837B8D","#C75127","#D58F5C","#7A65A5","#E4AF69","#3B1B53","#CDDEB7","#612A79","#AE1F63","#E7C76F","#5A655E","#CC9900","#99CC00","#A9A9A9","#CC9900","#99CC00","#33CC00","#00CC33","#00CC99","#0099CC","#0A47FF","#4775FF","#FFC20A","#FFD147","#990033","#991A00","#996600","#809900","#339900","#00991A","#009966","#008099","#003399","#1A0099","#660099","#990080","#D60047","#FF1463","#00D68F","#14FFB1","#5773CC","#FFB900","#D43F3A","#EEA236","#5CB85C","#46B8DA","#357EBD","#9632B8","#B8B8B8","#800000","#767676","#FFA319","#8A9045","#155F83","#C16622","#8F3931","#58593F","#350E20","#D6D6CE","#FFB547","#ADB17D","#5B8FA8","#D49464","#B1746F","#8A8B79","#725663","#CC8214","#616530","#0F425C","#9A5324","#642822","#3E3E23","#FED439","#709AE1","#8A9197","#D2AF81","#FD7446","#D5E4A2","#197EC0","#F05C3B","#46732E","#71D0F5","#370335","#075149","#C80813","#91331F","#1A9993","#FD8CC1","#FF6F00","#C71000","#008EA0","#8A4198","#5A9599","#FF6348","#84D7E1","#FF95A8","#3D3B25","#ADE2D0","#1A5354","#3F4041","#FAFD7C","#82491E","#24325F","#B7E4F9","#FB6467","#526E2D","#E762D7","#E89242","#FAE48B","#A6EEE6","#917C5D","#69C8EC","#CC0C00","#5C88DA","#84BD00","#FFCD00","#7C878E","#00B5E2","#00AF66","#FF410D","#6EE2FF","#F7C530","#95CC5E","#D0DFE6","#F79D1E","#748AA6","#4500AD","#2700D1","#6B58EF","#8888FF","#C7C1FF","#D5D5FF","#FFC0E5","#FF8989","#FF7080","#FF5A5A","#EF4040","#D60C00","#E7F0FA","#C9E2F6","#95CBEE","#0099DC","#4AB04A","#FFD73E","#EEC73A","#E29421","#E29421","#F05336","#CE472E","#FFEBEE","#FFCDD2","#EF9A9A","#E57373","#EF5350","#F44336","#E53935","#D32F2F","#C62828","#B71C1C","#FCE4EC","#F8BBD0","#F48FB1","#F06292","#EC407A","#E91E63","#D81B60","#C2185B","#AD1457","#880E4F","#F3E5F5","#E1BEE7","#CE93D8","#BA68C8","#AB47BC","#9C27B0","#8E24AA","#7B1FA2","#6A1B9A","#4A148C","#EDE7F6","#D1C4E9","#B39DDB","#9575CD","#7E57C2","#673AB7","#5E35B1","#512DA8","#4527A0","#311B92","#E8EAF6","#C5CAE9","#9FA8DA","#7986CB","#5C6BC0","#3F51B5","#3949AB","#303F9F","#283593","#1A237E","#E3F2FD","#BBDEFB","#90CAF9","#64B5F6","#42A5F5","#2196F3","#1E88E5","#1976D2","#1565C0","#0D47A1","#E1F5FE","#B3E5FC","#81D4FA","#4FC3F7","#29B6F6","#03A9F4","#039BE5","#0288D1","#0277BD","#01579B","#E0F7FA","#B2EBF2","#80DEEA","#4DD0E1","#26C6DA","#00BCD4","#00ACC1","#0097A7","#00838F","#006064","#E0F2F1","#B2DFDB","#80CBC4","#4DB6AC","#26A69A","#009688","#00897B","#00796B","#00695C","#004D40","#E8F5E9","#C8E6C9","#A5D6A7","#81C784","#66BB6A","#4CAF50","#43A047","#388E3C","#2E7D32","#1B5E20","#F1F8E9","#DCEDC8","#C5E1A5","#AED581","#9CCC65","#8BC34A","#7CB342","#689F38","#558B2F","#33691E","#F9FBE7","#F0F4C3","#E6EE9C","#DCE775","#D4E157","#CDDC39","#C0CA33","#AFB42B","#9E9D24","#827717","#FFFDE7","#FFF9C4","#FFF59D","#FFF176","#FFEE58","#FFEB3B","#FDD835","#FBC02D","#F9A825","#F57F17","#FFF8E1","#FFECB3","#FFE082","#FFD54F","#FFCA28","#FFC107","#FFB300","#FFA000","#FF8F00","#FF6F00","#FFF3E0","#FFE0B2","#FFCC80","#FFB74D","#FFA726","#FF9800","#FB8C00","#F57C00","#EF6C00","#E65100","#FBE9E7","#FFCCBC","#FFAB91","#FF8A65","#FF7043","#FF5722","#F4511E","#E64A19","#D84315","#BF360C","#EFEBE9","#D7CCC8","#BCAAA4","#A1887F","#8D6E63","#795548","#6D4C41","#5D4037","#4E342E","#3E2723","#FAFAFA","#F5F5F5","#EEEEEE","#E0E0E0","#BDBDBD","#9E9E9E","#757575","#616161","#424242","#212121","#ECEFF1","#CFD8DC","#B0BEC5","#90A4AE","#78909C","#607D8B","#546E7A","#455A64","#37474F","#263238")

segmentsbend <- function(xpos1,y1,y2,strand,refstart,refend,TRA=FALSE){ ## sc softclipping > |  <
    if(!TRA) xdiff = 1/350*(par('usr')[2]-par('usr')[1]) else xdiff = 1/150*(par('usr')[2]-par('usr')[1])
    if(strand=="+"){
        if(refstart==xpos1) segments(xpos1-0.5,y1,xpos1-0.5,y2,col="red",lwd=1) else{
	if(refend==xpos1) {
		segments(xpos1+0.5-xdiff,y1,xpos1+0.5,(y1+y2)/2,col="red",lwd=1)
                segments(xpos1+0.5-xdiff,y2,xpos1+0.5,(y1+y2)/2,col="red",lwd=1)
	}
	}
    }else{
        if(refend==xpos1) segments(xpos1+0.5,y1,xpos1+0.5,y2,col="red",lwd=1) else{
	if(refstart==xpos1){
	    segments(xpos1-0.5+xdiff,y1,xpos1-0.5,(y1+y2)/2,col="red",lwd=1)
	    segments(xpos1-0.5+xdiff,y2,xpos1-0.5,(y1+y2)/2,col="red",lwd=1)
        }
   }
}
}

trirect <- function(xleft,ybottom,xright,ytop, strand, col, start,end,TRA=FALSE,border="white",...){ # rect
    if(!TRA) xdiff = 1/350*(par('usr')[2]-par('usr')[1]) else xdiff = 1/150*(par('usr')[2]-par('usr')[1])
    d = 0.15
    for(i in 1:length(xleft)){
	xleft[i] = ifelse(xleft[i]<start,start,xleft[i])
        xright[i] = ifelse(xright[i]>end,end,xright[i])
        if(strand[i]=="+"){
	    #xi = c(xleft[i],xleft[i],xright[i]-xdiff,xright[i],xright[i]-xdiff,xleft[i])
	    #yi = c(ybottom[i],ytop[i],ytop[i],(ybottom[i]+ytop[i])/2,ybottom[i],ybottom[i])
	    xi = c(xleft[i],xleft[i],xright[i]-xdiff,xright[i]-xdiff,xright[i],xright[i]-xdiff,xright[i]-xdiff,xleft[i])
	    yi = c(ybottom[i],ytop[i],ytop[i],ytop[i]+d,(ybottom[i]+ytop[i])/2,ybottom[i]-d,ybottom[i],ybottom[i])
        }else{
	    #xi = c(xleft[i]+xdiff,xleft[i],xleft[i]+xdiff,xright[i],xright[i],xleft[i]+xdiff)
	    #yi = c(ybottom[i],(ybottom[i]+ytop[i])/2,ytop[i],ytop[i],ybottom[i],ybottom[i])
	    xi = c(xleft[i]+xdiff,xleft[i]+xdiff,xleft[i],xleft[i]+xdiff,xleft[i]+xdiff,xright[i],xright[i],xleft[i]+xdiff)
	    yi = c(ybottom[i],ybottom[i]-d,(ybottom[i]+ytop[i])/2,ytop[i]+d,ytop[i],ytop[i],ybottom[i],ybottom[i])
        }
	polygon(x=xi,y=yi,col = col[i],border=col[i],...)
	if(strand[i]=="+") segments(xleft[i],ybottom[i],xleft[i],ytop[i],col=border) else
	    segments(xright[i],ybottom[i],xright[i],ytop[i],col=border)
    }
}
rectarrow <- function(xleft,ybottom,xright,ytop,strand,col,start,end,...){  ## may need rm 
    xdiff = 1/15*(par('usr')[2]-par('usr')[1])
    #rect(xleft,ybottom,xright,ytop,col=col,border=NA)
    arrowscol = rgb(230/255,230/255,230/255)
    #segments(xleft, (ybottom+ytop)/2, xright, (ybottom+ytop)/2,col = arrowscol)
    for(i in 1:length(xleft)){
        xleft[i] = ifelse(xleft[i]<start,start,xleft[i])
        xright[i] = ifelse(xright[i]>end,end,xright[i])
	if(xright[i]-xleft[i] < 3*xdiff) {
	    trirect(xleft=xleft[i],ybottom=ybottom[i],xright=xright[i],ytop=ytop[i], strand=strand[i], col=col[i])
	}else{
	    rect(xleft[i],ybottom[i],xright[i],ytop[i],col=col[i],border=NA)
	    xa1 = seq(xleft[i],xright[i]-xdiff,by=xdiff)
            xa2 = seq(xleft[i]+xdiff,xright[i],by=xdiff)
            ya = (ybottom[i]+ytop[i])/2
	    codei = ifelse(strand[i]=="+",2,1)
            arrows(xa1, ya, xa2, ya,col=arrowscol,code=codei,lwd=1,length=0.1,angle=20)
	}
    }
}

getrandomreads = function(data,n=40){
    if(max(data$Readsorder)<=n) return(data)
    randomorder = sample(max(data$Readsorder), n )
    data = data[data$Readsorder%in%randomorder,,drop=F]
    data$Readsorder = as.numeric(as.factor(data$Readsorder)) # reset order
    return(data)
}

drawigv = function(data,MinInsert=10,MinINS=40,Insert_col = rgb(118/255,24/255,220/255),MinDEL=40,
		sample=NULL, TRA=FALSE, splitreads=NULL, left=TRUE, nsamples=NULL, start, end,
		normcol = rgb(222/255,222/255,222/255),svcolor = rgb(176/255,176/255,176/255),
		peakcol=rgb(192/255,192/255,192/255),rnacolor=rgb(247/255,184/255,129/255)){
	library("RColorBrewer")
	data <- data[order(data$Readsorder),]
	data$Readsorder <- trunc(data$Readsorder)

	data0 <- data # for draw depth
	if(nsamples == 1) data <- getrandomreads(data,n=50) else
	    if(nsamples == 2) data <- getrandomreads(data,n=40) else
	        if(nsamples >= 3) data <- getrandomreads(data,n=25) 
	#peak_h <- 2 # depth distribution
	d <- 0.15
	refstart = data$RefStart 
	refend = data$RefEnd
	strand = data$Strand
	colorindex = data$Color
	Reads = data$ReadsID
	type = data$Type

	if(TRA) {
		#split_colorL = rainbow(length(splitreads),s=0.6,v=0.6)
		#if(length(split_colorL)>7) split_colorL = split_colorL[c(unlist(sapply(1:7,seq,length(split_colorL),by=7)))]
		split_color = colorRampPalette(brewer.pal(9, "Paired"))(length(splitreads))
		if(length(split_color)>7) split_color = split_color[c(unlist(sapply(1:7,seq,length(split_color),by=7)))]
	}
	if(max(colorindex)>0){
		#Reads_color = colorRampPalette(brewer.pal(9, "Paired"))(max(colorindex))
		Reads_color = mycolors[1:max(colorindex)]
		#if(length(Reads_color)>7) Reads_color = Reads_color[c(unlist(sapply(1:7,seq,length(Reads_color),by=7)))]
	}
	## Reads_col  ## normal 0 cigar -1 splitmapping(Reads_color,split_color[TRA])
	Reads_col = c()
	for(i in 1:length(colorindex)){
		j = colorindex[i]
		if(TRA){ # TRA
			if(!is.na(match(Reads[i],splitreads))){
				Reads_col = c(Reads_col,split_color[match(Reads[i],splitreads)])
			}else{
				if(j==0) {
					Reads_col=c(Reads_col,normcol) 
			        }else if(j== -1) {Reads_col=c(Reads_col,svcolor) 
				}else if(j==-2) {Reads_col=c(Reads_col,rnacolor)
				}else{ 
					Reads_col=c(Reads_col,Reads_color[j])
				}
			}
		}else{
			if(j==0){ Reads_col=c(Reads_col,normcol) 
			}else if(j== -1) {
				Reads_col=c(Reads_col,svcolor) 
			}else if(j==-2) {
				Reads_col=c(Reads_col,rnacolor) 
			}else{ 
			        Reads_col=c(Reads_col,Reads_color[j])
			}
		}

	}

	## base region
	Data <- data[,-c(1:match('ReadsID',colnames(data)))] # dim(Data)[2] == 0  
	#region = as.numeric(colnames(Data))   ## Absolute value
	region <- start:end
	
	ylab <- ifelse(is.null(sample),'',sample)
	cex.lab = 1.6
	par(cex.lab=cex.lab,xpd=F)
	ylim1 = ifelse(max(data$Readsorder)<5, 5, max(data$Readsorder))
	ylimL = ifelse(end-start<=210,-ylim1-1 - 2,-ylim1-1 - 0.5)
	peak_h = 1/10*ylim1
	plot(1,type="n",xlab="",ylab=ylab, axes=F, xlim=c(min(region)-1,max(region)+1), ylim=c(ylimL,peak_h), xaxs="i", yaxs="i")
	yrow = data$Readsorder
	###rect(xleft = refstart, ybottom = -yrow-1+d, xright = refend, ytop = -yrow-d, col = Reads_col, border = NA) # need triangle_rect
        trirect(xleft = refstart-0.5, ybottom = -yrow-1+d, xright = refend+0.5, ytop = -yrow-d, col = Reads_col, strand = strand,start=start-0.5,end=end+0.5,TRA=TRA)
	##rectarrow(xleft = refstart, ybottom = -yrow-1+d, xright = refend, ytop = -yrow-d, col = Reads_col, strand = strand,start=start,end=end)


	# decorate reads like ins del dup ...
	for(i in 1:nrow(data)){ 
		y1 = -yrow[i] - 1 + d
		y2 = -yrow[i] - d
		DELsplit = FALSE
		INSsplit = FALSE
		DUPsplit = FALSE
		if(i>1){
			if(Reads[i]==Reads[i-1]){
				R1type = strsplit(type[i],"@@")[[1]]
				R2type = strsplit(type[i-1],"@@")[[1]]
				intersectR1R2 = intersect(R1type,R2type)
				if ( any(grepl('DEL--',intersectR1R2)) ) DELsplit = TRUE
				if ( any(grepl('INS--',intersectR1R2)) ) INSsplit = TRUE
				if ( any(grepl('DUP--',intersectR1R2)) ) DUPsplit = TRUE
			}
		}
		
		typei = strsplit(type[i],"@@")[[1]]
		# cigar del ins sc
		for( j in typei){
		    if(grepl('del--',j)){
	                pos_len = as.numeric(gsub('.*:', '', j))
		        pos = gsub('.*_(.*):.*', '\\1', j)
			posl = as.numeric(gsub('-.*', '', pos)) + 1 + 0.5
			posr = as.numeric(gsub('.*-', '', pos))    - 0.5
                        rect(posl, y1, posr, y2, col = "white", border="white")
			segments(posl,-yrow[i]-0.5,posr,-yrow[i]-0.5,col = Reads_col[i], lty=2, lwd=2)
                        text((posl+posr)/2, -yrow[i]-0.5,pos_len, cex=1, col = rgb(118/255,24/255,220/255), font = 2)
		    }
		    if(grepl('ins--',j)){
			pos_len = as.numeric(gsub('.*:', '', j))
		        pos = gsub('.*_(.*):.*', '\\1', j)
			posl = as.numeric(gsub('-.*', '', pos)) 
			posr = as.numeric(gsub('.*-', '', pos)) 
			if(strand[i]=="+"){
			    rect(posl - 0.5, y1, posl + 0.5 + strwidth(pos_len), y2, col = Insert_col,border=Insert_col)
			    text(posl, -yrow[i]-0.5, pos_len, cex=1, adj=0, col="white", font=2)
			}else{
			    rect(posr + 0.5, y1, posr - 0.5 - strwidth(pos_len), y2, col = Insert_col,border=Insert_col)
			    text(posr, -yrow[i]-0.5, pos_len, cex=1, adj=1, col="white", font=2)
			}
		    }
		    if(grepl('sc--',j)){
			pos_len = as.numeric(gsub('.*:', '', j))
		        pos = gsub('.*_(.*):.*', '\\1', j)
			posl = as.numeric(gsub('-.*', '', pos)) 
			posr = as.numeric(gsub('.*-', '', pos))
			segmentsbend(xpos1=posl,y1,y2,strand=strand[i],refstart=refstart[i],refend=refend[i],TRA=TRA)
		    }
		    if(grepl('splic--',j)){
		        pos_len = as.numeric(gsub('.*:', '', j))
		        pos = gsub('.*_(.*):.*', '\\1', j)
			posl = as.numeric(gsub('-.*', '', pos)) + 1 + 0.5
			posr = as.numeric(gsub('.*-', '', pos))    - 0.5
			rect(posl, y1, posr, y2, col = "white", border="white")
			segments(posl,-yrow[i]-0.5,posr,-yrow[i]-0.5,col = Reads_col[i], lty=1, lwd=2)
		    }
		}

			
	         if(end-start<=210&ncol(Data)>0){
	             BaseRegion = colnames(Data)
		     BaseRegion = as.numeric(gsub('_.*','',BaseRegion))
		     #if(grepl('_',BaseRegion[1])) RefBase = sapply(strsplit(BaseRegion,"_"),' [', 2)
		     readbaselist = unlist(Data[i,])
		     basecol <- list()
		     basecol[['A']] = 'green'
		     basecol[['T']] = 'red'
		     basecol[['G']] = 'orange'
		     basecol[['C']] = 'blue'
		     basecol[['N']] = 'black'

		     for(j in 1:ncol(Data)){
			xregionj = BaseRegion[j]
			if((readbaselist[j]==""|readbaselist[j]=="+")) next
			if(readbaselist[j]=="-") { #del
			    rect(xregionj-0.5,y1,xregionj+0.5,y2,col="white",border="white")
			    segments(xregionj - 0.5, -yrow[i] - 0.5, xregionj + 0.5,-yrow[i]-0.5,lwd=2,col = svcolor)
			}else{ 
			    #rect(region[j]-0.5,-yrow-1+d,region[j]+0.5,-yrow-d,col=Reads_col[i],border=Reads_col[i])
			    basej = readbaselist[j]
			    if (nchar(basej)==1){
				text(xregionj,-yrow[i]-0.5,basej, cex=150/(end-start),col=basecol[[basej]], family="mono",font=2)
			    }else{ # indel
				basej1 = strsplit(basej,"")[[1]][1]
			        #text(xregionj,-yrow[i]-0.5,"I",cex=1,col="purple")
			        rect(xregionj-0.5,y1,xregionj+0.5,y2,col="#A020F080",border=NA)  # purple
				text(xregionj,-yrow[i]-0.5,basej1,cex=150/(end-start),col = basecol[[basej1]],family="mono",font=2)
			    }
			 }
			}
		}
		### --------------------------- split mapping  -------------------------------#####
		## DEL split mapping  length-labels
		if(DELsplit){
		    typedeli = intersectR1R2[grepl('DEL--',intersectR1R2)]
		    dellen = gsub('.*:','',typedeli)
		    delpos = gsub('.*_(.*):.*','\\1',typedeli)
	            delx1 = as.numeric(gsub('-.*', '',delpos))+1
	            delx2 = as.numeric(gsub('.*-', '',delpos))
		    segments(delx1,-yrow[i]-0.5,delx2,-yrow[i]-0.5,col=Reads_col[i],lty=2,lwd=2)
		    text((delx1+delx2)/2,-yrow[i]-0.5,dellen,cex=1,col=rgb(118/255,24/255,220/255),font=2)
		}
		## INS split mapping length-labels
		if(INSsplit){
			typeinsi = intersectR1R2[grepl('INS--',intersectR1R2)]
			inslen = gsub('.*:','',typeinsi)
			inspos = gsub('.*_(.*):.*','\\1',typeinsi)
			insx1 = as.numeric(gsub('-.*', '',inspos))
			insx2 = as.numeric(gsub('.*-', '',inspos))

			rect(insx1-0.5,y1,insx1+0.5+strwidth(inslen),y2,col=Insert_col,border=Insert_col)
			text(insx1,-yrow[i]-0.5,inslen,cex=1,adj=0,col="white",font=2)
		}
		## DUP split mapping  length-labels
		if(DUPsplit){
			typedupi = intersectR1R2[grepl('DUP--',intersectR1R2)]
			duplen = gsub('.*:','',typedupi)
			duppos = gsub('.*_(.*):.*','\\1',typedupi)
			dupx1 = as.numeric(gsub('-.*', '',duppos))
			dupx2 = as.numeric(gsub('.*-', '',duppos))
			text((dupx1+dupx2)/2,-yrow[i],duplen, cex=1.5,col="black",font=4)
		}
	}
	## ref base
	par(xpd=T)
	if(end-start<=210&ncol(Data)>0){
	    BaseRegion = colnames(Data)
	    if(grepl('_',BaseRegion[1])) {
		basecol <- list()
	        basecol[['A']] = 'green'
		basecol[['T']] = 'red'
	        basecol[['G']] = 'orange'
		basecol[['C']] = 'blue'
		basecol[['N']] = 'black'
		RefBase = gsub('.*_','',BaseRegion)
		RefBase.col = as.character(sapply(RefBase,function(x)basecol[[x]]))
	        BaseRegion = as.numeric(gsub('_.*','',BaseRegion))
		segments(start-0.5,-ylim1-1-0.5,end+0.5,-ylim1-1-0.5,lwd=2)
	        text(BaseRegion,-ylim1-1-1-0.1,RefBase,cex= 150/(end-start),col = RefBase.col, family="mono",xpd=T,font=2)
		text(par('usr')[1]-1,-ylim1-1-1-0.1,"RefBase",adj=1,font=2,xpd=T,cex=cex.lab)
	    }
	}#else{
	 #   segments(start-0.5,-ylim1-1-0.3,end+0.5,-ylim1-1-0.3,col="grey90")
	#}
	## kurtosis depth or coverage
	y_peakbottom = 0
	if(end-start<=210&ncol(Data)>0){ # d ?
	    peakdata = data0[,-c(1:match('ReadsID',colnames(data0)))]
	    peakdata = ifelse(is.na(peakdata)|peakdata=="-"|peakdata==""|peakdata=="+",0,1)
	    peak = apply(peakdata,2,sum)
	}else{
	    peakdata = t(apply(data.frame(refstart,refend),1,function(x){ma=match(start:end,x[1]:x[2]);ifelse(is.na(ma),0,1)}))
	    peak = apply(peakdata,2,sum)
	    type1 = unlist(strsplit(data$Type,"@@"))
	    type1 = type1[grepl('del--',type1)]
	    if(length(type1)>0){
	        start1 = as.numeric(gsub('del--.*_(.*)-.*','\\1',type1))
	        end1 = as.numeric(gsub('del--.*-(.*):.*','\\1',type1))
	    
	        delpeakdata = t(apply(data.frame(start1,end1),1,function(x){ma=match(start:end,x[1]:x[2]);ifelse(is.na(ma),0,1)}))
                delpeak = apply(delpeakdata,2,sum)
	        peak = peak - delpeak
	        peak[peak<0] = 0
	    }
	  
	}
	type2 = unlist(strsplit(data$Type,"@@"))
	type2 = type2[grepl('splic--',type2)]
	if(length(type2)>0){
	    start1 = as.numeric(gsub('splic--.*_(.*)-.*','\\1',type2))
	    end1 = as.numeric(gsub('splic--.*-(.*):.*','\\1',type2))
	    splicpeakdata = t(apply(data.frame(start1,end1),1,function(x){ma=match(start:end,x[1]:x[2]);ifelse(is.na(ma),0,1)}))
	    splicpeak = apply(splicpeakdata,2,sum)
	    head(splicpeak)
	    peak = peak - splicpeak
	    peak[peak<0] = 0
	}
	maxpeak = max(peak)
	Depth <- paste0("Max\nDepth:", maxpeak)
	peak = peak/maxpeak * peak_h
	segments(start-0.5,y_peakbottom,end+0.5,y_peakbottom,col="grey90",lwd=1) # peak bottom line
	rect(region-0.5,y_peakbottom,region+0.5,y_peakbottom+peak,col=peakcol,border=NA) # col="grey",border="#563624" # #8B7E66

	
	if(left){
		#segments(par('usr')[1]-1/3*(par('usr')[2]-par('usr')[1]),par('usr')[4],par('usr')[1],par('usr')[4],lwd=2,lty=2)
		text(par('usr')[1],y_peakbottom+peak_h/2,Depth,adj=1,cex=cex.lab,font=2)
	}else{
		text(par('usr')[2],y_peakbottom+peak_h/2,Depth,adj=0,cex=cex.lab,font=2)
	}
}

colrgb = function(color,alpha=0.7){
	tmp = col2rgb(color)/255
	RGB = apply(tmp,2,function(x)rgb(x[1],x[2],x[3],alpha))
	return(RGB)
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
	if(cytoband==""){
		cyto = data.frame()
	}else{
		cyto = Readcytoband(file=cytoband,chrom=chrom)
		bei = (max(cyto$end)-min(cyto$start))/(end-start)
		cyto$start = cyto$start/bei +start
		cyto$end = cyto$end/bei +start
		cstart = start/bei + start
		cend = end/bei +start
	}
	
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
        rect(cstart,4-0.2,cend,5+0.2,col="#DC143C",border="#DC143C",lwd=1.5) # red
	par(xpd=T)
	if(left) text(par('usr')[1],4.5,genome,adj=1,cex=2,font=2)
	}


}


GetregionIntersect = function(start1,end1,start2,end2){
        TF = !( end1<start2 | start1>end2)
        return(TF)
}

getannot = function(scriptdir,genome,chrom,start,end,genepred){
	annot = list()
	databasedir = paste0(dirname(scriptdir),"/database")
	if(file.exists(genepred)) {
	    if(grepl('.gz$',genepred)) refgene = gettextf("gzip -dc %s",genepred) else
		    refgene = refgene
	}else{
	    refgenefile = gettextf("%s/%s/%s.genePred.gz",databasedir,genome,genome)
	    if(file.exists(refgenefile)) refgene =  gettextf("gzip -dc %s",refgenefile) else annot[['Refgene']] = data.frame()
	}

	if(genome=="hg19"|genome=="hg38"){
        	segdup = gettextf("gzip -dc %s/%s/%s_genomicSuperDups.bed.gz",databasedir,genome,genome)
        	rmsk = gettextf("gzip -dc %s/%s/%s_rmsk.bed.gz",databasedir,genome,genome)
        	cytoband = paste0(databasedir,"/",genome,"/cytoBand.txt")
	}else{
		#annot[['Refgene']] = data.frame()
		annot[['Rmsk']] = data.frame()
		annot[['Segdup']] = data.frame()
		annot[['cytoband']] = ''
		return(annot)
	}
	Segdup = fread(cmd=segdup,header=F,sep="\t")
	Segdup = as.data.frame(Segdup)
	Segdup[,1] = gsub('chr','',Segdup[,1])
	Segdup = Segdup[Segdup[,1]==chrom,]
	Segdup = Segdup[GetregionIntersect(Segdup[,2],Segdup[,3],start,end),,drop=F]

	Rmsk = fread(cmd=rmsk,sep="\t",header=F)
	Rmsk = as.data.frame(Rmsk)
	Rmsk[,1] = gsub('chr','',Rmsk[,1])
	Rmsk = Rmsk[Rmsk[,1]==chrom,]
	Rmsk = Rmsk[GetregionIntersect(Rmsk[,2],Rmsk[,3],start,end),,drop=F]

	Refgene = fread(cmd=refgene,sep="\t",header=F)
	# transcript,chr,strand,start,end,utr1,utr2,exonnum,exonstart,exonend,0,genename
	Refgene = as.data.frame(Refgene)
	Refgene[,2] = gsub('chr','',Refgene[,2])
	Refgene = Refgene[Refgene[,2]==chrom,]
	Refgene = Refgene[GetregionIntersect(Refgene[,4],Refgene[,4],start,end),,drop=F]
	Refgene = Refgene[order(Refgene$V12,Refgene$V5-Refgene$V4,decreasing=T),,drop=F]
	Refgene = Refgene[!duplicated(Refgene$V12),,drop=F] # get max transcript region

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
	rect(rmskintersect[,2],down,rmskintersect[,3],up,col="#00BFFF",border="#00BFFF")
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

rgb2deeprgb <- function(col,alpha=0.6){
   return(apply(col2rgb(col)/255*alpha,2,function(x)rgb(x[1],x[2],x[3])))
}

drawgene <- function(refgeneintersect,rna=FALSE,TRA=FALSE,y1=8,y2=9,ty1=7.5,ty2=7){ #arrowscol = rgb(120/255,120/255,120/255)
	mygenecol = c("#3386AE", "#8DCD6F", "#889D5A", "#ED5051", "#F58F57", "#FE880F", "#CAB2D6","#A6CEE3")
        if(!TRA) text(par('usr')[1],(y1+y2)/2,"RefGene",adj=1,cex=1.5,font=2,xpd=T)
	xdiff = 1/15*(par('usr')[2]-par('usr')[1])
	refgenen = nrow(refgeneintersect)
        if(refgenen>0){
	    transcriptstart = as.numeric(refgeneintersect[,4])+1
	    transcriptend = as.numeric(refgeneintersect[,5])
	    strand = refgeneintersect[,3]
	    utr1pos = as.numeric(refgeneintersect[,6])
	    utr2pos = as.numeric(refgeneintersect[,7])+1
            exonstart = as.numeric(unlist(strsplit(gsub(',$','',refgeneintersect[,9]),",")))+1
	    exonend = as.numeric(unlist(strsplit(gsub(',$','',refgeneintersect[,10]),",")))
	    exonnum = as.numeric(refgeneintersect[,8])
	    #if(refgenen<=8){
	    exoncol = rep(rep(mygenecol[1:refgenen],length=length(exonnum)),exonnum)
	    genecol = rep(mygenecol[1:refgenen],length=length(exonnum))
	    deepgenecol = rgb2deeprgb(genecol)
	    #}else{
	    #    exoncol = mygenecol[1]
	    #    genecol = rep(mygenecol[1],refgenen)
	    #}
	    par(xpd=F)
	    #segments(transcriptstart,(y1+y2)/2,transcriptend,(y1+y2)/2,col=deepgenecol,lwd=2) # line
	    rect(exonstart,y1,exonend,y2,col=exoncol,border=exoncol)
	    for(j in 1:nrow(refgeneintersect)){
		if(utr1pos[j]!=utr2pos[j]){ # UTR
		    utr1 = rect(utr2pos[j],y1,transcriptend[j],y2,col="white",border="white")
		    utr2 = rect(transcriptstart[j],y1,utr1pos[j],y2,col="white",border="white")
		    # 7:end  start->6
		    rect(utr2pos[j],y1+0.2*(y2-y1),transcriptend[j],y2-0.2*(y2-y1),col=genecol[j],border=genecol[j]) # + 3UTR / - 5UTR
		    rect(transcriptstart[j],y1+0.2*(y2-y1),utr1pos[j],y2-0.2*(y2-y1),col=genecol[j],border=genecol[j]) # + 5UTR / - 3UTR
		}
                ###  arrows
		xleft = ifelse(transcriptstart[j]<start,start,transcriptstart[j])
		xright = ifelse(transcriptend[j]>end,end,transcriptend[j])
	        if(xright-xleft > xdiff) {
		    xa1 = seq(xleft,xright-xdiff,by=xdiff)
		    xa2 = seq(xleft+xdiff,xright,by=xdiff)
		    ya = (y1+y2)/2
		   codej = ifelse(strand[j]=="+",2,1)
		   arrows(xa1, ya, xa2, ya,col=deepgenecol[j],code=codej,lwd=2,length=0.1,angle=20)
	        }
	    }
	    segments(transcriptstart,(y1+y2)/2,transcriptend,(y1+y2)/2,col=deepgenecol,lwd=2) # draw line
		### genename text
	        genejname = refgeneintersect[,12]
		textpos = apply(overlap(transcriptstart,transcriptend,start,end),1,sum)/2
                Text = genejname

		textcex = 0.8
		textwidth = strwidth(Text,cex=textcex)
		for(k in 1:length(textpos)){
		    textkregion = c(textpos[k]-textwidth[k]/2,textpos[k]+textwidth[k]/2)
		    if(k==1) {
	                text(textpos[k],ty1, Text[k],cex=1.2,font=2)
                        textregion = textkregion
                    }else{
                      if(GetregionIntersect(textregion[1],textregion[2],textkregion[1],textkregion[2])){
			    text(textpos[k],ty2, Text[k],cex=1.2,font=2)
			    textregion = textregion
		      }else{
			    text(textpos[k],ty1, Text[k],cex=1.2,font=2)
                            textregion = c(min(textregion[1],textkregion[1]),max(textregion[2],textkregion[2]))
		      }
                      }
                }
	    }
}


drawannot = function(xlim,Annot,start,end,TRA=FALSE,RNA=FALSE){
	plot(-8:10,xlim=xlim,type="n",xlab="",ylab="",xaxs="i",yaxs="i",axes=F)
	## refgene [6.5,10]
	refgeneintersect = Annot[['Refgene']]
	par(xpd=T)
	abline(h=10,col="grey60")
	abline(h=9.9,col="grey60")
        ## gene
	if(RNA){
	    drawgene(refgeneintersect,rna=TRUE,TRA=TRA,y1=2,y2=7,ty1=-2,ty2=-5)
	}else{
	    drawgene(refgeneintersect,rna=FALSE,TRA=FALSE,y1=8,y2=9,ty1=7.5,ty2=7)

	## --- rmsk LINE5.5 SINE2.5 Simple_repeat-0.5 Other-3.5 Segdup-6.5 
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
	if(!TRA)text(par('usr')[1],5.5,"LINE",adj=1,cex=1.5,font=2)
	if(!TRA)text(par('usr')[1],2.5,"SINE",adj=1,cex=1.5,font=2)
	if(!TRA)text(par('usr')[1],-0.5,"Simple_repeat",adj=1,cex=1.5,font=2)
	if(!TRA)text(par('usr')[1],-3.5,"Other",adj=1,cex=1.5,font=2)
	#if(!TRA)text(par('usr')[1],5.5,"Rmsk",adj=1,cex=2)

	##--- segdup
	Segdupintersect = Annot[['Segdup']]
	if(nrow(Segdupintersect)>0){
        	par(xpd=F)
        	rect(Segdupintersect[,2],-7,Segdupintersect[,3],-6,col="#3CB371",border="#3CB371")
	}
	par(xpd=T)
	if(!TRA)text(par('usr')[1],-6.5,"Segdup",adj=1,cex=1.5,font=2)
	}
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



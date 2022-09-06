#!/usr/bin/env Rscript
library('getopt')
library('data.table')

spec = matrix(c(
        'help' ,   'h', 0, "logical",
        'input',   'i', 1, "character",
        'main' ,   'm', 1, "character",
        'samples', 's', 1, "character",
        'outpng',  'o', 1, "character",
	'genome',  'ge',1, "character",
	'main2',   'm2',0, "character",
	'input2',  'i2',0, "character",
	'splitreads', 'ss',0,"character",
	'RNA',     'r',0, "logical",
	'genepred', 'p',0,"character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);
# define usage function
print_usage <- function(spec=NULL){
        cat(getopt(spec, usage=TRUE));
        cat("
Usage example:
1) Rscript SVhawkeye.r -i sample1.bam,sample2.bam -m 1:101212-1191211 -s sample1,sample2  -o /path/outdir/ --genome hg19 
2) Rscript SVhawkeye.r -i sample1.bam.out,sample2.bam.out -m 1:101212-1191211 -s sample1,sample2  -o /path/outdir/ --genome hg19 --input2 sample1.bam.out_2,sample2.bam.out_2 --main2 13:98271-102121

Options:
--help          NULL            get this help;
--input       character       the input Sampleregionpysamout mark "," [forced];
--outpng      character       the outputfile format:png [forced];
--main        character       the region name [forced];
--genome      character       hg38 / hg19 [forced];
--samples     character       sample1,sample2,...[forced];
--main2       character       the region name [when type==TRA/fusion/BND ];
--input2      character       the input Sampleregionpysamout mark "," [when type==TRA/fusion/BND]
--splitreads  character       splitreads mark "," [when type==TRA/fusion/BND]
--RNA         logical         if RNA draw iso-seq reads else draw dnareads(sv...)
--genepred    character       gtfTogenePred from gtf, default:hg19/hg19.genePred.gz 
\n")
        q(status=1);
}
if ( !is.null(opt$help) |is.null(opt$input)) { print_usage(spec) }

scriptdir = dirname(get_Rscript_filename())
sourcefile = paste0(scriptdir,"/source.r")
source(sourcefile)
databasedir = paste0(dirname(scriptdir),"/database")

input = opt$input
main = opt$main
samples = opt$samples
outpng = opt$outpng
if(is.null(opt$genome)) genome = '' else genome = opt$genome
if(is.null(opt$genepred)) {
  if(genome=="hg19"|genome=="hg38") genepred = sprintf("%s/%s/%s.genePred.gz",databasedir,genome,genome) else
     genepred = NA 
}else {
  genepred = opt$genepred
}

if(is.null(opt$main2)|is.null(opt$input2)) TRA = FALSE else TRA = TRUE

options(stringsAsFactors=F)

## deal parameter
input = strsplit(input,",")[[1]]
n = length(input)
samples = strsplit(samples,",")[[1]]
Region = as.numeric(strsplit(strsplit(main,":")[[1]][2],"-")[[1]])
chrom = strsplit(main,":")[[1]][1]
chrom = gsub('chr','',chrom) # need ?
start = Region[1]
end = Region[2]
if(is.null(opt$RNA)) RNA = FALSE else RNA=TRUE
if(TRA){
	splitreads = list()
	splitreadsfile = strsplit(opt$splitreads,",")[[1]]
	for(i in 1:n) splitreads[[i]] = scan(splitreadsfile[i],what="")
}

## annot data: cytoband,ref,rmsk,segdup
Annot = getannot(scriptdir,genome,chrom,start,end,genepred)

##  draw
height = 24+n*6
xlim = c(Region[1]-1,Region[2]+1)

res = ifelse(n>=4,200,300)
if(grepl('.png$',outpng)) png(outpng,height=height,width=40,units="cm",res=res)
if(grepl('.pdf$',outpng)) pdf(outpng,height=min(height/3,200),width=40/3)
if(RNA){
    if(TRA) layout(matrix(1:((n+2)*2),nc=2,byrow=F),heights=c(0.6,rep(2,n),0.4)) else
	    layout(matrix(1:(n+2),nc=1),heights=c(0.6,rep(2,n),0.4))
}else{
    if(genome!="hg19"&genome!="hg38"){ # no annot 
	if(TRA) layout(matrix(1:((n+1)*2),nc=2,byrow=F),heights=c(0.6,rep(2,n))) else
		layout(matrix(1:(n+1),nc=1),heights=c(0.6,rep(2,n)))
    }else{  # annot
	if(TRA) layout(matrix(1:((n+2)*2),nc=2,byrow=F),heights=c(0.6,rep(2,n),1.1)) else
		layout(matrix(1:(n+2),nc=1),heights=c(0.6,rep(2,n),1.1))
}
}
## 1. up  # other species cytobank:blank
par(mar=c(4,8.5,0,3))
drawcytoband(cytoband=Annot[['cytoband']],chrom=chrom,start=start,end=end,xlim=xlim,
    main=main,genome=genome,left=TRUE)
## 2. middle
for(i in 1:n){
	filei = input[i]
	df = fread(filei,header=T,sep="\t",fill=T)
	df = as.data.frame(df)
	#df[is.na(df)]="+"
	data = df[order(df$ReadsID,df$QueryStart),,drop=F]
	sample = gsub('.bam|.sort.bam|.merged|.merge|.fastq','',samples[i])
	par(mar=c(0.5,8.5,2,3))
	if(nrow(data)>0){
		if(is.null(opt$splitreads)) { 
			drawigv(data,sample=sample, TRA = TRA,left=TRUE,nsamples=n,start=start,end=end)  
		}else { # TRA
			drawigv(data,sample=sample, TRA = TRA,splitreads=splitreads[[i]],left=TRUE,nsamples=n,start=start,end=end)
		}
	}else{ # no mapping
		plot(1:10,xlim=xlim,type="n",ylab=sample,xlab="",xaxs="i",yaxs="i",axes=F)
		text(sum(xlim)/2,5,"No Mapping Reads!",cex=2)
	}
}
#####  3. bottom annot  ####
if(genome=="hg19"|genome=="hg38"){
	par(mar=c(0,8.5,0,3))
	drawannot(xlim=xlim,Annot=Annot,start=start,end=end,RNA=RNA)
}
## finish if not TRA
## ----------------------------------------- ##

## TRA right

if(TRA){
	main = opt$main2
	input = opt$input2
	input = strsplit(input,",")[[1]]
	Region = as.numeric(strsplit(strsplit(main,":")[[1]][2],"-")[[1]])
	xlim = c(Region[1]-1,Region[2]+1)
	chrom = strsplit(main,":")[[1]][1]
	chrom = gsub('chr','',chrom)
	start = Region[1]
	end = Region[2]
	Annot = getannot(scriptdir,genome,chrom,start,end,genepred)
	# 1.up
	par(mar=c(4,3,0,8.5))
	drawcytoband(cytoband=Annot[['cytoband']],chrom=chrom,start=start,end=end,xlim=xlim,main=main,genome=genome,left=FALSE)
	# 2.middle
	for(i in 1:n){
        	filei = input[i]
        	df = fread(filei,header=T,sep="\t",fill=T)
        	df = as.data.frame(df)
        	df[is.na(df)]="+"
        	data = df[order(df$ReadsID,df$QueryStart),,drop=F]
        	sample = gsub('.bam|.sort.bam|.merged|.merge|.fastq','',samples[i])
        	par(mar=c(0.5,3,2,8.5))
        	if(nrow(data)>0){
                	drawigv(data,sample="", TRA = TRA,splitreads=splitreads[[i]],left=FALSE,nsamples=n,start=start,end=end)
        	}else{
			print('???')
                	plot(1:10,xlim=xlim,type="n",ylab=sample,xlab="",xaxs="i",yaxs="i",axes=F)
                	text(sum(xlim)/2,5,"No Mapping Reads!",cex=2,xpd=T)
			print(par('usr'))
        	}
	}
	## 3.bottom
	if(genome=="hg19"|genome=="hg38"){
		par(mar=c(0,3,0,8.5))
		drawannot(xlim=xlim,Annot=Annot,TRA=TRA,start=start,end=end,RNA=RNA)
	}
}



## over
dev.off()






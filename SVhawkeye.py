#!/usr/bin/env python2
import os
import re
import sys
import time
import pysam
import argparse
from glob import glob
from multiprocessing import Pool

scriptdir = os.path.dirname(sys.argv[0]) if os.path.dirname(sys.argv[0])!="" else "."
sys.path.append(scriptdir+'/script/')
import sv_vcf

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
	pass
"""
def filter_reads(reads, min_mapq, min_baseq, no_del, no_dup):
	if min_mapq > 0:
		reads = [r for r in reads if r.alignment.mapq >= min_mapq]
	if min_baseq > 0:
		reads = [r for r, q in zip(reads, baseq(reads))
		if q is not None and q >= min_baseq]
	if no_del:
		reads = nodel(reads)
	if no_dup:
		reads = nodup(reads)
	return reads
"""
def pysamdeal(bf,chrom,start,end,min_mapq, min_baseq=0,no_del=False, no_dup=False):
	dout = {}
	#dout.setdefault('ref',{})
	dout.setdefault('reads',{})
	dout.setdefault('strand',{})
	dout.setdefault('ref_start',{})
	dout.setdefault('ref_end',{})
	dout.setdefault('query_start',{})
	dout.setdefault('query_end',{})
	dout.setdefault('reads_length',{})
	dout.setdefault('mapq',{})
	numreads = []
	#for col in bf.pileup(chrom,start,end,truncate=True,stepper="nofilter"): #truncate=True  ignore_overlaps=True flag_filter=0
	for col in bf.pileup(reference=chrom, start=start, end=end, min_base_quality=0, min_mapping_quality=min_mapq,truncate=True):
		reads = col.pileups
		#reads = filter_reads(col.pileups, min_mapq, min_baseq, no_del, no_dup)
		#ref = fafile.fetch(chrom, col.pos, col.pos+1).upper()
		#dout['ref'][col.pos] = ref
		for read in reads:
			#dout.setdefault(read.alignment.query_name+'ref',[]).append(col.pos)
			#if not read.is_del:dout.setdefault(read.alignment.query_name+'read',[]).append(read.query_position)
			if len(numreads)<100 and read.alignment.query_name not in numreads:
				numreads.append(read.alignment.query_name)
			if len(numreads)>=100 and read.alignment.query_name not in numreads:
				continue
			if read.alignment.query_name in numreads:
				readsname = read.alignment.query_name+'__'+str(read.alignment.query_alignment_start)
				dout['reads'].setdefault(readsname,{})
				if read.is_del:
					base = '-'
				elif read.indel >0 :
					base = read.alignment.seq[read.query_position:(read.query_position+read.indel+1)]
				else :
					base = read.alignment.seq[read.query_position]
				dout['reads'][readsname][col.pos]=base
			
				if read.alignment.is_reverse:
					strand = '-'
				else:
					strand = "+"
				dout['strand'][readsname] = strand
				dout['ref_start'][readsname] = read.alignment.reference_start
				dout['ref_end'][readsname] = read.alignment.reference_end
				dout['query_start'][readsname] = read.alignment.query_alignment_start
				dout['query_end'][readsname] = read.alignment.query_alignment_end
				dout['reads_length'][readsname] = read.alignment.query_length
				dout['mapq'][readsname] = read.alignment.mapping_quality
	return dout

def getpysamout(Bamfile,chrom,start,end,outfile,minMapq,TRA):
	bf = pysam.AlignmentFile(Bamfile, 'rb')
	dout = pysamdeal(bf, chrom, start, end, minMapq)
	region = range(start,end)
	out = ('RefStart\tRefEnd\tQueryStart\tQueryEnd\tReadsLen\tmapq\tStrand\tColor\tType\tReads\t'+
		'\t'.join([str(i) for i in range(start,end)]) +'\n')
	preReads = dout['reads'].keys()
	Reads = []  ## no __  readsname     ['bb','bb','cc']
	dupreads = [] ## no __             ['aa','bb']
	dictdup = {} ## dict no __   values has __    {'aa':['aa__32','aa__1234','aa__568'],'bb':['bb_56','bb_345']}
	for i in preReads:
		dictdup.setdefault(i.split('__')[0],[]).append(i)
		if i.split('__')[0] in Reads:
			Reads.append(i.split('__')[0])
			dupreads.append(i.split('__')[0])
		else:
			Reads.append(i.split('__')[0])
	dreadstype = {}  ## has _
	dupreads = list(set(dupreads))
	for i in dupreads:
		for j in dictdup[i]:
			dreadstype.setdefault(j,[])
		dupreadsi = sorted(dictdup[i], key = lambda x:int(x.split('__')[1]) )    ## ['aa__32','aa__568','aa__1234']
		for j in range(len(dupreadsi)-1):    ## 0,1
			k = j + 1
			readj = dupreadsi[j]  ## has __ aa__32
			readk = dupreadsi[k]  ## hsa __ aa__568
			startj = dout['ref_start'][readj]
			endj = dout['ref_end'][readj]
			startk = dout['ref_start'][readk]
			endk = dout['ref_end'][readk]
			## ------------------ INV -------------------- #
			if dout['strand'][readj] != dout['strand'][readk]:
				#if dout['strand'][readj]=="+":
				#	dlen = abs(endj-endk)+1
				#else:
				#	dlen = abs(startj-startk)+1
				dlen = abs(startj-startk)+1 if j==1 else abs(endj-endk)+1
				dreadstype.setdefault(readj,[]).append("INV--"+str(dlen))
				dreadstype.setdefault(readk,[]).append("INV--"+str(dlen))
			else:
				if dout['query_start'][readk] - dout['query_end'][readj]>=100:
					dlen = dout['query_start'][readk] - dout['query_end'][readj]-1
					dreadstype.setdefault(readj,[]).append("INS--"+str(dlen))
					dreadstype.setdefault(readk,[]).append("INS--"+str(dlen))
				if not (startj>endk or endj<startk):
					dlen = min(endk,endj)-max(startk,startj) +1
					if dlen>=40:
						dreadstype.setdefault(readj,[]).append("DUP--"+str(dlen))
						dreadstype.setdefault(readk,[]).append("DUP--"+str(dlen))
				else:
					#if ((startj>endk and dout['strand'][dictdup[i][j]] == "+") or
					#	(endj<startk and dout['strand'][dictdup[i][j]] == "-")):
					if startj>endk:
						dlen = abs(endj-startk)+1
						if dlen>=40:
							dreadstype.setdefault(readj,[]).append("DUP--"+str(dlen))
							dreadstype.setdefault(readk,[]).append("DUP--"+str(dlen))
					else:
						#if dout['strand'][dictdup[i][j]] == "+":
						#	dlen = abs(startk-endj)+1
						#else:
						#	dlen = abs(startj-endk)+1
						dlen = abs(startk-endj)+1
						if dlen>=40:
							dreadstype.setdefault(readj,[]).append("DEL--"+str(dlen))
							dreadstype.setdefault(readk,[]).append("DEL--"+str(dlen))
		for j in dictdup[i]:
			tmp = {}
			for pt in dreadstype[j]:tmp.setdefault(pt[0:5],[]).append(pt[5:])
			dreadstype[j] = '@@'.join([pt+','.join(tmp[pt]) for pt in tmp])

	#Reads = [i.split('_')[0] for i in preReads]
	#dupreads = [i for i in Reads if Reads.count(i)>1]
	dupindex = {}
	n = 0
	#row = 0
	for i in sorted(dout['reads']):
		readsname = i
		Read = i.split('__')[0]   ## read name
		Readslist = [dout['reads'][i][j] if j in dout['reads'][i] else 'N' for j in region]
		BASEchar = ''.join([ j[0] for j in Readslist])   ## read chuan
		Type = ''
		Bnormal = 1
		if Read in dupreads:
			if Read not in dupindex:
				n+=1
				dupindex[Read]=str(n)
			Type += dreadstype[readsname]+'@@'
			Bnormal = 0
		if max([len(k) for k in re.split('A|T|G|C|N',BASEchar)]) >= 40 :
			if Read not in dupindex:
				dupindex[Read]='-1'
			delindex = ( [str((m.start()+m.end())/2)+':'+str(m.end()-m.start())
				for m in re.finditer('----------------------------------------+', BASEchar)])
			Type += "DEL"+'_'+','.join(delindex)+'@@'
			Bnormal = 0
		if max([len(k) for k in Readslist]) >= 40:
			if Read not in dupindex:
				dupindex[Read]='-1'
			INSindex = [str(m+1)+':'+str(len(Readslist[m])) for m in range(len(Readslist)) if len(Readslist[m])>=40]
			Type += 'INS'+'_'+','.join(INSindex)+'@@'
			Bnormal = 0
		if Bnormal:
			if ( (dout['ref_start'][i]>(end-0.25*(end-start)) or 
				dout['ref_end'][i]<(start+0.25*(end-start)) or 
				(dout['ref_end'][i]-dout['ref_start'][i])<0.25*(end-start) )
				and not TRA ) :
				continue
			dupindex[Read]= '0'
			Type = "normal"
		out += ( str(dout['ref_start'][i])+'\t'+str(dout['ref_end'][i])+'\t'+str(dout['query_start'][i])+
			'\t'+str(dout['query_end'][i])+'\t'+str(dout['reads_length'][i])+'\t'+str(dout['mapq'][i])+
			'\t'+dout['strand'][i]+'\t'+dupindex[Read]+'\t'+Type.strip('@@')+'\t'+Read+'\t'+
			'\t'.join([dout['reads'][i][j] if j in dout['reads'][i] else '+' for j in region]) +'\n' )
	bf.close()
	print >>open(outfile,"w"),out.strip()

	return [i.split('__')[0] for i in dout['reads']]


def MakeDir(Dir):
	if not os.path.exists(Dir):
		os.makedirs(Dir)

def isfile(File):
	if not os.path.isfile(File):
		print('Error:%s file is not existed!' %(File))
		sys.exit(1)
	else :
		return 1

def vcf2bed(vcfFile,Bedfile,extendDist):
	fileout = ''
	with open(vcfFile) as io:
		for line in io:
			if '#' in line or line=="\n":
				continue
			line = line.strip()
			sv_record = sv_vcf.sv_vcf_record(line)
			start1 = max(int(sv_record.pos1)-extendDist,0)
			end1 = int(sv_record.pos2)+extendDist
			if sv_record.svtype=="TRA":
				start1 = start1
				start2 = max(int(sv_record.pos2)-extendDist,0)
				end2 = end1
				end1 = int(sv_record.pos1)+extendDist
				bedrow = [sv_record.chrom1,str(start1),str(end1),sv_record.svtype,sv_record.chrom2,str(start2),str(end2)]
				fileout += '\t'.join(bedrow)+'\n'
			elif sv_record.svtype == "INS":
				end1 = int(sv_record.pos1)+extendDist
				bedrow = [sv_record.chrom1,str(start1),str(end1),sv_record.svtype]
				fileout += '\t'.join(bedrow)+'\n'
			else:
				if int(sv_record.pos1) > int(sv_record.pos2) :
					continue
				bedrow = [sv_record.chrom1,str(start1),str(end1),sv_record.svtype]
				fileout += '\t'.join(bedrow)+'\n'
	print >> open(Bedfile,"w"),fileout.strip()

def bed2pysamout(bedlist,bams,pysambeddir,figuredir,minMapq,genome,outfmt):
	line0 = bedlist
	TRA = 1 if len(line0)==7  else 0
	chrom1  = line0[0]
	start0 = line0[1]
	end0   = line0[2]
	start1  = max(int(line0[1]),0)
	end1 = int(line0[2])
	TRA = 0
	if len(line0)==7:
		chrom2 = line0[4]
		start2 = max(int(line0[5]),0)
		end2 = int(line0[6])
		TRA = 1
	if end1 - start1 >500000:
		chrom2 = chrom1
		start2 = end1 - 5000
		end2 = end1
		end1 = start1 + 5000
		TRA = 1
	Sample = []
	bedoutfile = []
	for bam in bams:
		Bamfile = bam
		sample = os.path.basename(bam)
		Sample.append(sample)
		outfile = pysambeddir+'/'+sample+'_'+chrom1+'_' + start0 + '_' + end0
		bedoutfile.append(outfile)
		R1 = getpysamout(Bamfile,chrom1,start1,end1,outfile,minMapq,TRA=TRA)
		main1 = '%s:%s-%s' %(chrom1,str(start1),str(end1))
		samples = ','.join(Sample)
		inputpysamout = ','.join(bedoutfile)
		#if outfmt=="pdf":outpng = '%s/%s_%s_%s.pdf' %(figuredir,chrom1,start0,end0)
		outpng = '%s/%s_%s_%s.pdf' %(figuredir,chrom1,start0,end0) if outfmt=="pdf" else '%s/%s_%s_%s.png' %(figuredir,chrom1,start0,end0)
		if TRA:
			outfile2 = outfile + '_2'
			R2 = getpysamout(Bamfile,chrom2,start2,end2,outfile2,minMapq,TRA=TRA)
			splitmappingreads = [k for k in R1 if k in R2]
			print >>open(outfile+'_ss',"w"),'\n'.join(splitmappingreads)
			main2 = '%s:%s-%s' %(chrom2,str(start2),str(end2))
			inputpysamout2 = ','.join([k+'_2' for k in bedoutfile])
			splitrfile = ','.join([k+'_ss' for k in bedoutfile])
			work = ("%s/script/SVhawkeye.r --input %s --main %s --samples %s --outpng %s --genome %s --main2 %s --input2 %s --splitreads %s"
				%(os.path.realpath(os.path.dirname(sys.argv[0])), inputpysamout, main1, samples, outpng,genome, main2, inputpysamout2, splitrfile) )
		else:
			work = ("%s/script/SVhawkeye.r --input %s --main %s --samples %s --outpng %s --genome %s"
				%(os.path.realpath(os.path.dirname(sys.argv[0])), inputpysamout, main1, samples, outpng, genome) )
	return work


def runpysamout(bedfile,bams,pysamfile,pysambeddir,figuredir,minMapq,genome,outfmt,threads):
	pool = Pool(processes=threads)
	results = []
	bedf = open(bedfile,"r").readlines()
	for line in bedf:
		if line=="\n" or '#' in line:continue
		line0 = line.strip().split('\t')
		if int(line0[1]) > int(line0[2]):continue
		bedlist = line0
		results.append(pool.apply_async(bed2pysamout,(bedlist,bams,pysambeddir,figuredir,minMapq,genome,outfmt)))
	pool.close()
	pool.join()
	#work = [r.get() for r in enumerate(results)]
	work = [r.get() for r in results]
	print >> open(pysamfile,"w"),'\n'.join(work)

def _TF(x):
	y = True if x==True or x=="True" or x=="T" or x=="1" else False
	return y

def main(args):
	bams = args[0].bams
	bams = bams.split(',')
	for i in bams: isfile(i)
	bedvcf = args[0].bedvcf
	isfile(bedvcf)
	bedvcf = os.path.realpath(bedvcf)
	genome = args[0].genome
	outdir = args[0].outdir
	thread = max(int(args[0].thread)+1,1)
	minMapq = int(args[0].quanlty)
	script = outdir+'/script'
	MakeDir(outdir)
	MakeDir(script)
	outdir = os.path.realpath(outdir)

	if args[0].infmt == "vcf":
		bed = outdir+'/input.bed'
		dcoordinate = vcf2bed(bedvcf,bed,int(args[0].extend))
	else:
		bed = bedvcf
	pysambeddir = outdir+'/bedpysamout'
	figuredir = outdir+'/figure'
	MakeDir(pysambeddir)
	pysamscript = script+'/Rigvfrompysam.sh'
	outfmt = args[0].outfmt



	if _TF(args[0].runpysam):
		runpysamout(bed,bams,pysamscript,pysambeddir,figuredir,minMapq,genome,outfmt,thread)
	if _TF(args[0].Drawpicture) and isfile(pysamscript):
		MakeDir(figuredir)
		os.system("perl %s/script/multi_cpu.pl %s %s" %(os.path.dirname(sys.argv[0]),str(thread), pysamscript))
	if _TF(args[0].createvcf) and args[0].infmt == "vcf" : 
		import sv_genotype
		newvcf = sv_genotype.vcf2Genotype(bedvcf,int(args[0].extend),thread)
		newvcf.createvcf(bams,outdir)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		formatter_class=HelpFormatter,
		description='''
		python2 Rpysamigv.py -i M609.0.bam,M609.1.bam,M609.2.bam -g hg19 -b igv.bed -o test
		python2 Rpysamigv_v2.py -i C553-T.bam,C553-N.bam -b test1.vcf --format vcf -d 1000 -o test1 -g hg19
		
		testdir : /export/home/xiaoyh/Pipeline/IGVscreenshot/test
	''')
	parser.add_argument('-i', '--bams', metavar='FILE', required=True,help='set the input bam file. mark=","')
	parser.add_argument('-g', '--genome', metavar='hg19/hg38', default="hg19",help='set reference,support:hg19/hg38')
	parser.add_argument('-b', '--bedvcf', metavar='FILE', required=True,help='set the input bed or vcf file.')
	parser.add_argument('-o', '--outdir', metavar='Dir', required=True,help='set output dirname.')
	parser.add_argument('-t', '--thread', metavar=int,default=0, help='Number of additional threads to use [0]')
	parser.add_argument('-q', '--quanlty', metavar='num',default=20, help='set reads mapping quanlty for filter,default:Q20 ')
	parser.add_argument('-d', '--extend', metavar=int,default=1000, help='set region extend length,if format==vcf,default:1000bp')
	parser.add_argument('-f', '--infmt', metavar='vcf/bed',default='bed', help='set input format:vcf or bed;default:bed')
	parser.add_argument('-p', '--runpysam', metavar='True/False',default=True, help='does run pysam for reads out?')
	parser.add_argument('-D', '--Drawpicture', metavar='True/False',default=True, help='does draw IGV picture')
	parser.add_argument('-c', '--createvcf', metavar='True/False',default=True, help='does create new vcf? recall sv')
	parser.add_argument('-fo', '--outfmt', metavar='png/pdf',default='png', help='set out picture format,default:png')
	

	args = parser.parse_known_args()
	main(args)

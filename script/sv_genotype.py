import os
import re
import sv_vcf
from multiprocessing import Pool


def samplegenotype(dsvinfo,sample):
	chrom1 = dsvinfo[0]
	start0 = dsvinfo[1]
	end0 = dsvinfo[2]
	Type = dsvinfo[3]
	length = dsvinfo[-1]
	file1 = sample+'_'+chrom1+'_' + str(start0) + '_' + str(end0)
	if not os.path.isfile(file1):return '.:.:.:.'
	if len(dsvinfo) == 8 or end0-start0>500000:
		files = file1+'_ss'     ## duan dian support reads
		DVtmp = open(files,"r").readlines()
		DVtmp = list(set([i for i in DVtmp if i!="\n"]))
		DAf1 = open(file1,"r").readlines()
		DA1 =  list(set([i.strip().split('\t')[9] for i in DAf1 if 'Reads' not in i]))
		DAf2 = open(file1+'_2',"r").readlines()
		DA2 =  list(set([i.strip().split('\t')[9] for i in DAf2 if 'Reads' not in i]))
		DA = DA1 if len(DA1)<len(DA2) else DA2
		DV = DVtmp
		DO = []
	else:
		f=open(file1,"r").readlines()
		reads = []
		target = []
		othertaget = []
		for row in range(len(f)):
			if row == 0 : continue
			line = f[row].strip('\n').split('\t')
			reads.append(line[9])
			if Type == "DUP" or Type=="INV":
				if Type in line[8]:
					line[8] = re.sub('@@.*','',re.sub('.*('+Type+'.*)','\\1',line[8]))
					info = line[8].split('--')[1]
					Bool = 0
					for j in info.split(','):
						if int(j) >=length*0.7 and int(j) <= length*1.3:
							Bool = 1
						else:
							othertaget.append(line[9])
					if Bool :
						target.append(line[9])
			if Type == "DEL":
				if 'DEL_' in line[8]:
					line[8] = re.sub('@@.*','',re.sub('.*(DEL_.*)','\\1',line[8]))
					info = line[8].split('_')[1]
					Bool = 0
					for j in info.split(','):
						if int(j.split(':')[1]) >=length*0.7 and int(j.split(':')[1]) <= length*1.3:
							Bool =1
						else:
							othertaget.append(line[9])
					if Bool:
						target.append(line[9])
				if 'DEL--' in line[8]:
					line[8] = re.sub('@@.*','',re.sub('.*(DEL--.*)','\\1',line[8]))
					info = line[8].split('--')[1]
					Bool = 0
					for j in info.split(','):
						if int(j) >=length*0.7 and int(j) <= length*1.3:
							Bool = 1
						else:
							othertaget.append(line[9])
					if Bool :
						target.append(line[9])
			if Type == "INS":
				if 'INS_' in line[8]:
					line[8] = re.sub('@@.*','',re.sub('.*(INS_.*)','\\1',line[8]))
					info = line[8].split('_')[1]
					Bool = 0
					for j in info.split(','):
						if int(j.split(':')[1]) >=length*0.7 and int(j.split(':')[1]) <= length*1.3:
							Bool =1
						else:
							othertaget.append(line[9])
					if Bool:
						target.append(line[9])
				if 'INS--' in line[8]:
					line[8] = re.sub('@@.*','',re.sub('.*(INS--.*)','\\1',line[8]))
					info = line[8].split('--')[1]
					Bool = 0
					for j in info.split(','):
						if int(j) >=length*0.7 and int(j) <= length*1.3:
							Bool = 1
						else:
							othertaget.append(line[9])
					if Bool:
						target.append(line[9])
				if 'DUP' in line[8]:
					line[8] = re.sub('@@.*','',re.sub('.*(DUP.*)','\\1',line[8]))
					info = line[8].split('--')[1]
					Bool = 0
					for j in info.split(','):
						if int(j) >=length*0.7 and int(j) <= length*1.3:
							Bool = 1
						else:
							othertaget.append(line[9])
					if Bool:
						target.append(line[9])
		DV = target
		DA = reads
		DO = othertaget
		DO = list(set(othertaget)-set(DV))
	DV = len(list(set(DV)))
	DA = len(list(set(DA)))
	DO = len(list(set(DO)))
	if DA == 0:
		genotype = '.:.:.:.'
	else:
		AF = DV / float(DA)
		if AF<0.1 or DV<2 :
			genotype = '0/0' if 'X' not in chrom1 and 'Y' not in chrom1 else '0'
			genotype = genotype+':'+str(DV)+':'+str(DA)+":"+str(DO)
		elif AF<0.8:
			genotype = '0/1' if 'X' not in chrom1 and 'Y' not in chrom1 else '1'
			genotype = genotype+':'+str(DV)+':'+str(DA)+":"+str(DO)
		else:
			genotype = '1/1' if 'X' not in chrom1 and 'Y' not in chrom1 else '1'
			genotype = genotype+':'+str(DV)+':'+str(DA)+":"+str(DO)
	return genotype

def Samplesgenotype(vout,dsvinfo,samples):
	for i in range(len(samples)):
		Sample = samples[i]
		Sgenotype = samplegenotype(dsvinfo,Sample)
		vout.append(Sgenotype)
	return '\t'.join(vout)


class vcf2Genotype(object):
	def __init__(self,vcfFile,extendDist,threads):
		self.vcfFile = vcfFile
		self.extendDist = extendDist
		self.threads = threads
		self.Samplesgenotype = Samplesgenotype
	
	def vcf2bed(self):
		extendDist = self.extendDist
		dcoordinate = dict()
		with open(self.vcfFile,"r") as io:
			for line in io:
				if '#' in line or line=="\n":continue
				line = line.strip()
				sv_record = sv_vcf.sv_vcf_record(line)
				start1 = max(int(sv_record.pos1)-extendDist,0)
				end1 = int(sv_record.pos1)+extendDist if sv_record.svtype == "INS" else int(sv_record.pos2)+extendDist
				svid = sv_record.id
				svlen = sv_record.svlen
				if start1 > end1 and sv_record.svtype!="TRA":continue
				
				if sv_record.svtype=="TRA":
					start1 = start1
					start2 = max(int(sv_record.pos2)-extendDist,0)
					end2 = end1
					end1 = int(sv_record.pos1)+extendDist
					bedrow = [sv_record.chrom1,start1,end1,sv_record.svtype,sv_record.chrom2,start2,end2]
				else:
					bedrow = [sv_record.chrom1,start1,end1,sv_record.svtype]
				dcoordinate[svid] = bedrow+[svlen]
		return dcoordinate

	def createvcf(self,bams,outdir):
		newvcf = outdir+'/new.'+os.path.basename(self.vcfFile)
		out = ''
		samples = [os.path.basename(sample) for sample in bams]
		dsv = self.vcf2bed()
		os.chdir(outdir+'/bedpysamout/')
		pool = Pool(self.threads)
		results = []
		with open(self.vcfFile) as io:
			for line in io:
				line0 = line.strip().split('\t')
				if '#' in line:continue
				if out=="":
					out += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
					out += '##FORMAT=<ID=DA,Number=1,Type=Integer,Description="# number of high-quality reads(depth)">\n'
					out += '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# number of high-quality variant reads">\n'
					out += '##FORMAT=<ID=DO,Number=1,Type=Integer,Description="# number of high-quality other variant reads">\n'
					out += ('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + 
						'\t'.join([re.sub('.bam$','',sample) for sample in samples]) +'\n')
				vout = line0[0:8]+['GT:DV:DA:DO']
				SVID = line0[2]
				dsvinfo = dsv[SVID]
				results.append(pool.apply_async(self.Samplesgenotype,(vout,dsvinfo,samples)))
		io.close()
		pool.close()
		pool.join()
		out += '\n'.join([r.get() for r in results])
		print >>open(newvcf,"w"),out.strip()
				


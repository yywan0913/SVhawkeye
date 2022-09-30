import os
import re
import sv_vcf
from base import *
from multiprocessing import Pool


def samplegenotype(dsvinfo,sample,extendDist):
    chrom1 = dsvinfo[0]
    start0 = dsvinfo[1]
    end0 = dsvinfo[2]
    Type = dsvinfo[3]
    length = dsvinfo[-1]
    downlenbei = 0.7
    uplenbei = 1.3
    file1 = sample+'_'+chrom1+'_' + str(max(start0-int(extendDist),0)) + '_' + str(end0+int(extendDist))
    if not os.path.isfile(file1):return '.:.:.:.'
    if len(dsvinfo) == 8 : 
        files = file1+'_ss'     ## TRA support reads
        DVtmp = open(files,"r").readlines()
        DVtmp = list(set([i for i in DVtmp if i!="\n"]))
        DAf1 = open(file1,"r").readlines()
        DA1pos = [int(float(i.strip().split('\t')[11])) for i in DAf1[1:]]
        DA1 =  0 if DA1pos==[] else max(DA1pos)
        DAf2 = open(file1+'_2',"r").readlines()
        DA2pos = [int(float(i.strip().split('\t')[11])) for i in DAf2[1:]]
        DA2 =  0 if DA2pos==[] else max(DA2pos)
        DA = DA1 if DA1<DA2 else DA2
        DV = len(DVtmp)
        DO = 0
    else:
        f=open(file1,"r").readlines()
        reads = []
        target = []
        target_pos = []
        othertaget = []
        reads_pos = []
        for row in range(len(f)):
            if row == 0 : continue
            line = f[row].strip('\n').split('\t')
            readstypei = line[10].upper()
            readsi = line[12]
            reads.append(readsi)
            reads_pos.append(int(float(line[11])))

            if Type == "DUP" or Type=="INV" or Type=="DEL" or Type=="INS":
                if Type=="INS":readstypei = readstypei.replace('DUP','INS')
                if Type=="DUP":readstypei = readstypei.replace('INS','DUP')
                readstypei = readstypei.split('@@')
                typeinfo = [ k for k in readstypei if Type in k]
                for j in typeinfo: # j: DEL--chr6_start-end:len
                    svlengthj = j.split(':')[-1]
                    svstartj = int(j.split('_')[1].split('-')[0])
                    svendj = int(j.split('-')[-1].split(':')[0])
                    if Type=="DEL" or Type=="INV":
                        if intersectregionlen(start0,end0,svstartj,svendj) <0.5*max(end0-start0,svendj-svstartj):continue
                    else:
                        if not intersectregion(start0,end0,svstartj,svendj):continue
                    if Type == "INV":
                        downlenbei = 0.5
                        uplenbei = 2
                    if int(svlengthj) >=length*downlenbei and int(svlengthj) <= length*uplenbei:
                        target.append(readsi)
                        target_pos.append(int(float(line[11])))
                    else:
                        othertaget.append(int(float(line[11])))
                if Type == "DUP" or Type == "INS":
                    scinfo = [k for k in readstypei if 'SC' in k]
                    for j in scinfo:
                        svlengthj = j.split(':')[-1]
                        svstartj = int(j.split('_')[1].split('-')[0])
                        if intersectregion(start0,end0,svstartj-50,svstartj+50):
                            if int(svlengthj) >=length*0.3 and int(svlengthj) <= length*1.3:
                                target.append(readsi)
                                target_pos.append(int(float(line[11])))


        DV = len(list(set(target)))
        DA = max(reads_pos)
        DO = len(list(set(othertaget)-set(target_pos)))

    if DA == 0:
        genotype = '.:.:.:.'
    else:
        AF = DV / float(DA)
        if AF<0.1 or DV<2 :
            genotype = '0/0'
            #genotype = '0/0' if 'X' not in chrom1 and 'Y' not in chrom1 else '0'
            genotype = genotype+':'+str(DV)+':'+str(DA)+":"+str(DO)
        elif AF<0.8:
            genotype = '0/1'
            #genotype = '0/1' if 'X' not in chrom1 and 'Y' not in chrom1 else '1'
            genotype = genotype+':'+str(DV)+':'+str(DA)+":"+str(DO)
        else:
            genotype = '1/1'
            #genotype = '1/1' if 'X' not in chrom1 and 'Y' not in chrom1 else '1'
            genotype = genotype+':'+str(DV)+':'+str(DA)+":"+str(DO)	
        #genotype = ';'.join(DV_reads)+'#'+genotype+'#'+';'.join(DO_reads)
    ##print (file1+'\t'+Type+'\t'+str(length)+'\t'+genotype+'\t'+';'.join(list(set(DV_reads)))+'\t'+';'.join(list(set(DA_reads))) )
    return genotype

def Samplesgenotype(vout,dsvinfo,samples,extendDist):
    for i in range(len(samples)):
        Sample = samples[i]
        Sgenotype = samplegenotype(dsvinfo,Sample,extendDist)
        vout.append(Sgenotype)
    return '\t'.join(vout)


class vcf2Genotype(object):
    def __init__(self,vcfFile,extendDist,threads):
        self.vcfFile = vcfFile
        self.extendDist = extendDist
        self.threads = threads
        self.Samplesgenotype = Samplesgenotype
	
    def vcf2svinfo(self):
        #extendDist = self.extendDist
        dcoordinate = dict()
        with open(self.vcfFile,"r") as io:
            for line in io:
                if '#' in line or line=="\n":continue
                line = line.strip()
                sv_record = sv_vcf.sv_vcf_record(line)
                start1 = max(int(sv_record.pos1),0)
                end1 = int(sv_record.pos1) if sv_record.svtype == "INS" else int(sv_record.pos2)
                svid = sv_record.id
                svlen = sv_record.svlen
                if start1 > end1 and sv_record.svtype!="TRA":continue
				
                if sv_record.svtype=="TRA":
                    start1 = start1
                    start2 = max(int(sv_record.pos2),0)
                    end2 = end1
                    end1 = int(sv_record.pos1)
                    bedrow = [sv_record.chrom1,start1,end1,sv_record.svtype,sv_record.chrom2,start2,end2]
                else:
                    bedrow = [sv_record.chrom1,start1,end1,sv_record.svtype]
                dcoordinate[svid] = bedrow+[svlen]
        return dcoordinate

    def createvcf(self,bams,outdir):
        newvcf = outdir+'/new.'+os.path.basename(self.vcfFile)
        out = ''
        head = ''
        samples = [os.path.basename(sample) for sample in bams]
        dsv = self.vcf2svinfo()
        os.chdir(outdir+'/bedpysamout/')
        pool = Pool(self.threads)
        results = []
        with open(self.vcfFile) as io:
            for line in io:
                if '#' in line:
                    head+=line
                    continue
                line = line.strip().split('\t')
                if out=="":
                    out = head
                    out += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                    out += '##FORMAT=<ID=DA,Number=1,Type=Integer,Description="# number of high-quality reads(depth)">\n'
                    out += '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# number of high-quality variant reads">\n'
                    out += '##FORMAT=<ID=DO,Number=1,Type=Integer,Description="# number of high-quality other same variant reads">\n'
                    out += ('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + 
                        '\t'.join([re.sub('.bam$','',sample) for sample in samples]) +'\n')
                vout = line[0:8]+['GT:DV:DA:DO']
                SVID = line[2] 
                dsvinfo = dsv[SVID]
                results.append(pool.apply_async(self.Samplesgenotype,(vout,dsvinfo,samples,self.extendDist)))
        io.close()
        pool.close()
        pool.join()
        out += '\n'.join([r.get() for r in results])
        writefile(out.strip(),newvcf)
				


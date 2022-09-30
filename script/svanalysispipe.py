import os,sys
SCRIPTDIR, SCRIPTNAME=os.path.split(os.path.abspath(sys.argv[0]))
mainscript = os.path.realpath(SCRIPTDIR+'/script')
sys.path.append(mainscript)
import pysam
from optparse import OptionParser
from base import *
from pysam2svprequel import *
from sv_vcf import *
from multiprocessing import Pool

class SVToView(object):
    def __init__(self,):
        self.mainscript = mainscript
    
    def GetOpt(self,genotyping=0):
        parser = OptionParser();
        parser.add_option('-i', '--bams', action="store", dest="bams",metavar='FILE', help='set the input bam file. mark=","; [required:True]')
        parser.add_option('-g', '--genome', action="store",dest="genome",default="hg19", help='set reference,support:hg19/hg38;while,other genome also can draw but no annotation; [default:hg19]')
        parser.add_option('-b', '--bedvcf', action="store", dest="bedvcf",metavar='FILE', help='set the input bed or vcf file; [required=True]')
        parser.add_option('-r', '--reffa', action="store",dest="reffa",metavar='FILE', default=None, help='set the reference fasta file of inputbam, when region<210bp and which can dispaly ref base; [default:None]')
        parser.add_option('-o', '--outdir', action="store",dest="outdir",metavar='Dir', default="./", help='set output dirname; [default: ./]')
        parser.add_option('-t', '--thread', action="store",dest="thread",metavar=int, default=0, help='Number of additional threads to use ; [default:0]')
        parser.add_option('-q', '--quanlty', action="store", dest="quanlty",metavar='num',default=20, help='set reads mapping quanlty for filter; [default: 20 (means Q20)]')
        parser.add_option('-I', '--identity', action="store", dest="identity",metavar=float,default=0.6, help='set min identity of mapping reads for filter; [default: 0.6]')
        parser.add_option('-d', '--extend', action="store",dest="extend",metavar=int, default=1000, help='set region extend length; [default:1000bp]')
        parser.add_option('-f', '--infmt', action="store", dest="infmt",metavar='vcf/bed',default='bed', help='set input format:vcf or bed; [default:bed]')
        parser.add_option('-F', '--outfmt', action="store", dest="outfmt",metavar='png/pdf',default='png', help='set out picture format; [default:png]')
        parser.add_option('-l', '--sv_min_length', action="store", dest="sv_min_length",metavar=int, default=50, help='Minimum length of SV to be reported; [default:50]')
        (options, args) = parser.parse_args()
        if not options.bams or not options.bedvcf:
            parser.print_help();
            sys.exit(0)
        self.infmt = options.infmt
        if genotyping and (self.infmt != "vcf" and options.bedvcf[-4:]!=".vcf"):
            print("Warning:\n    sv_genotyping: sv recall must input vcf file!\n")
            parser.print_help();
            sys.exit(0)
        self.bams = options.bams.split(',')
        for i in self.bams:isfile(i)
        self.reffa = options.reffa
        self.bedvcf = os.path.realpath(options.bedvcf)
        self.genome = options.genome
        self.outdir = os.path.realpath(options.outdir)
        self.thread = max(int(options.thread),1)
        self.minMapq = int(options.quanlty)
        self.SV_min_length = int(options.sv_min_length)
        self.infmt = options.infmt
        self.outfmt = options.outfmt
        self.script = self.outdir+'/script'
        self.extend = int(options.extend)
        MakeDir(self.outdir)
        MakeDir(self.script)
        self.inbed = self.outdir+'/input.bed'
        if self.infmt == "vcf":
            self.dcoordinate = self.vcf2bed()
        else:
            self.dcoordinate = self.Bed2bed()
        self.pysambeddir = self.outdir+'/bedpysamout'
        self.figuredir = self.outdir+'/figure'
        MakeDir(self.pysambeddir)
        self.pysamscript = self.script+'/Rigvfrompysam.sh'

    def pysamToinfo(self,Bamfile,chrom,start,end,outfile):
        bf = pysam.AlignmentFile(Bamfile, 'rb')
        pyref = None if self.reffa==None else pysam.FastaFile(self.reffa)
        allreads = bf.fetch(chrom,start,end) ## add min_map
        dsv={}
        n = 0
        for read in allreads:
            n+=1
            if n>2000:continue
            if read.mapq <self.minMapq :continue
            dsv = resolve_read(read,start,end).cigar2sv(dsv) ## add min_map self.minMapq

        dsv = resolve_split_reads(dsv).get_sp_sv(self.SV_min_length)

        if end-start<=210:
            Reads = dsv2info(dsv,bf,chrom,start,end,outfile,pyref,min_mapq=self.minMapq,loci=True)
        else:
            Reads = dsv2info(dsv,bf,chrom,start,end,outfile,pyref,min_mapq=self.minMapq,loci=False)
        return(Reads)

    def bed2pysamout(self,bedlist): # single row
        outfile_prefix = self.dcoordinate['\t'.join(bedlist)] # chr_start_end.e1000
        TRA = 1 if len(bedlist)==7  else 0
        chrom1  = bedlist[0]
        start0 = bedlist[1] # str
        end0   = bedlist[2] 
        start1  = max(int(bedlist[1]),0) # int
        end1 = int(bedlist[2])
        if TRA:
            chrom2 = bedlist[4]
            start2 = max(int(bedlist[5]),0)
            end2 = int(bedlist[6])
        Sample = []
        bedoutfile = []
        main1 = '%s:%s-%s' %(chrom1,str(start1),str(end1))
        for bam in self.bams:
            Bamfile = bam
            sample = os.path.basename(Bamfile)
            Sample.append(sample)
            outfile = self.pysambeddir+'/'+sample+'_'+chrom1+'_' + start0 + '_' + end0
            bedoutfile.append(outfile)
            R1 = self.pysamToinfo(Bamfile,chrom1,start1,end1,outfile)
            if TRA:
                outfile2 = outfile + '_2'
                R2 = self.pysamToinfo(Bamfile,chrom2,start2,end2,outfile2)
                splitmappingreads = [k for k in R1 if k in R2] # R1 R2
                writefile('\n'.join(splitmappingreads),outfile+'_ss')
        samples = ','.join(Sample)
        inputpysamout = ','.join(bedoutfile)
        outpng = '%s/%s.pdf' %(self.figuredir,outfile_prefix) if self.outfmt=="pdf" else '%s/%s.png' %(self.figuredir,outfile_prefix)
        if TRA:
            main2 = '%s:%s-%s' %(chrom2,str(start2),str(end2))
            inputpysamout2 = ','.join([k+'_2' for k in bedoutfile])
            splitrfile = ','.join([k+'_ss' for k in bedoutfile])
            work = ("Rscript %s/SVhawkeye.r --input %s --main %s --samples %s --outpng %s --genome %s --main2 %s --input2 %s --splitreads %s"
                    %(self.mainscript, inputpysamout, main1, samples, outpng,self.genome, main2, inputpysamout2, splitrfile) )
        else:
            work = ("Rscript %s/SVhawkeye.r --input %s --main %s --samples %s --outpng %s --genome %s"
                    %(self.mainscript, inputpysamout, main1, samples, outpng, self.genome) )
        return work # for draw  -->runpysamout

    def runpysamout(self):
        pool = Pool(processes=self.thread)
        results = []
        for line in self.dcoordinate:
            if line=="\n" or '#' in line:continue
            line = line.strip().split('\t')
            if int(line[1]) > int(line[2]):continue
            bedlist = line
            results.append(pool.apply_async(self.bed2pysamout,(bedlist,)))

        pool.close()
        pool.join()
        work = [r.get() for r in results]
        writefile('\n'.join(work),self.pysamscript)


    def vcf2bed(self):
        dinput = {}
        fileout = ''
        with open(self.bedvcf) as io:
            for line in io:
                if '#' in line or line=="\n":
                    continue
                line = line.strip()
                sv_record = sv_vcf_record(line)  ## from sv_vcf
                start1 = max(int(sv_record.pos1) - self.extend,0)
                end1 = int(sv_record.pos2)+ self.extend
                if sv_record.svtype=="TRA":
                    start1 = start1
                    start2 = max(int(sv_record.pos2)-self.extend,0)
                    end2 = end1
                    end1 = int(sv_record.pos1)+self.extend
                    bedrow = [sv_record.chrom1,str(start1),str(end1),sv_record.svtype,sv_record.chrom2,str(start2),str(end2)]
                    fileout += '\t'.join(bedrow)+'\n'
                    dinput['\t'.join(bedrow)] = sv_record.chrom1+'_'+str(sv_record.pos1)+'_'+str(sv_record.pos2)+'.e'+str(self.extend)
                elif sv_record.svtype == "INS":
                    end1 = int(sv_record.pos1)+self.extend
                    bedrow = [sv_record.chrom1,str(start1),str(end1),sv_record.svtype]
                    fileout += '\t'.join(bedrow)+'\n'
                    dinput['\t'.join(bedrow)] = sv_record.chrom1+'_'+str(sv_record.pos1)+'_'+str(sv_record.pos1)+'.e'+str(self.extend)
                else:
                    if int(sv_record.pos1) > int(sv_record.pos2):
                        continue
                    bedrow = [sv_record.chrom1,str(start1),str(end1),sv_record.svtype]
                    fileout += '\t'.join(bedrow)+'\n'
                    dinput['\t'.join(bedrow)] = sv_record.chrom1+'_'+str(sv_record.pos1)+'_'+str(sv_record.pos2)+'.e'+str(self.extend)
        writefile(fileout.strip(),self.inbed)
        return dinput

    def Bed2bed(self): 
        dinput = {}
        fileout = ''
        f=open(self.bedvcf,"r").readlines()
        for line in f:
            line = line.strip().split('\t')
            if '#' in line[0] or len(line)<3:
                continue
            if len(line)<=5: # chr start end --> chr start-1000 end+1000 : char_start_end.e1000
                start1 = str(max(int(line[1])-self.extend,0))
                end1 = str(int(line[2])+self.extend)
                fileout += '\t'.join([line[0],start1,end1])+'\n'
                dinput['\t'.join([line[0],start1,end1])] = line[0]+'_'+line[1]+'_'+line[2]+'.e'+str(self.extend)
            else: # chr1 start1-1000 end1+1000 TRA chr2-1000 start2+1000 end2 --> ..:chr1 start1 end2.e1000
                chrom1 = line[0]
                start1 = str(max(int(line[1])-self.extend,0))
                end1 = str(int(line[2])+self.extend)
                chrom2 = line[4]
                start2 = str(max(int(line[5])-self.extend,0))
                end2 = str(int(line[6])+self.extend)
                fileout += '\t'.join([line[0],start1,end1,'TRA',chrom2,start2,end2])+'\n'
                dinput['\t'.join([line[0],start1,end1,'TRA',chrom2,start2,end2])] = line[0]+'_'+line[1]+'_'+line[4]+'.e'+str(self.extend)
        writefile(fileout.strip(),self.inbed)
        return dinput

    def runsvbrowser(self):
        self.runpysamout()
        MakeDir(self.figuredir)
        pool = Pool(processes=self.thread)
        fio = open(self.pysamscript,"r").readlines()
        for line in fio:
            pool.apply_async(Runcmd,(line.strip('\n'),))
        pool.close()
        pool.join()
    def runsvgenotyping(self):
        self.runpysamout()
        import sv_genotype
        newvcf = sv_genotype.vcf2Genotype(self.bedvcf,int(self.extend),self.thread)
        newvcf.createvcf(self.bams,self.outdir)



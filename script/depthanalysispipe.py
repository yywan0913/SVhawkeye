import os,sys
mainscript = os.path.dirname(os.path.realpath(sys.argv[0]))+'/script'
sys.path.append(mainscript)
import pysam
import math
from optparse import OptionParser
from base import *
from multiprocessing import Pool

class regionToView(object):
    def __init__(self,):
        self.drawdepth = mainscript+'/depth2distribution.r'
    
    def GetOpt(self,genotyping=0):
        parser = OptionParser();
        parser.add_option('-i', '--bams', action="store", dest="bams",metavar='FILE', help='set the input bam file. mark=","; [required:True]')
        parser.add_option('-r', '--region', action="store",dest="region",metavar='chr:from-to', default=None, help='set region like chr1:1000-3000')
        parser.add_option('-o', '--outdir', action="store",dest="outdir",metavar='Dir', default="./", help='set output dirname; [default: ./]')
        parser.add_option('-t', '--thread', action="store",dest="thread",metavar=int, default=0, help='Number of additional threads to use ; [default:0]')
        parser.add_option('-q', '--quanlty', action="store", dest="quanlty",metavar='num',default=0, help='set reads mapping quanlty for filter; [default: 0]')
        #parser.add_option('-I', '--identity', action="store", dest="identity",metavar=float,default=0.6, help='set min identity of mapping reads for filter; [default: 0.6]')
        parser.add_option('-F', '--outfmt', action="store", dest="outfmt",metavar='png/pdf',default='png', help='set out picture format; [default:png]')
        (options, args) = parser.parse_args()
        if not options.bams or not options.region:
            parser.print_help();
            sys.exit(0)
        self.bams = options.bams.split(',')
        for i in self.bams:isfile(i)
        self.region = options.region
        if ':' not in self.region and '-' not in self.region:
            print("region: chr:from-to")
            parser.print_help();
            sys.exit(0)
        self.outdir = os.path.realpath(options.outdir)
        self.thread = max(int(options.thread)+1,1)
        self.minMapq = int(options.quanlty)
        self.outfmt = options.outfmt
        if self.outfmt!="png" and self.outfmt !="pdf":
            print("Temporary outfmt only support png and pdf format!")
            parser.print_help();
            sys.exit(0)
        MakeDir(self.outdir)
        self.inbed = self.outdir+'/step.bed'

        self.beddepthdir = self.outdir+'/beddepthout'
        self.figuredir = self.outdir+'/figure'
        MakeDir(self.beddepthdir)
        MakeDir(self.figuredir)
        self.chrom = self.region.split(':')[0]
        self.start = int(self.region.split(':')[1].split('-')[0])
        self.end = int(self.region.split(':')[1].split('-')[1])

    def region2bed(self):
        regionlen = self.end-self.start
        windows = 1 if len(str(regionlen))<=3 else int('1'+'0'*(len(str(regionlen))-3))
        stepregion = [(self.start+i*windows,min(self.start+(i+1)*windows, self.end)) for i in range(int(math.ceil(regionlen/float(windows))))]
        fo = open(self.inbed,"w")
        for i in stepregion:
            fo.write("{}\t{}\t{}\n".format(self.chrom,str(i[0]),str(i[1])))
        fo.close()

    def getdepthfrombam(self):
        cmd = []
        regionname = self.region.replace(':','_').replace('-','_')
        alldepthfile = []
        pool = Pool(processes=self.thread)
        results = []
        for i in self.bams:
            baminame = os.path.basename(i).replace('.bam','')
            depthfilei = self.beddepthdir+'/'+baminame+'.'+regionname+'.depth'
            alldepthfile.append(depthfilei)
            cmdi = "samtools bedcov -Q {} {} {} -j >{}".format(str(self.minMapq), self.inbed,i,depthfilei)
            pool.apply_async(Runcmd,(cmdi,))

        pool.close()
        pool.join()
        if self.outfmt=="pdf":
            Rcmd = "Rscript {} {} {}".format(self.drawdepth,','.join(alldepthfile),self.figuredir+'/'+regionname+'.pdf')
        if self.outfmt=="png":
            Rcmd = "Rscript {} {} {}".format(self.drawdepth,','.join(alldepthfile),self.figuredir+'/'+regionname+'.png')
        os.system(Rcmd)

    def run(self):
        self.region2bed()
        self.getdepthfrombam()



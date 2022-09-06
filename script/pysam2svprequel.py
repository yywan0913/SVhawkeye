import sys,os
sys.path.append(os.path.dirname(os.path.realpath(sys.argv[0]))+'/script')
from  regioncluster import *
## all region reads base  when loci=True
def getlocibase(bf,chrom,start,end,pyref=None,min_mapq=20, min_baseq=0,no_del=False, no_dup=False):
    dbasepos = {}
    dbasepos.setdefault('ref',{})
    if pyref != None:
        reffa = pyref.fetch(chrom, start-1, end).upper()
        for i in range(end-start+1):
            dbasepos['ref'][start+i] = reffa[i]
    #for col in bf.pileup(chrom,start,end,truncate=True,stepper="nofilter"): #truncate=True  ignore_overlaps=True flag_filter=0
    for col in bf.pileup(reference=chrom, start=start-1, end=end, min_base_quality=0, min_mapping_quality=min_mapq,truncate=True):
        reads = col.pileups # allreads
        for read in reads:
            if read.alignment.is_reverse:
                strand = '-'
                query_start = read.alignment.query_length - read.alignment.query_alignment_end + 1
                query_end = read.alignment.query_length - read.alignment.query_alignment_start + 1
            else:
                strand = "+"
                query_start = read.alignment.query_alignment_start
                query_end = read.alignment.query_alignment_end
            readsname = read.alignment.query_name
            readsnamelabel = read.alignment.query_name+'__'+str(query_start)
            dbasepos.setdefault(readsname,{})
            dbasepos[readsname].setdefault(readsnamelabel,{})
            if read.is_del:
                base = '-'
            elif read.indel >0 :
                base = read.alignment.seq[read.query_position:(read.query_position+read.indel+1)]
            else :
                base = read.alignment.seq[read.query_position]
            dbasepos[readsname][readsnamelabel][col.pos+1] = base
	
    return dbasepos

###  cigar 2 SVinfo
## [(4,15),(0,11),(1,2),(2,1)] 
## M~0=BAM_CMATCH; I~1=BAM_CINS; D~2=BAM_CDEL; 
## S~4=BAM_CSOFT_CLIP; H~5=BAM_CHARD_CLIP
## N~3 P~6 =~7 X~8 ignord
class resolve_read(object):  # 0-base
    def __init__(self,pysam_read,start,end):
        self.cigar = pysam_read.cigar
        self.reference_start = pysam_read.reference_start
        self.read = pysam_read
        #self.NM = pysam_read.get_tag('NM')
        self.start = start
        self.end = end

    def _get_cigar(self):
        dict_cigar = {}
        absolute_pos = self.reference_start
        l_m = 0; l_d = 0; l_i = 0; l_o = 0
        for i in range(len(self.cigar)):
            if self.cigar[i][0] not in [1,4,5]:
                absolute_pos+=self.cigar[i][1]

            if self.cigar[i][0]==0: # match
                l_m += self.cigar[i][1]
                l_o += 1

            elif self.cigar[i][0]==4 or self.cigar[i][0]==5: #softclipping
                soft_cliping_len = self.cigar[i][1]
                if absolute_pos>self.end or absolute_pos<self.start : continue
                if soft_cliping_len >=50:
                    if i == 0:
                        if not self.read.is_reverse:
                            dict_cigar.setdefault('ssc',{}) # start softclipping
                            dict_cigar['ssc'][absolute_pos]=soft_cliping_len
                        else:
                            dict_cigar.setdefault('esc',{})
                            dict_cigar['esc'][absolute_pos]=soft_cliping_len
                    else:
                        if not self.read.is_reverse:
                            dict_cigar.setdefault('esc',{}) # end softclipping
                            dict_cigar['esc'][absolute_pos]=soft_cliping_len
                        else:
                            dict_cigar.setdefault('ssc',{})
                            dict_cigar['ssc'][absolute_pos]=soft_cliping_len

            elif self.cigar[i][0]==1: # ins
                ins_len = self.cigar[i][1]
                l_i += ins_len
                l_o += 1
                if absolute_pos>self.end or absolute_pos<self.start : continue
                if ins_len>=50:   ## AAGC(TTTTTT)GCGA  # absolute_pos=4,ins_len=6
                    dict_cigar.setdefault('ins',{})
                    dict_cigar['ins'][absolute_pos]=ins_len

            elif self.cigar[i][0]==2: # del
                del_len = self.cigar[i][1]
                l_d += del_len
                l_o += 1
                #if absolute_pos>self.end or absolute_pos<self.start : continue
                if not Overlap(self.start,self.end,absolute_pos,absolute_pos+del_len):continue
                if del_len >=50:   ## AATT-----GGACGA : absolute_pos=10,del_len=5
                    dict_cigar.setdefault('del',{})
                    dict_cigar['del'][absolute_pos]=del_len
                #else:
                #	if self.end-self.start<=5000:
                #		dict_cigar.setdefault('inDel',{})
                #		dict_cigar['inDel'][absolute_pos]=del_len						
            elif self.cigar[i][0]==3: # rna N splicing
                splic_len = self.cigar[i][1]
                #if absolute_pos>self.end or absolute_pos<self.start : continue
                #if not Overlap(self.start,self.end,absolute_pos,absolute_pos+splic_len):continue
                dict_cigar.setdefault('splic',{})
                dict_cigar['splic'][absolute_pos]=splic_len
            else:
                pass

	# gap_compressed_identity
        #try:
        #    self.NM = pysam_read.get_tag('NM')
        #except:
        #    self.NM = l_m + l_o
        #identity = 1 - float(self.NM - l_d - l_i + l_o)/(l_m + l_o);
        #dict_cigar['identity'] = identity
        dict_cigar['identity'] = 'NA'

        return dict_cigar		
		
    def cigar2sv(self,dsv):
        #dsv = {}
        dict_cigar = self._get_cigar()
        chromosome = self.read.reference_name
        readsname = self.read.query_name
        if self.read.is_reverse:
            strand = '-'
            query_start = self.read.query_length - self.read.query_alignment_end + 1
            query_end = self.read.query_length - self.read.query_alignment_start + 1
        else:
            strand = '+'
            query_start = self.read.query_alignment_start
            query_end = self.read.query_alignment_end

        readsnamelabel = self.read.query_name+'__'+str(query_start)
        dsv.setdefault(readsname,{})
        dsv[readsname].setdefault(readsnamelabel,{})

        Dsvname = dsv[readsname][readsnamelabel]  ## readsname:readsnamelabel:(strand,query_start,...,chromosome) 
        Dsvname['strand'] = strand
        Dsvname['query_start'] = query_start
        Dsvname['query_end'] = query_end
        Dsvname['ref_start'] = self.read.reference_start
        Dsvname['ref_end'] = self.read.reference_end
        Dsvname['reads_length']=self.read.query_length
        Dsvname['identity'] = dict_cigar['identity']
        Dsvname['mapq'] = self.read.mapping_quality
        Dsvname['chromosome'] = self.read.reference_name
        Dsvname.setdefault('labelinfo',[])

        if len(dict_cigar.keys())>1:
            if 'ssc' in dict_cigar:
                #Dsvname.setdefault('sc',{})
                Dsvname['ssc'] = list(dict_cigar['ssc'].items())[0]
            if 'esc' in dict_cigar:
                #Dsvname.setdefault('sc',{})
                Dsvname['esc'] = list(dict_cigar['esc'].items())[0]
            if 'ins' in dict_cigar:
                for ipos in dict_cigar['ins']:
                    svlabels = "%s--%s_%s-%s:%s" %('ins', chromosome, str(ipos), str(ipos+1), str(dict_cigar['ins'][ipos]))
                    Dsvname.setdefault('labelinfo',[]).append(svlabels)
            if 'del' in dict_cigar:
                for ipos in dict_cigar['del']:
                    svlen = dict_cigar['del'][ipos]
                    svlabels = "%s--%s_%s-%s:%s" %('del', chromosome, str(ipos-svlen), str(ipos), str(svlen))
                    Dsvname.setdefault('labelinfo',[]).append(svlabels)
            if 'inDel' in dict_cigar:
                for ipos in dict_cigar['inDel']:
                    svlen = dict_cigar['inDel'][ipos]
                    svlabels = "%s--%s_%s-%s:%s" %('inDel', chromosome, str(ipos-svlen), str(ipos), str(svlen))
                    Dsvname.setdefault('labelinfo',[]).append(svlabels)
            if 'splic' in dict_cigar:
                for ipos in dict_cigar['splic']:
                    svlen = dict_cigar['splic'][ipos]
                    svlabels = "%s--%s_%s-%s:%s" %('splic', chromosome, str(ipos-svlen), str(ipos), str(svlen))
                    Dsvname.setdefault('labelinfo',[]).append(svlabels)
        return dsv

# split_reads 2 SVinfo
class resolve_split_reads(object):
	def __init__(self,dsv):
		self.dsv = dsv

	def _singlereadsc(self,dictsingle):
		if 'ssc' in dictsingle:
			pos,sclen = dictsingle['ssc']
			svlabels = "%s--%s_%s-%s:%s" %('sc',dictsingle['chromosome'],str(pos),str(pos+1),str(sclen))
			dictsingle.setdefault('labelinfo',[]).append(svlabels)
		if 'esc' in dictsingle:
			pos,sclen = dictsingle['esc']
			svlabels = "%s--%s_%s-%s:%s" %('sc',dictsingle['chromosome'],str(pos),str(pos+1),str(sclen))
			dictsingle.setdefault('labelinfo',[]).append(svlabels)
		return dictsingle

	def _split_reads2sv(self,dictdup,sv_min_length): # dictdup = {readsname_x1:{ref_start:111,mapq:20,...},reads_x2:{ref_end:221,mapq:20,...},reads_x3:{...}}
		dupreadsi = sorted(dictdup.keys(), key = lambda x:int(x.split('__')[1]) )
		for j in range(len(dupreadsi)-1):    ## 0,1
			k = j + 1
			readj = dupreadsi[j]  ## has __ ;  name__32
			readk = dupreadsi[k]  ## hsa __ aa__568
			startj = dictdup[readj]['ref_start']
			endj = dictdup[readj]['ref_end']
			strandj = dictdup[readj]['strand']
			chrj = dictdup[readj]['chromosome']
			query_start_j = dictdup[readj]['query_start']
			query_end_j = dictdup[readj]['query_end']

			startk = dictdup[readk]['ref_start']
			endk = dictdup[readk]['ref_end']
			strandk = dictdup[readk]['strand']
			chrk = dictdup[readk]['chromosome']
			query_start_k = dictdup[readk]['query_start']
			query_end_k = dictdup[readk]['query_end']
			

			## ------------------ SC --------------------- #
			if j == 0:
				if 'sc' in dictdup[readj]:
					if 'ssc' in dictdup[readj]['sc']:
						jpos,sclen = dictdup[readj]['sc']['ssc']
						svlabels = "%s--%s_%s-%s:%s" %('SC',chrj,str(jpos),str(jpos+1),str(sclen))
						dictdup[readj].setdefault('labelinfo',[]).append(svlabels)
			if k == len(dupreadsi)-1:
				if 'sc' in dictdup[readk]:
					if 'esc' in dictdup[readk]['sc']:
						kpos,sclen = dictdup[readk]['sc']['esc']
						svlabels = "%s--%s_%s-%s:%s" %('SC',chrk,str(kpos),str(kpos+1),str(sclen))
						dictdup[readk].setdefault('labelinfo',[]).append(svlabels)
			## ------------------ TRA -------------------- #

			## ------------------ INV -------------------- #
			if strandj != strandk:
				if strandj=="+":
					svlen = abs(endj - endk) + 1
					svpos = (endj, endk) if endj < endk else (endk, endj)
				else:
					svlen = abs(startj - startk) + 1
					svpos = (startj, startk) if startj < startk else (startk, startj)

				if svlen >= sv_min_length:
					svlabels = "%s--%s_%s-%s:%s" %('INV',chrj,str(svpos[0]),str(svpos[1]),svlen)
					dictdup[readj].setdefault('labelinfo',[]).append(svlabels)
					dictdup[readk].setdefault('labelinfo',[]).append(svlabels)
			else:
				## ------------------ INS -------------------- #
				if query_start_k - query_end_j >= 100:
					#if abs(startk - endj) < query_start_k - query_end_j:
					svlen = query_start_k - query_end_j - 1
					if strandj == "+":
						svlabels = "%s--%s_%s-%s:%s" %('INS',chrk,str(startk),str(startk+1),svlen)
					else:
						svlabels = "%s--%s_%s-%s:%s" %('INS',chrk,str(endk-1),str(endk),svlen)
					#dictdup[readj].setdefault('sv',{})
					#dictdup[readk].setdefault('sv',{})
					dictdup[readj].setdefault('labelinfo',[]).append(svlabels)
					dictdup[readk].setdefault('labelinfo',[]).append(svlabels)

				## ------------------ DUP -------------------- #
				if not (startj > endk or endj < startk):
					if strandj == "+":
						svlen = abs(endj - startk) + 1
						svpos = (startk, endj) if startk < endj else (endj, startk)
					else:
						svlen = abs(endk - startj) + 1
						svpos = (startj, endk) if startj < endk else (endk, startj)
					if svlen >= sv_min_length:
						svlabels = "%s--%s_%s-%s:%s" %('DUP',chrj,str(svpos[0]),str(svpos[1]),svlen)
						dictdup[readj].setdefault('labelinfo',[]).append(svlabels)
						dictdup[readk].setdefault('labelinfo',[]).append(svlabels)

				## ------------------ DEL or DUP -------------------- #
				else:
					if startj>endk:
						if strandj == "+":
							svlen = abs(endj-startk)+1
							svpos = (startk, endj) if startk < endj else (endj, startk)

							if svlen>=sv_min_length:
								svlabels = "%s--%s_%s-%s:%s" %('DUP',chrj,str(svpos[0]),str(svpos[1]),svlen)
								dictdup[readj].setdefault('labelinfo',[]).append(svlabels)
								dictdup[readk].setdefault('labelinfo',[]).append(svlabels)
						else:
							svlen = abs(startj-endk)+1
							svpos = (startj, endk) if startj < endk else (endk, startj)

							if svlen >= sv_min_length:
								svlabels = "%s--%s_%s-%s:%s" %('DEL',chrj,str(svpos[0]),str(svpos[1]),svlen)
								dictdup[readj].setdefault('labelinfo',[]).append(svlabels)
								dictdup[readk].setdefault('labelinfo',[]).append(svlabels)
					else:
						if strandj == "+":
							svlen = abs(startk-endj)+1
							svpos = (startk, endj) if startk < endj else (endj, startk)

							if svlen>=sv_min_length:
								svlabels = "%s--%s_%s-%s:%s" %('DEL',chrj,str(svpos[0]),str(svpos[1]),svlen)
								dictdup[readj].setdefault('labelinfo',[]).append(svlabels)
								dictdup[readk].setdefault('labelinfo',[]).append(svlabels)
						else:
							svlen = abs(startj-endk)+1
							svpos = (startj, endk) if startj < endk else (endk, startj)

							if svlen >= sv_min_length:
								svlabels = "%s--%s_%s-%s:%s" %('DUP',chrj,str(svpos[0]),str(svpos[1]),svlen)
								dictdup[readj].setdefault('labelinfo',[]).append(svlabels)
								dictdup[readk].setdefault('labelinfo',[]).append(svlabels)

		return dictdup	

	def get_sp_sv(self,sv_min_length):
		for readname in self.dsv:
			readname_sp = list(self.dsv[readname].keys())
			if len(readname_sp)==1:
				dsinglereads = self.dsv[readname][readname_sp[0]]
				self.dsv[readname][readname_sp[0]] = self._singlereadsc(dsinglereads)
				continue
			dictdup = self.dsv[readname]
			self.dsv[readname] = self._split_reads2sv(dictdup,sv_min_length)
		return self.dsv

def Overlap(start1,end1,start2,end2):
	if not (start2>end1 or start1>end2):
		return True
	else:
		return False

## merge cigar and split_reads and base if loci=True
def dsv2info(dsv,bf,chrom,start,end,outfile,pyref=None,min_mapq=20,loci=True):
    dreadsnamelabelorder = RegionClusterOrder().SVeraseOverlapIntervals(dsv) # {'name3__12':'4','name2__4':'3','name4__66':1,'name8_33':2}...
    # head
    outhead1 = ['Chr','RefStart','RefEnd','QueryStart','QueryEnd','ReadsLen','Mapq','Identity','Strand','Color','Type','Readsorder','ReadsID']
    if loci:
        dbasepos = getlocibase(bf,chrom,start,end,pyref=pyref,min_mapq=min_mapq, min_baseq=0)
        if pyref!=None:
            outhead2 = [str(pos)+'_'+dbasepos['ref'][pos] for pos in range(start,end+1)]
        else:
            outhead2 = [str(pos) for pos in range(start,end+1)]
        out = ['\t'.join(outhead1+outhead2)]
    else:
        out = ['\t'.join(outhead1)]
    # content
    reads_n = 0
    for i in dsv: # i = readsname
        if len(dsv[i].keys())==1: 
            colorid = '-1'
        else:  # split reads
            reads_n += 1
            colorid = str(reads_n)
        for j in dsv[i]: # j = readsname_start
            #if j == "sv":continue
            dsvtmp = dsv[i][j]
            if loci:dbasetmp = dbasepos[i][j]
            if dsvtmp['labelinfo'] == []:
                readlabels = 'normal'
                colorid = '0'
            else:
                if any(['splic--' in km for km in dsvtmp['labelinfo']]):
                    colorid = '-2'
                readlabels = '@@'.join(dsvtmp['labelinfo'])
            readsorderj = dreadsnamelabelorder[j]
    
            #basic info and sv
            out1 = [dsvtmp['chromosome'],str(dsvtmp['ref_start']),str(dsvtmp['ref_end']),str(dsvtmp['query_start']),
                         str(dsvtmp['query_end']), str(dsvtmp['reads_length']), str(dsvtmp['mapq']), str(dsvtmp['identity']), 
                         dsvtmp['strand'], colorid, readlabels, readsorderj, i]
            if loci:
                # base info
                out2 = [dbasetmp[pos] if pos in dbasetmp else '+' for pos in range(start,end+1)]
                out.append('\t'.join(out1+out2))
            else:
                out.append('\t'.join(out1))

    fout = open(outfile, 'w')
    fout.write('\n'.join(out)+'\n')
    fout.close()
    return list(dsv.keys())

def Tosnpinfo(dsv,bf,chrom,start,end,outfile,pyref=None,min_mapq=20):
    dreadsnamelabelorder = RegionClusterOrder().SVeraseOverlapIntervals(dsv)
    outhead1 = ['Chr','RefStart','RefEnd','QueryStart','QueryEnd','ReadsLen','Mapq','Identity','Strand','Color','Type','Readsorder','ReadsID']
    dbasepos = getlocibase(bf,chrom,start,end,pyref=pyref,min_mapq=min_mapq, min_baseq=0)
    if pyref!=None:
        outhead2 = [str(pos)+'_'+dbasepos['ref'][pos] for pos in range(start,end+1)]
    else:
        outhead2 = [str(pos) for pos in range(start,end+1)]
    out = ['\t'.join(outhead1+outhead2)]
    reads_n = 0
    for i in dsv:
        if len(dsv[i].keys())==1:
            colorid = '0'
        else:
            reads_n += 1
            colorid = str(reads_n)
        for j in dsv[i]:
            dsvtmp = dsv[i][j]
            dbasetmp = dbasepos[i][j]
            readsorderj = dreadsnamelabelorder[j]
            out1 = [dsvtmp['chromosome'],str(dsvtmp['ref_start']),str(dsvtmp['ref_end']),str(dsvtmp['query_start']),
                    str(dsvtmp['query_end']), str(dsvtmp['reads_length']), str(dsvtmp['mapq']), str(dsvtmp['identity']),
                    dsvtmp['strand'], colorid, 'normal',readsorderj, i]
            out2 = [dbasetmp[pos] if pos in dbasetmp else '+' for pos in range(start,end+1)]
            out.append('\t'.join(out1+out2))
    fout = open(outfile, 'w')
    fout.write('\n'.join(out)+'\n')
    fout.close()



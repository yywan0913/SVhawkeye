def getlocibase(bf,chrom,start,end,pyref=None,min_mapq=20, min_baseq=0,no_del=False, no_dup=False):
	dbasepos = {}
	dbasepos.setdefault('ref',{})
	if pyref != None:
		reffa = pyref.fetch(chrom, start-1, end).upper()
		for i in range(end-start+1):
			dbasepos['ref'][start+i] = reffa[i]
		
	#for col in bf.pileup(chrom,start,end,truncate=True,stepper="nofilter"): #truncate=True  ignore_overlaps=True flag_filter=0
	for col in bf.pileup(reference=chrom, start=start, end=end+1, min_base_quality=0, min_mapping_quality=min_mapq,truncate=True):
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
			dbasepos[readsname][readsnamelabel][col.pos] = base
	
	return dbasepos

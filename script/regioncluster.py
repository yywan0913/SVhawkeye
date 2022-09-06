def intersectregion(start1,end1,start2,end2):
    if start2 > end1 or end2 < start1 :
        return False
    else:
        return True

class RegionClusterOrder:
    def eraseOverlapIntervals(self, intervals):
        if not intervals:
            return 0
        intervals.sort(key=lambda x:(x[1],x[0]))
		#[[2,4],[3,4],[1,3],[0,3],[2,3]]
		#to [[1, 3], [0, 3], [2, 3], [2, 4], [3, 4]]
        #right_end =intervals[0][1]
        dorder = {}
        N = 0    
        while intervals:
            N+=1
            dorder.setdefault(str(N),[]).append(intervals[0])
            diu = []
            n = len(intervals)
            right_end =intervals[0][1]
            for i in range(1, n):
                if intervals[i][0] >= right_end:
                    right_end =intervals[i][1]
                    dorder.setdefault(str(N),[]).append(intervals[i])
                else:
                    diu.append(intervals[i])
            intervals = diu
        return dorder

    def SVeraseOverlapIntervals(self, dictsv):
        if not dictsv:
            return 0
        svdictlist = []
        for i in dictsv.values():
            svdictlist.append({ list(i.keys())[0] : [ list(i.values())[0]['ref_start'], list(i.values())[0]['ref_end']]})
        svdictlist = sorted(svdictlist,key=lambda x:(list(x.values())[0][1],list(x.values())[0][0]))
        dorder = {}
        #realorder = {}
        N = 0
        while svdictlist:
            N+=1
            firstname = list(svdictlist[0].keys())[0].split('__')[0]
            if len(dictsv[firstname])==1: # not split mapping reads
                #dorder.setdefault(str(N),[]).append(svdictlist[0]) # first fill
                secondname = list(svdictlist[0].keys())[0]
                dorder[secondname] = str(N)
                right_end = list(svdictlist[0].values())[0][1]
                diu = []
            else: # split mapping reads
                #if dictsv[firstname][list(dictsv[firstname])[0]]['strand'] == "+":
                dsvchildname = sorted(dictsv[firstname],key=lambda x:int(x.split('__')[-1])) # namelable order: __7 __10 __ 22
                #else:
                #    dsvchildname = sorted(dictsv[firstname],key=lambda x:int(x.split('__')[-1]), reverse=True )
                m = 0
                for k in range(len(dsvchildname)):
                    secondname = dsvchildname[k]
                    if k==0:
                        dorder[secondname] = str(N)
                        secondname1 = secondname
                    else:
                        secondname2 = secondname
                        secondname1label = set(dictsv[firstname][secondname1]['labelinfo'])
                        secondname2label = set(dictsv[firstname][secondname2]['labelinfo'])
                        intersect = set(secondname1label).intersection(set(secondname2label))
                        read1start = dictsv[firstname][secondname1]['ref_start']
                        read1end = dictsv[firstname][secondname1]['ref_end']
                        read2start = dictsv[firstname][secondname2]['ref_start']
                        read2end = dictsv[firstname][secondname2]['ref_end']
                        if len(intersect)==0:
                            N+=1
                            dorder[secondname2] = str(N)
                        else:
                            if any(['DUP' in j for j in intersect]):
                                N+=1
                                dorder[secondname2] = str(N)
                            elif intersectregion(read1start,read1end,read2start,read2end):
                                N++1
                                dorder[secondname2] = str(N)
                            else:
                                m += 0.01
                                dorder[secondname2] = str(N+m)
                        secondname1 = secondname2
                right_end = 1000000000
                #del svdictlist[0]
                diu = []
                
            n = len(svdictlist)
            for i in range(1, n):
                firstname = list(svdictlist[i].keys())[0].split('__')[0]
                if len(dictsv[firstname])>1: #split mapping
                    diu.append(svdictlist[i])
                    continue
                lefti = list(svdictlist[i].values())[0][0]  # start
                if lefti >= right_end:
                    right_end = list(svdictlist[i].values())[0][1]
                    #dorder.setdefault(str(N),[]).append(svdictlist[i])
                    secondname = list(svdictlist[i].keys())[0]
                    dorder[secondname] = str(N)
                else:
                    diu.append(svdictlist[i])
            svdictlist = diu
        return dorder

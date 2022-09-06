import os
import sys



def intersectregion(start1,end1,start2,end2):
    if start2 > end1 or end2 < start1 :
        return False
    else:
        return True


def MakeDir(Dir):
    if not os.path.exists(Dir):
        os.makedirs(Dir)

def isfile(File):
    if not os.path.isfile(File):
        print('Error:%s file is not existed!' %(File))
        sys.exit(1)
    else :
        return 1


def writefile(strings,outfile):
    f=open(outfile,"w")
    f.write(strings+'\n')
    f.close()

def _TF(x):
    y = True if x==True or x=="True" or x=="T" or x=="1" else False
    return y

def Runcmd(cmd):
    os.system(cmd)

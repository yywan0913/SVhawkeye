#!/usr/bin/env python
import os,sys
from sys import version_info

def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

if version_info.major == 2:
    import copy_reg
    import types
    copy_reg.pickle(types.MethodType, _pickle_method)

def Help(argv):
    print("python %s <command> [option]"%os.path.basename(argv[0]));
    command = '''
commands:
    sv_browse             fast draw Structural variation or snp-inDel as IGV.
    snpindel_browse       fast draw snp or indel variation as IGV
    sv_genotyping         recall sv of existing input vcf file.
    rna_browse            display isoform structure from iso-seq
    regiondepth_browse    display depth distribution of your region

    '''
    print(command);
    sys.exit(0)


if __name__ == '__main__':
    if len(sys.argv) < 2 or (sys.argv[1] != 'sv_browse' and sys.argv[1] != 'sv_genotyping' 
            and sys.argv[1]!='regiondepth_browse' and sys.argv[1]!="rna_browse" and sys.argv[1]!="snpindel_browse"):
        Help(sys.argv)
        sys.exit(0);
    if sys.argv[1] == 'sv_browse': 
        from script.svanalysispipe import *
        pipe1 = SVToView();
        pipe1.GetOpt();
        pipe1.runsvbrowser();
    if sys.argv[1] == "snpindel_browse":
        from script.snpanalysispipe import *
        pipe2 = SNPToView();
        pipe2.GetOpt();
        pipe2.runsnpbrowser();
    if sys.argv[1] == 'sv_genotyping':
        from script.svanalysispipe import *
        pipe3 = SVToView();
        pipe3.GetOpt(genotyping=1);
        pipe3.runsvgenotyping();
    if sys.argv[1] == 'regiondepth_browse':
        pipe4 = '...'
    if sys.argv[1] == 'rna_browse':
        from script.isoformanalysispipe import *
        pipe5 = RNAToView();
        pipe5.GetOpt();
        pipe5.runrnabrowser();

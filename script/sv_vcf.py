"""
A Universal Stucture Variant VCF parsing module
tested vcf: sniffles vcf, nanosv vcf, picky vcf
shared INFO ID are: SVTYPE, END, SVLEN
RE(reads evidence): sniffles, picky; nano SV: RT(2d,template,complement)
BND shared format: N[ref:pos2[

BND format:
    N]chr6:25647927]
STAT  REF  ALT   Meaning
s1    s    t[p[  piece extending to the right of p is joined after t
s2    s    t]p]  reverse comp piece extending left of p is joined after t
s3    s    ]p]t  piece extending to the left of p is joined before t
s4    s    [p[t  reverse comp piece extending right of p is joined before t
"""
import logging


class bnd(object):
    """
    just a bnd classifier
    """
    def __init__(self,bnd_string):
        self.bnd_string = bnd_string
        # N[chr1:123456[, N]chr1:123456], ]chr1:123456]N, [chr1:123456[N
        if '[' in self.bnd_string:
            # there should be two '['(right_braket) in a BND string
            first_right_braket = self.bnd_string.index('[')
            next_right_braket = self.bnd_string.index('[', first_right_braket+1)
            if first_right_braket != 0:
                self.stat = "s1"
            else:
                self.stat = "s4"
            self.pos = self.bnd_string[first_right_braket+1:next_right_braket]
        elif ']' in self.bnd_string:
            # there should be two ']'(left_braket) in a BND string
            first_left_braket = self.bnd_string.index(']')
            next_left_braket = self.bnd_string.index(']', first_left_braket+1)
            if first_left_braket != 0:
                self.stat = "s2"
            else:
                self.stat = "s3"
            self.pos = self.bnd_string[first_left_braket+1:next_left_braket]
        self.chrom = self.pos.split(":")[0]
        self.pos_num = self.pos.split(":")[1]


class sv_vcf_record(object):
    def __init__(self,record):
        self.record = record

        fields = record.strip().split("\t")
        # checked that ID is unique for sniffles, nanosv, picky vcf
        (self.chrom1,self.pos1,self.id,self.ref,self.alt,self.qual,self.filter,
                self.info,self.format) = fields[:9]
        self.sample_formats = fields[9:]

        # info dict
        self.info_dict = {}
        info_list = self.info.split(";")
        for i in info_list:
            if "=" in i:
                info_id,info_value = i.split("=")
                self.info_dict[info_id] = info_value
            else:
                self.info_dict[i] = i

        self.svtype = self.info_dict["SVTYPE"]
        # BND, svtype, pos2
        if self.svtype == "BND":
            bnd_pos = bnd(self.alt) # BND position
            self.chrom2 = bnd_pos.chrom
            self.pos2 = bnd_pos.pos_num
            if (self.chrom1 == self.chrom2 and (bnd_pos.stat == "s2"
                    or bnd_pos.stat == "s4")): # INV
                self.svtype = "INV"
            elif self.chrom1 != self.chrom2:
                self.svtype = "TRA"
            else:
                # raise RuntimeError("bad line {}".format(record)) not raise this error now
                pass
        elif self.svtype == "TRA": # sniffles TRA (BND not specified)
            self.chrom2 = self.info_dict["CHR2"]
            self.pos2 = self.info_dict["END"]
        else:
            self.chrom2 = self.chrom1
            self.pos2 = self.info_dict["END"]
            # exchange pos1 and pos2, if pos1 > pos2
            if int(self.pos1) > int(self.pos2):
                tmp = self.pos1
                self.pos1 = self.pos2
                self.pos2 = tmp

        # svlen
        try:
            self.svlen = abs(float(self.info_dict["SVLEN"]))
        except KeyError:
            # INS, TRA do not have SVLEN attribute
            self.svlen = 0 if self.chrom1!=self.chrom2 else int(self.pos2)-int(self.pos1)

        # RE(number of read evidence)
        if "RE" in self.info_dict: # sniffles and picky
            self.re = self.info_dict["RE"]
        elif "RT" in self.info_dict: # nanosv
            # RT=2d,template,compementary; no matter what kind of reads they
            # are, add them up
            self.re = sum([int(i) for i in self.info_dict["RT"].split(",")])
            self.re = str(self.re)
        else:
            self.re = "NA"
            #logging.warning("Can not get RE(support reads num) "
            #"for {}".format(record))


    @property
    def sv_dict(self):

        _sv_dict = {}

        key_string = "{}_{}_{}-{}_{}".format(
                self.svtype,
                self.chrom1,
                self.pos1,
                self.chrom2,
                self.pos2)
        _sv_dict[key_string] = self

        return _sv_dict


def main():
    #test
    pass

if __name__ == "__main__":
    main()


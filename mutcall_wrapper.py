#!/usr/bin/env python

import os
import sys
import argparse
import tempfile
import subprocess
import string
import shutil
import traceback
from multiprocessing import Pool

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def cmd_caller(i):
    ar, cwd = i
    cmd, out = ar
    print "Running CMD", cmd
    #subprocess.check_call(cmd, shell=True, cwd=cwd) 
    return out

class MutCallerWrapper:
    def __init__(self, args):
        self.args = args

    def getTemplateDict(self):
        return dict(
            TOOL_DIR=os.path.join(BASE_DIR, "tools"),
            NORMAL_BAM=self.args.normal,
            TUMOR_BAM=self.args.tumor,
            COSMIC_VCF=self.args.cosmic_vcf,
            DBSNP_VCF=self.args.dbsnp_vcf,
            REF_SEQ=self.args.refseq,
            BLOCK_SIZE=self.args.blocksize,
            OUT_DIR=self.args.workdir)

    def check(self, params):
        raise Exception("Implement check to see if your tool can be run here")

    def run(self):
        params = self.getTemplateDict()

        self.check(params)

        p = Pool(self.args.cpus)
        values = p.map(cmd_caller, list( (a, params['OUT_DIR']) for a in self.run_map(params) ))
    
        cmd = self.run_reduce(params, values)
        if cmd is not None:
            print "Running CMD", cmd
            #subprocess.check_call(cmd, shell=True, cwd=params['OUT_DIR']) 

    def run_map(self, params):
        raise Exception("Implement a command line generator for calculations here")

    def run_reduce(self, params, values):
        raise Exception("Implement a command line for finalizing calculations here")



class Mutect(MutCallerWrapper):
    def run_map(self,params):
        
        template = """java -jar ${TOOL_DIR}/muTect-1.1.4.jar \
--analysis_type MuTect \
--reference_sequence ${REF_SEQ} \
--cosmic ${COSMIC_VCF} \
--dbsnp ${DBSNP_VCF} \
--input_file:normal ${NORMAL_BAM} \
--input_file:tumor ${TUMOR_BAM} \
--out ${OUT_DIR}/stats.txt \
--coverage_file ${OUT_DIR}/coverage.txt
"""
        cmd = string.Template(template).substitute( params )
        yield cmd, None

    #def run_reduce(self,params):
    #    return None

    def check(self, params):
        if not os.path.exists(os.path.join(params['TOOL_DIR'], "MuSEv0.9.8.6")):
            raise Exception("Can't find MuSEv0.9.8.6")

        if params['COSMIC_VCF'] is None:
            raise Exception("Missing COSMIC VCF")

        if params['DBSNP_VCF'] is None:
            raise Exception("Missing dnSNP VCF")


class Muse(MutCallerWrapper):
    
    def check(self, params):
        if not os.path.exists(os.path.join(params['TOOL_DIR'], "muTect-1.1.4.jar")):
            raise Exception("Can't find muTect-1.1.4.jar")
    
    def run_map(self, params):
        seq_map = {}
        with open( params['REF_SEQ'] + ".fai" ) as handle:
            for line in handle:
                tmp = line.split("\t")
                seq_map[tmp[0]] = long(tmp[1])

        template = "${TOOL_DIR}/MuSEv0.9.8.6 call -P MuSEv0.9.8.6 -p 0.05 -b 0.0001 -B -f ${REF_SEQ} ${TUMOR_BAM} ${NORMAL_BAM} -l ${INTERVAL_FILE}"
        
        counter = 0
        for seq in seq_map:
            l = seq_map[seq]
            for i in xrange(1, l, params['BLOCK_SIZE']):
                interval="%s:%s-%s" % (seq, i, min(i+params['BLOCK_SIZE']-1, l))
                interval_file = os.path.join(params['OUT_DIR'], "interval.%s" % (counter))
                counter += 1
                with open(interval_file, "w") as handle:
                    handle.write("%s\n" % interval) 
                cmd = string.Template(template).substitute( dict(params, INTERVAL_FILE=interval_file ) )
                yield (cmd, interval_file)
    
    def run_reduce(self,params, values):
        print "Reduce:", values
        return None

method_callers = {
    'mutect' : Mutect,
    'muse' : Muse
}



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-name", help="Name of the Sample", default="SAMPLE")
    parser.add_argument("--normal", help="Normal Tissue BAM file", required=True)
    parser.add_argument("--tumor", help="Tumor Tissue BAM file", required=True)
    parser.add_argument("--cosmic-vcf", help="COSMIC mutations in VCF")
    parser.add_argument("--dbsnp-vcf", help="dbSNP in VCF format")
    parser.add_argument("--method", choices=['mutect', 'muse'], required=True)
    parser.add_argument("--refseq")
    parser.add_argument("--blocksize", type=long, default=5000000L)
    parser.add_argument("--outdir", default="out")
    parser.add_argument("--cpus", type=int, default=8)    
    parser.add_argument("--clean-on-fail", action="store_true", default=False)
    args = parser.parse_args()

    if os.path.exists(args.outdir):
        raise Exception("Output directory already exists")

    args.workdir = os.path.abspath(tempfile.mkdtemp(dir="./", prefix="variant_calling_%s_" % (args.sample_name)))

    try:
        method = method_callers[args.method](args)
        method.run()
        shutil.move(args.workdir, args.outdir)
    except:
        traceback.print_exc()
        print "Failure"
        if args.clean_on_fail:
            shutil.rmtree(args.workdir)
        sys.exit(1)
    sys.exit(0)

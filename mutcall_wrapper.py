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
import vcf

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def cmd_caller(i):
    block, ar, cwd = i
    cmd, out = ar
    if isinstance(cmd,basestring):
        #print "Running CMD", cmd
        p = subprocess.Popen(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            with open(os.path.join(cwd, "%s.error" % (block)), "w") as handle:
                handle.write(stderr)
    else:
        cmd(block)
    return out


def fai_chunk(path, blocksize):
    seq_map = {}
    with open( path ) as handle:
        for line in handle:
            tmp = line.split("\t")
            seq_map[tmp[0]] = long(tmp[1])
    
    for seq in seq_map:
        l = seq_map[seq]
        for i in xrange(1, l, blocksize):
            yield (seq, i, min(i+blocksize-1, l)) 

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
        values = p.map(cmd_caller, list( (i, a, params['OUT_DIR']) for i, a in enumerate(self.run_map(params)) ), 1)
    
        cmd = self.run_reduce(params, values)
        if cmd is not None:
            if isinstance(cmd,basestring):
                print "Running CMD", cmd
                subprocess.check_call(cmd, shell=True, cwd=params['OUT_DIR']) 
            else:
                cmd(params)
                
    def run_map(self, params):
        raise Exception("Implement a command line generator for calculations here")

    def run_reduce(self, params, values):
        raise Exception("Implement a command line for finalizing calculations here")



class Mutect(MutCallerWrapper):
    def run_map(self,params):
        
        template = """java -Xmx4g -XX:ParallelGCThreads=2 -jar ${TOOL_DIR}/muTect-1.1.5.jar \
--analysis_type MuTect \
--intervals ${INTERVAL} \
--reference_sequence ${REF_SEQ} \
--cosmic ${COSMIC_VCF} \
--dbsnp ${DBSNP_VCF} \
--input_file:normal ${NORMAL_BAM} \
--input_file:tumor ${TUMOR_BAM} \
--out ${OUT_DIR}/stats.${BLOCK_NUM}.txt \
--coverage_file ${OUT_DIR}/coverage.${BLOCK_NUM}.txt \
--vcf ${OUT_DIR}/out.${BLOCK_NUM}.vcf"""
        for i, block in enumerate(fai_chunk( params['REF_SEQ'] + ".fai", params['BLOCK_SIZE'] ) ):
            if block[0] != "hs37d5":
                cmd = string.Template(template).substitute( dict(params, BLOCK_NUM=i, INTERVAL="%s:%s-%s" % (block[0], block[1], block[2]) ) )
                yield cmd, "%s/out.%s.vcf" % (params['OUT_DIR'], i)

    def run_reduce(self,params, values):
        print "Reduce:", values
        
        def r(block):
            vcf_writer = None
            for a in values:
                vcf_reader = vcf.Reader(filename=a)
                if vcf_writer is None:
                    vcf_writer = vcf.Writer(open(os.path.join(params['OUT_DIR'], "out.MuTect.vcf"), "w"), vcf_reader)
                for record in vcf_reader:
                    vcf_writer.write_record(record)                   
        return r


    def check(self, params):
        if not os.path.exists(os.path.join(params['TOOL_DIR'], "muTect-1.1.5.jar")):
            raise Exception("Can't find muTect-1.1.5.jar")

        if params['COSMIC_VCF'] is None:
            raise Exception("Missing COSMIC VCF")

        if params['DBSNP_VCF'] is None:
            raise Exception("Missing dnSNP VCF")


class Muse(MutCallerWrapper):
    
    def check(self, params):
        if not os.path.exists(os.path.join(params['TOOL_DIR'], "MuSEv0.9.8.6")):
            raise Exception("Can't find MuSEv0.9.8.6")
    
    def run_map(self, params):
        seq_map = {}
        with open( params['REF_SEQ'] + ".fai" ) as handle:
            for line in handle:
                tmp = line.split("\t")
                seq_map[tmp[0]] = long(tmp[1])

        template = "${TOOL_DIR}/MuSEv0.9.8.6 call -P varcall -p 0.05 -b 0.0001 -B -f ${REF_SEQ} ${TUMOR_BAM} ${NORMAL_BAM} -l ${INTERVAL_FILE}"
        
        counter = 0
        for seq in seq_map:
            l = seq_map[seq]
            for i in xrange(1, l, params['BLOCK_SIZE']):
                interval="%s:%s-%s" % (seq, i, min(i+params['BLOCK_SIZE']-1, l))
                interval_file = os.path.join(params['OUT_DIR'], "interval.%s" % (counter))
                with open(interval_file, "w") as handle:
                    handle.write("%s\n" % interval) 
                cmd = string.Template(template).substitute( dict(params, INTERVAL_FILE=interval_file ) )
                yield (cmd, "%s/varcall_interval.%s.MuSE.txt" % (params['OUT_DIR'], counter))
                counter += 1
    
    def run_reduce(self,params, values):
        print "Reduce:", values
        
        def r(block):
            with open(os.path.join(params['OUT_DIR'], "out.MuSE.txt"), "w") as handle:
                for a in values:
                    with open(a) as inhandle:
                        for line in inhandle:
                            handle.write(line)       
        return r

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

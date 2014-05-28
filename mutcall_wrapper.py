#!/usr/bin/env python

import os
import sys
import argparse
import tempfile
import subprocess
import string
import shutil
import traceback

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

class MutationWrapper:
    def __init__(self, args):
        self.args = args

    def getTemplateDict(self):
        return dict(
            BASE_DIR=BASE_DIR,
            NORMAL_BAM=self.args.normal,
            TUMOR_BAM=self.args.tumor,
            COSMIC_VCF=self.args.cosmic_vcf,
            DBSNP_VCF=self.args.dbsnp_vcf,
            REF_SEQ=self.args.refseq,
            OUT_DIR=self.args.workdir)

    def check(self, params):
        raise Exception("Implement check to see if your tool can be run here")

    def run(self):
        params = self.getTemplateDict()

        self.check(params)

        for cmd in self.run_map(params):
            print "Running CMD", cmd

        cmd = self.run_reduce(params)
        if cmd is not None:
            print "Running CMD", cmd

    def run_map(self, params):
        raise Exception("Implement a command line generator for calculations here")

    def run_reduce(self, params):
        raise Exception("Implement a command line for finalizing calculations here")



class Mutect(MutationWrapper):
    def run_map(self,params):
        template = """java -jar ${BASE_DIR}/muTect-1.1.4.jar
--analysis_type MuTect
--reference_sequence ${REF_SEQ}
--cosmic ${COSMIC_VCF}
--dbsnp ${DBSNP_VCF}
--input_file:normal ${NORMAL_BAM}
--input_file:tumor ${TUMOR_BAM}
--out ${OUT_DIR}/stats.txt
--coverage_file ${OUT_DIR}/coverage.txt
"""
        cmd = string.Template(template).substitute( params )
        yield cmd

    #def run_reduce(self,params):
    #    return None

    def check(self, params):
        if not os.path.exists(os.path.join(params['BASE_DIR'], "muTect-1.1.4.jar")):
            raise Exception("Can't find muTect-1.1.4.jar")

method_callers = {
    'mutect' : Mutect
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
    parser.add_argument("--outdir", default="out")
    parser.add_argument("--clean-on-fail", action="store_true", default=False)

    args = parser.parse_args()

    if os.path.exists(args.outdir):
        raise Exception("Output directory already exists")

    args.workdir = tempfile.mkdtemp(dir="./", prefix="variant_calling_%s_" % (args.sample_name))

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

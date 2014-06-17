#!/bin/bash

wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/b37_cosmic_v54_120711.vcf
wget http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/dbsnp_132_b37.leftAligned.vcf.gz
gunzip dbsnp_132_b37.leftAligned.vcf.gz

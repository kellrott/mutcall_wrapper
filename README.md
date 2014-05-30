
Mutation Caller Wrapper

A script to wrap multiple mutation callers so that they use the same command line 
format. 

Change which mutation caller to use with the '--method' flag.

Example:

```
./mutcall_wrapper.py --refseq /pod/podstore/projects/PAWG/reference/genome.fa.gz \
--tumor /pod/pstore/projects/ICGC/e8f56d0f-eee4-4def-a43a-dec91f4382a1/f39e3105-2509-4dbd-80f4-71b5d417e27b/PAWG.568cd4c8-c91d-483d-ba55-62b76677871f.bam \
--normal /pod/pstore/projects/ICGC/e8f56d0f-eee4-4def-a43a-dec91f4382a1/6dfb4980-a8ef-4f57-90b9-490f170e95ad/PAWG.e4693c31-87c9-44ca-934a-29f6149a9a4f.bam \
--method muse
```


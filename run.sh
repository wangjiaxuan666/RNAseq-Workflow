#!/bin/bash

#---------Wang Jiaxuan 2022.03.08------------
# if run_RNAseq run failured but don't know why
# can use the code below this line
#--------------------------------------------


# java \
# -jar WDL/cromwell-58.jar \
# run WDL/rna_seq.wdl \
# -i WDL/input.RNAseq.json \
# -m WDL/output.RNAseq.json

# 但是这样的你的结果会在/data3/Group7/wangjiaxuan/workflow/bulk-rna-seq/cromwell-executions/RNAseq/e5a3dcfc-440d-4efc-a2b3-5b8cdffde16a/call-genome_map/shard-0/execution/ERR188044.sorted.bam
# 类似于这里面
./run_RNAseq -i WDL/input.RNAseq.json


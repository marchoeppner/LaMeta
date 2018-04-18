#!/bin/bash
for every sample in $allsets
do
echo \$sample >> tmp1
done

awk '{printf " -1 " \$3 " -2 " \$4 " -r " \$5}' tmp1 > $outcontigs

#echo "$MEGAHIT \$(cat tmp | tr -d '\n') --num-cpu-threads ${task.cpus} --presets meta-large -o megahit_out --mem-flag 2 --verbose
#mv megahit_out/final.contigs.fa $outcontigs
#mv megahit_out/log $megahitlog"

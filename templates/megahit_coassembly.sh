#!/bin/bash
echo $left_decon > out
echo $right_decon >> out
echo $unpaired_decon >> out

awk '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = \$i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' out > tmp1

awk '{printf " -1 " \$1 " -2 " \$2 " -r " \$3}' tmp1 > tmp

$MEGAHIT \$(cat tmp | tr -d '\n') --num-cpu-threads ${task.cpus} --presets meta-large -o megahit_out --mem-flag 2 --verbose
mv megahit_out/final.contigs.fa $outcontigs
mv megahit_out/log $megahitlog

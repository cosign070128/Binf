#!/bin/bash

for i in `ls *_extract`
do
  name=`ls $i | awk -F '_ex' '{print $1}'`
  echo $name
  awk '{print $3}' ./$i | sed '1,2d' | awk '$0=NR"F\t"$0' | sed "s/^/>"$name"_/g" | sed 's/\t/\n/g' > ./"$name".fasta
  awk '{print $5}' ./$i | sed '1,2d' | awk '$0=NR"R\t"$0' | sed "s/^/>"$name"_/g" | sed 's/\t/\n/g' >> ./"$name".fasta
done

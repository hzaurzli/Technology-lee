#!/bin/bash
for file in *.fna
do
echo "processing $file"
blastn -db ./89k/89k -query $file -out ${file%.*}.out -outfmt 6 -evalue 1e-5 -max_target_seqs 10000
done

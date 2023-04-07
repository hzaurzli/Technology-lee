#!/bin/bash

source ~/miniconda3/bin/activate blast;
for file in *.faa
do
echo "processing $file"
blastp -db /home/rzli/ssuis/card/index/protein_ARGs_length -query $file -out ${file%.*}.ARG -outfmt 6 -evalue 1e-5
done;
conda deactivate;

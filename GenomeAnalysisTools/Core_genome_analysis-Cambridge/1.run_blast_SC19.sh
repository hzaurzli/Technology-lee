for file in *.fna
do
echo "processing $file"
blastn -db ./sc19 -query $file -out ${file%.*}.sc19 -outfmt 6 -evalue 1e-5 -max_target_seqs 100000 # "-max_target_seqs 100000" is vary important, if not set this parameter, some contigs will not mapping all potential reference genes, particularly for the completed genome.
done

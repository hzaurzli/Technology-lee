for file in *.fna
do
echo "processing $file"
blastn -db ./Ssuis_Serotyping_len -query $file -out ${file%.*}.Sero -outfmt 6 -evalue 1e-5
done

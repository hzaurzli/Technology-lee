for file in *.out
do
echo "processing $file"
python 100100_blast_filter.py -i $file -o ${file%.*}.filter
done

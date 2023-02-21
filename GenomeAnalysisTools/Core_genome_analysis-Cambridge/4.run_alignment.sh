for file in *.fa
do
echo "processing $file"
/Users/zougeng/Documents/Software/clustalo -i $file -o ${file%.*}.aln -v #install the clustalo to your own device
done

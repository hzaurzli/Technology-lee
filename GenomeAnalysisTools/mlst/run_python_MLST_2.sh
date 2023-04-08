for file in *.filter
do
echo "processing $file"
python2 ST_analysis.py -i $file
done

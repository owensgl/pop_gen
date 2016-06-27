ls | grep processed.bam > list.txt
rm output.txt
while read name
do
echo -n "$name ">> output.txt
samtools flagstat $name  >> output.txt
done < list.txt

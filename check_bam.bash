ls | grep final.bam | grep -v bai > list.txt
rm output.txt
while read name
do
echo -n "$name ">> output.txt
samtools flagstat $name  > tmp.txt
cat tmp.txt | head -n 1 | cut -f 1 -d " " | tr '\n' '\t' >> output.txt
cat tmp.txt | head -n 5 | tail -n 1 | cut -f 1 -d " " >> output.txt
done < list.txt

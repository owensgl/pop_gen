ls | grep processed.bam > samplelist.txt
rm stats.txt
while read sample
do
samtools flagstat $sample > tmp.txt
echo -n -e "$sample\t" >> stats.txt
cut -f 1 -d ' ' tmp.txt | head -n 1 | xargs echo -n >> stats.txt
echo -n -e "\t" >> stats.txt
cut -f 1 -d ' ' tmp.txt | sed '1,2d' | head -n 1 | xargs echo -n >> stats.txt
echo -n -e "\n" >> stats.txt
done < samplelist.txt

rep=$1
while read n
do
#Make locations of junctions
Rscript /home/owens/bin/Hel_hybrids/make_junctions.R $n $rep

#Make Fake data
cat helianthus.est.full.cm.hmp | perl /home/owens/bin/pop_gen/make_artificial_hybrids.pl Helianthus.hybrid.specieslist.txt junction_locations_${n}cm.${rep}.txt > Helianthus.artificialhybrids.${n}cm.${rep}.hmp

#Make reverse fake data
head  -n 1 Helianthus.artificialhybrids.${n}cm.${rep}.hmp > tmp.head.$n.${rep}
sed '1d' Helianthus.artificialhybrids.${n}cm.${rep}.hmp | tac > tmp.body.$n.${rep}
cat tmp.head.$n.${rep} tmp.body.$n.${rep} > Helianthus.artificialhybrids.${n}cm.${rep}.rev.hmp
rm tmp.head.$n.${rep}
rm tmp.body.$n.${rep}

#Run HMM forward
cat Helianthus.artificialhybrids.${n}cm.${rep}.hmp | perl /home/owens/bin/pop_gen/hmp2hybviterbi_forward_v2.pl Helianthus.hybrid.specieslist.txt Helianthus.est.full.cm.parentalfreq.txt > Helianthus.artificialhybrids.${n}cm.${rep}.viterbi.foward.txt

#Run HMM reverse
cat Helianthus.artificialhybrids.${n}cm.${rep}.rev.hmp | perl /home/owens/bin/pop_gen/hmp2hybviterbi_reverse_v2.pl Helianthus.hybrid.specieslist.txt Helianthus.est.full.cm.parentalfreq.txt | tac> Helianthus.artificialhybrids.${n}cm.${rep}.viterbi.rev.txt

#Merge likelihoods
paste Helianthus.artificialhybrids.${n}cm.${rep}.viterbi.foward.txt Helianthus.artificialhybrids.${n}cm.${rep}.viterbi.rev.txt | perl /home/owens/bin/pop_gen/merge_forrev_viterbi.pl | perl /home/owens/bin/pop_gen/measure_accuracy.pl Helianthus.artificialhybrids.${n}cm.${rep}.hmp > tmp.${n}.${rep}.txt

echo -e  "\t$n" >> tmp.${n}.${rep}.txt
cat tmp.${n}.${rep}.txt >> accuracy.txt

rm tmp.${n}.${rep}.txt
rm junction_locations_${n}cm.${rep}.txt
rm Helianthus.artificialhybrids.${n}cm.${rep}.hmp
rm Helianthus.artificialhybrids.${n}cm.${rep}.rev.hmp
rm Helianthus.artificialhybrids.${n}cm.${rep}.viterbi.foward.txt
rm Helianthus.artificialhybrids.${n}cm.${rep}.viterbi.rev.txt

done < junction_sizes.txt

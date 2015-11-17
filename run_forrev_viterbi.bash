
cat helianthus.est.full.cm.hmp | perl /home/owens/bin/pop_gen/hmp2hybviterbi_forward.pl Helianthus.hybrid.specieslist.txt > Helianthus.hybrid.viterbi.0.1.cm.forward.txt

cat helianthus.est.full.cm.rev.hmp | perl /home/owens/bin/pop_gen/hmp2hybviterbi_reverse.pl Helianthus.hybrid.specieslist.txt | tac> Helianthus.hybrid.viterbi.cm.rev.txt

paste Helianthus.hybrid.viterbi.0.1.cm.forward.txt Helianthus.hybrid.viterbi.cm.rev.txt | perl /home/owens/bin/pop_gen/merge_forrev_viterbi.pl > Helianthus.hybrid.viterbi.cm.merged.txt

cat Helianthus.hybrid.viterbi.cm.merged.txt | perl /home/owens/bin/pop_gen/find_viterbi_windows.pl > Helianthus.hybrid.viterbi.cm.windows.txt

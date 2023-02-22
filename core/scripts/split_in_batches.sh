#! /bin/bash

plasmid_paths=plasmid_paths.txt
batch_size=100
plasmid_nr=`wc -l < $plasmid_paths`
counter=1
counter_lim=$(($plasmid_nr / $batch_size))

mkdir batches

echo
echo ">>> Spliting filepaths into batches ... <<<"
while [ $counter -le $counter_lim ]; do

	echo -en "\tBatch $counter"
	from=$[$counter*$batch_size - $batch_size + 1]
	to=$[$counter*$batch_size]
	echo " -- ( from: $from, to: $to )"
	sed -n -e ${from},${to}p $plasmid_paths > batches/batch_$counter

	counter=$[$counter + 1]

done

echo -en "\tLast batch $counter"
from=$[$counter*$batch_size - $batch_size + 1]
to=$(($plasmid_nr % $batch_size + $from - 1))
echo " -- ( from: $from to: $to )"
sed -n -e ${from},${to}p $plasmid_paths > batches/batch_$counter

# get all combinations of batches:
echo
touch batches_comb.txt
echo ">>> Making combinations of batches in batches_comb.txt ... <<<"
for i in $(seq 1 1 $counter); do
	for j in $(seq $i 1 $counter); do
		echo "batch_$i batch_$j" >> batches_comb.txt
	done
done

echo "DONE!"
echo
echo
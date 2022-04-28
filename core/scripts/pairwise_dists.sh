#! /bin/bash

batches_comb=batches_comb.txt
bindash=/path/to/bindash

counter=1
mkdir bindash

echo
echo "Sketching batches"
lines=`ls batches/* | wc -l`

for batch in batches/*; do

	$bindash sketch --listfname="$batch" --outfname=bindash/batch_${counter}_full_sketch --nthreads=20 --minhashtype=-1
	counter=$[$counter + 1]

	echo
	echo
	echo "~~~~~~~~~~~~"
	echo " $counter / $lines"
	echo "~~~~~~~~~~~~"

done
echo
echo "DONE!!!"
echo
echo
echo "Estimating distances ....."
calc(){ awk 'BEGIN { printf "%.3f", '$*' }'; }
counter=0
lines=`wc -l < $batches_comb`

while read line; do
	batch_1=`echo $line | cut -d' ' -f1`
	batch_2=`echo $line | cut -d' ' -f2`
	$bindash dist  bindash/${batch_1}_full_sketch bindash/${batch_2}_full_sketch --nthreads=20 --mthres=1e9 >> bindash_out.tsv

	progress=`calc $counter/$lines*100`
	counter=$[$counter + 1]

	progress=`calc $counter/$lines*100`
	echo -ne "`date`\t$progress%\n"

done < $batches_comb

echo
echo "DONE!!"
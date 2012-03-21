#!/bin/sh
path_to_clustal=./src/clustalw2
path_to_source=./src/PfamTest.fasta
path_to_result=./src/PfamTest.aln
path_to_result_dnd=./src/PfamTest.dnd
path_to_sample=./sample/PfamTest.aln

echo "All right script started..."
echo `$path_to_clustal $path_to_source` > log.txt
if ! diff $path_to_sample $path_to_result; then
	echo "ERROR! Files NOT equal!"
	exit 1
else
	echo "Ok. Files equal"
	rm -f $path_to_result
	rm -f $path_to_result_dnd
fi
exit 0

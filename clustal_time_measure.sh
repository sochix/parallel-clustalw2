#!/bin/sh
path_to_clustal=./src/clustalw2
path_to_source=./src/PfamTest.fasta

echo "All right script started..."
start_time=`date +%s`
echo "Current time: " $start_time
echo `$path_to_clustal $path_to_source` > log.txt
echo "Program finished. Worked time: "
end_time=`date +%s`
echo $(($end_time - $start_time))
exit 0

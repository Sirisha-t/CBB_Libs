DATA=../../data/strep_2x.fasta
EXEC=../bin/mm_index

## Running the program
/usr/bin/time -f "time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory %MKB \ncpu %P" $EXEC --seq_type nucl --k 10 --w 3 --target_file $DATA --query_file $DATA

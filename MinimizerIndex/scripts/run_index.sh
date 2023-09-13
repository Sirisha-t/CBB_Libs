DNA_DATA=../../data/strep_2x.fasta
PROT_DATA=../../data/strep_2x.faa
EXEC=../bin/mm_index

## Running the program --DNA overlaps
/usr/bin/time -f "time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory %MKB \ncpu %P" $EXEC --seq_type nucl --k 15 --w 3 --outfile dna.overlaps --target_file $DNA_DATA --query_file $DNA_DATA


## Protein overlaps
/usr/bin/time -f "time result\ncmd:%C\nreal %es\nuser %Us \nsys  %Ss \nmemory %MKB \ncpu %P" $EXEC --seq_type prot --k 10 --w 3 --outfile protein.overlaps --target_file $PROT_DATA --query_file $PROT_DATA

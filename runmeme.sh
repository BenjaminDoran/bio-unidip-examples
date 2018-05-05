#!/bin/bash

# example command
# date && ./runmeme.sh meme /home/biodev/MA0852_2/MA0852.2SEQS1K.fasta -dna -oc /home/biodev/MA0852_2/raw -maxsize 6000000 -mod oops -nmotifs 3 -minw 6 -maxw 30 && date 

docker run --rm -v $(pwd):/home/biodev/ ddiez/meme "$@"



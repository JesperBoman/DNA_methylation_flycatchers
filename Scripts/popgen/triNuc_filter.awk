#u!/usr/bin/awk -f

# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #
# awk -v genome=$GENOME -f triNuc_filter.awk <List of positions> > <List of positions with dinucleotide annotation>
# ==========================================================================================================
# Search reference genome for dinucleotide context (here ancestral "CpG" sites) around focal position
# ==========================================================================================================
# Jesper Boman                      28 oct 2021
# ==========================================================================================================
# Example:
# Use in pipe like this (start with a file of counts for each allele at a position):
# awk -v genome=$GENOME -f triNuc_filter.awk
# # # # # # # # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # # # # # # # # #

{
cmd = "samtools faidx " genome " "$1":"$2-1"-"$2+1" | grep -v '>'"
cmd | getline triNuc 

if(triNuc~/CG/){
        print $0 "\t" "CpG_ans";}
else{print $0 "\t" "Other"}
close(cmd);
}

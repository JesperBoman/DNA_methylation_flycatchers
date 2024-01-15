#!/bin/bash -l


ml bioinfo-tools samtools/1.10



#Part I: Filter fixed differences for biallelic differences

wc -l fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper
#38913

awk '{if($4 ~/,/ || $5 ~/,/){next} else{print $0}}' fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper > fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper_biallel 

#37271

awk 'BEGIN{OFS="\t"}NR==1{print "SNP-ID","Chromosome","Position","Strand","Ref/SNP"} NR>1{print NR, $1, $2, 1,$3 "/" $4}' fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper_biallel > fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper_biallel_SNPsplit_format




fixeddiff="fixed_diff_chr_strict_OC_OP_alexn_vcf_ludoscript_Jesper_biallel_SNPsplit_format"




#Part II: Split the BAM files according to fixed differences

mkdir SNPsplitted_BAMs

declare -a arr=("Heart" "Kidney" "Liver" "Brain" "Testis")

for tissue in "${arr[@]}"
do

bamdir="/crex/proj/sllstore2017033/nobackup/work/jesper/DNA_methylation_flycatchers/methylseq_pipe/${tissue}/bismark_deduplicated"
ls "$bamdir" | grep "HYB" | grep -v "HYB04" | grep -v "v2" > file_list_${tissue}

while IFS= read -r inputbam
do
		perl SNPsplit_v0.3.2/SNPsplit -o ./SNPsplitted_BAMs --snp_file "$fixeddiff" $bamdir/$inputbam &
		if [ $(jobs -p | wc -w) -ge 10 ]; then
	        wait
      		fi	
	done < "file_list_${tissue}"
done
wait

mv $bamdir/*genome* /crex/proj/sllstore2017033/nobackup/work/jesper/DNA_methylation_flycatchers/asm_SNPs/SNPsplitted_BAMs/
mv $bamdir/*SNPsplit* /crex/proj/sllstore2017033/nobackup/work/jesper/DNA_methylation_flycatchers/asm_SNPs/SNPsplitted_BAMs/
mv $bamdir/*allele* /crex/proj/sllstore2017033/nobackup/work/jesper/DNA_methylation_flycatchers/asm_SNPs/SNPsplitted_BAMs/
mv $bamdir/*unassigned* /crex/proj/sllstore2017033/nobackup/work/jesper/DNA_methylation_flycatchers/asm_SNPs/SNPsplitted_BAMs/




#Part III: Extract methylation calls 

declare -a arr=("Heart" "Kidney" "Liver" "Brain" "Testis")

for tissue in "${arr[@]}"
do
ls SNPsplitted_BAMs | grep "genome" | grep "$tissue" > files_tmp

if [ $tissue == "Kidney" ]; then
grep -v "HYB02" files_tmp > files_tmp2
mv files_tmp2 files_tmp
fi

while IFS= read -r inputbam
do

/sw/bioinfo/bismark/0.22.1/rackham/bismark_methylation_extractor --paired-end --no_overlap --comprehensive --cytosine_report --bedGraph --report --gzip --samtools_path /sw/apps/bioinfo/samtools/1.10/rackham/bin/samtools -o Bismark_output_CpG_only --genome_folder /methylseq_pipe/Brain/reference_genome/BismarkIndex --multicore 20 --buffer_size 120G --ignore_r2 2  --ignore_3prime_r2 2 SNPsplitted_BAMs/$inputbam

done < "files_tmp"
done

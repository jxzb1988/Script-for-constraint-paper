#!/usr/bin/bash

dosage=$1

num=$2

gene=$3

num_causal=$4

index=$5

local_time=`date +%s`

R --vanilla --slave --input ${dosage} --num ${num_causal} --output  simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_same_${local_time}_${index} --list  SNPs_${gene}_choose_to_be_causal_${num}_same_${local_time}_${index}  < ./simulate_phenotype_same.R


perl -ne 'chomp;my @array=split;print "$array[1]\t$array[0]\t@array[2..$#array]\n";' simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_same_${local_time}_${index} |perl -ne 'chomp;my @array=split;my $str=join("\t",@array);print "$str\n";'   >adjust_simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_same_${local_time}_${index}

echo "It is Ok, now\n"

R --vanilla --slave --input adjust_simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_same_${local_time}_${index}  --output output_file_${gene}_${num}_same_${local_time}_${index} <./simulation_based_on_ZZZ3.R  >output_${gene}_${num}_same_${local_time}_${index}

perl ./compare_result.pl  simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_same_${local_time}_${index}_LD_r2 output_file_${gene}_${num}_same_${local_time}_${index} output_file_${gene}_${num}_same_${local_time}_${index}_condition  >final_result_${gene}_${num}_same_${local_time}_${index}




rm simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_same_${local_time}_${index}  adjust_simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_same_${local_time}_${index}


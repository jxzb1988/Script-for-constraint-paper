#!/usr/bin/bash

dosage=$1

num=$2

gene=$3

num_causal=$4

index=$5

local_time=`date +%s`

R --vanilla --slave --input ${dosage} --num ${num_causal} --output  simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_two_inverse_${local_time}_${index} --list  SNPs_${gene}_choose_to_be_causal_${num}_two_inverse_${local_time}_${index}  <./simulate_phenotype_two_inverse.R


perl -ne 'chomp;my @array=split;print "$array[1]\t$array[0]\t@array[2..$#array]\n";' simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_two_inverse_${local_time}_${index} |perl -ne 'chomp;my @array=split;my $str=join("\t",@array);print "$str\n";'   >adjust_simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_two_inverse_${local_time}_${index}



R --vanilla --slave --input adjust_simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_two_inverse_${local_time}_${index}  --output output_file_${gene}_${num}_two_inverse_${local_time}_${index} </nv/hp10/bzeng30/scratch/zengbiao/simulation_simulation_with_random_picked_genes/scripts/simulation_based_on_ZZZ3.R  >output_${gene}_${num}_two_inverse_${local_time}_${index}



perl ./compare_result.pl  simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_two_inverse_${local_time}_${index}_LD_r2 output_file_${gene}_${num}_two_inverse_${local_time}_${index} output_file_${gene}_${num}_two_inverse_${local_time}_${index}_condition >final_result_${gene}_${num}_two_inverse_${local_time}_${index}



rm simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_two_inverse_${local_time}_${index} adjust_simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_two_inverse_${local_time}_${index}


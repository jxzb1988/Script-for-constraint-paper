#!/usr/bin/bash

dosage=$1

num=$2

gene=$3

R --vanilla --slave --input ${dosage} --output simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num} --list  SNPs_${gene}_choose_to_be_causal_${num}  < ./simulate_phenotype.R


perl -ne 'chomp;my @array=split;print "$array[1]\t$array[0]\t@array[2..$#array]\n";' simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num} |perl -ne 'chomp;my @array=split;my $str=join("\t",@array);print "$str\n";'   >adjust_simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}



R --vanilla --slave --input adjust_simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}  --output output_file_${gene}_${num} <./simulation_based_on_ZZZ3.R  >output_${gene}_${num} 

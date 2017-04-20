#!/usr/bin/bash

dosage=$1

num=$2

gene=$3

num_causal=$4

index=$5

local_time=`date +%s`

R --vanilla --slave --input ${dosage} --num ${num_causal} --output  simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_same_${local_time}_${index}_two_causal --list  SNPs_${gene}_choose_to_be_causal_${num}_same_${local_time}_${index}_two_causal  --ld_file simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_same_${local_time}_${index}_ld_file_two_causal < ./simulate_phenotype_two.R



rm simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_same_${local_time}_${index}_ld_file_two_causal


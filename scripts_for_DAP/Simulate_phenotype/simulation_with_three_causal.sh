#!/usr/bin/bash

dosage=$1

num=$2

gene=$3

num_causal=$4

index=$5

ld_filter=$6

local_time=`date +%s`

R --vanilla --slave --input ${dosage} --num ${num_causal} --output  simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_same_${local_time}_${index}_one_causal  --list  SNPs_${gene}_choose_to_be_causal_${num}_same_${local_time}_${index}_one_causal  --ld_file simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_same_${local_time}_${index}_ld_file  --ld_filter ${ld_filter}< /nv/hp10/bzeng30/scratch/zengbiao/simulated_data_for_Liu_to_evaluate_multiple_eQTL_detection_methods/simulation_to_evaluate_DAP_method/detect_multiple_eQTL_simulation_with_3_causal_and_requirement_for_LD/simulation_data/scripts/simulate_phenotype_single.R


#rm  simulated_phenotype_dos_genotype_${gene}_with_SNP_name_${num}_same_${local_time}_${index}_ld_file


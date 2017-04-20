#!/usr/bin/perl



my $gene;

my %gene;

open(FHIN,"<$ARGV[0]");


while(<FHIN>)
 {
   chomp;my @array=split;
   $gene{$array[-1]}=$_;
 }
close FHIN;

my @key=keys %gene;

my $value_test=@key;

#print "number of genes are: $value_test\n";

my $value=int(rand($value_test));

my $index=$ARGV[2]."_".$ARGV[3];

#print "the chosen gene is: $key[$value], the index is: $value\n";   

my @array=split(/\s+/,$gene{$key[$value]});



#system "sh  /nv/hp10/bzeng30/data/my_script/condition_analysis/simulation_genotype_and_phenotype/extract_genotype.sh $array[-1] /nv/hp10/bzeng30/data2/microarray_annotation/illumina/illumina_annotation/annotation/analysis/location_for_illumina_probe_genes_final  $array[0]";


if(-e "/nv/hp10/bzeng30/data2/data_script_for_multi_eQTL_paper/simulation/real_genotype/genotype_for_each_genes/$array[-1]/combined_$array[-1]\_genotype_phenotype")
 {
   for(my $i=1;$i<=4;$i++)
    {
      system "sh ./simulation_with_one_causal.sh   /nv/hp10/bzeng30/data2/data_script_for_multi_eQTL_paper/simulation/real_genotype/genotype_for_each_genes/$array[-1]/combined_$array[-1]\_genotype_phenotype  $i  $array[-1]  1  $index
sh   ./simulation_with_two_causal.sh    /nv/hp10/bzeng30/data2/data_script_for_multi_eQTL_paper/simulation/real_genotype/genotype_for_each_genes/$array[-1]/combined_$array[-1]\_genotype_phenotype  $i  $array[-1] 2  $index
sh  ./simulation_with_three_causal.sh    /nv/hp10/bzeng30/data2/data_script_for_multi_eQTL_paper/simulation/real_genotype/genotype_for_each_genes/$array[-1]/combined_$array[-1]\_genotype_phenotype  $i  $array[-1] 3  $index";
    }
#  system "cat final*one*_$index >simulation_one_inverse_$index\_result";
#  system "cat final*same*_$index >simulation_same_inverse_$index\_result";
#  system "cat final*two*_$index >simulation_two_inverse_$index\_result";
#  system "rm *_$index";
#  system "rm *_$index\_*LD_r2";
#  system "rm *_$index\_*condition";
 }

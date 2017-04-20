#!/usr/bin/perl

my $pwd=`pwd`;
chomp $pwd;
#print "Working Dic is $pwd\n";
   my $file=$ARGV[0];
   my $file_pbs=$file;
   $file_pbs=$file.".pbs";
   my $file_log=$file;
   $file_log=$file.".log";
   my $file_sh=$file;
   $file_sh=$file.".sh";



#   my $command= "perl -e \'for(my \$i=1;\$i<=150;\$i++){system \"perl  \.\/simulate_random_picked_gene.pl gene_location_for_expressed_genes  4 $ARGV[0] \$i\";\}\' ; cat simulation_one_inverse_$ARGV[0]\_?*_result  >simulation_one_inverse_$ARGV[0]\_result; rm simulation_one_inverse_$ARGV[0]\_?*_result; cat simulation_one_inverse_$ARGV[0]\_?*_result  >simulation_one_inverse_$ARGV[0]\_result; rm simulation_one_inverse_$ARGV[0]\_?*_result; cat simulation_two_inverse_$ARGV[0]\_?*_result  >simulation_two_inverse_$ARGV[0]\_result; rm simulation_two_inverse_$ARGV[0]\_?*_result; cat simulation_same_inverse_$ARGV[0]\_?*_result  >simulation_same_inverse_$ARGV[0]\_result; rm simulation_same_inverse_$ARGV[0]\_?*_result";

if(!(-d $ARGV[0]))
 {
   system "mkdir $ARGV[0]";
 }

    my $command= "perl -e \'for(my \$i=1;\$i<=10;\$i++){system \"perl  ./simulate_random_picked_gene.pl /nv/hp10/bzeng30/data2/simulation_with_random_picked_genes_to_evaluate_influence_on_CLPP/gene_location_for_expressed_genes  4 $ARGV[0] \$i\";\}\' ; cat simulated_phenotype_dos_genotype_*_one_causal_CLPP_causl_combined_with_causal_infor >simulated_phenotype_dos_genotype_all_one_causal_CLPP_causl_combined_with_causal_infor; cat cat simulated_phenotype_dos_genotype_*_two_causal_CLPP_causl_combined_with_causal_infor >simulated_phenotype_dos_genotype_all_two_causal_CLPP_causl_combined_with_causal_infor; cat simulated_phenotype_dos_genotype_*_three_causal_CLPP_causl_combined_with_causal_infor >simulated_phenotype_dos_genotype_all_three_causal_CLPP_causl_combined_with_causal_infor; rm simulated_phenotype_dos_genotype_*_with_SNP_name_*_causal_CLPP_causl_combined_with_causal_infor";



   open(FHOU,"+>$pwd/$file/$file_pbs");
   print FHOU "#PBS -N order gene
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=12:00:00
#PBS -q iw-shared-6
#PBS -j oe
#PBS -o $file_log
cd $pwd/$file
sh $file_sh";
close FHOU;
  open(FHOU,"+>$pwd/$file/$file_sh");
  print FHOU "$command";
  close FHOU;
  system "msub $pwd/$file/$file_pbs";
#  print "command is $command\n";

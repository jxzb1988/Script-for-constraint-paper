#!/usr/bin/perl

my $pwd=`pwd`;
   my $file=$ARGV[0]."/".$ARGV[0];
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

    my $command= "perl -e \'for(my \$i=1;\$i<=30;\$i++){system \"perl  \.\.\/simulate_random_picked_gene.pl ../gene_location_for_expressed_genes  4 $ARGV[0] \$i\";\}\' ; cat simulation_one_inverse_$ARGV[0]\_?*_result  >simulation_one_inverse_$ARGV[0]\_result; rm simulation_one_inverse_$ARGV[0]\_?*_result; cat simulation_same_inverse_$ARGV[0]\_?*_result  >simulation_same_inverse_$ARGV[0]\_result; rm simulation_same_inverse_$ARGV[0]\_?*_result; cat simulation_two_inverse_$ARGV[0]\_?*_result  >simulation_two_inverse_$ARGV[0]\_result; rm simulation_two_inverse_$ARGV[0]\_?*_result;";



   open(FHOU,"+>$file_pbs");
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
  open(FHOU,"+>$file_sh");
  print FHOU "$command";
  close FHOU;
#  system "msub $file_pbs";
#  print "command is $command\n";

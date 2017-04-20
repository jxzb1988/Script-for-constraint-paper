#!/usr/bin/perl


my $folder=$ARGV[0];

#my $folder=$ARGV[1];


if($folder1=~/\/$/)
 {
   $folder1=~s/\/$//;
 }


if($folder=~/\/$/)
 {
   $folder=~s/\/$//;
 }

my @file=glob("$folder/[1-9]*/simulated_phenotype_dos_genotype_*");

#print @file;

#exit;

my $pwd=`pwd`;
foreach my $file1 (@file)
 {
   my @test=split(/\//,$file1);
   my $file=$test[-1];
   my $file_pbs=$file;
   $file_pbs=$file.".pbs";
   my $file_log=$file;
   $file_log=$file.".log";
   my $file_sh=$file;
   $file_sh=$file.".sh";
   my $command= "sh ./perform_DAP.sh  $file1 ";
   open(FHOU,"+>$file_pbs");
   print FHOU "#PBS -N order gene
#PBS -l nodes=1:ppn=2
#PBS -l mem=10gb
#PBS -l walltime=12:00:00
#PBS -q iw-shared-6
#PBS -j oe
#PBS -o $file_log
cd $pwd
sh $file_sh";
close FHOU;
  open(FHOU,"+>$file_sh");
  print FHOU "$command";
  close FHOU;
  system "msub $file_pbs";
#  print "command is $command\n";
 }

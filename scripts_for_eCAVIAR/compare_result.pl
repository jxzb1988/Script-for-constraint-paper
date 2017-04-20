#!/usr/bin/perl


my %snp;

open(FHIN,"<$ARGV[0]");
<FHIN>;
while(<FHIN>)
 {
   chomp;my @array=split;
   $snp{$array[0]}=$_;
 }
close FHIN;

my %result1;

my %result2;

open(FHIN,"<$ARGV[1]");

while(<FHIN>)
 {
   if(/^\s/)
    {
      next;
    }
   chomp;my @array=split;
   $result1{$array[1]}=$array[3];
   
 }
close FHIN;

open(FHIN,"<$ARGV[2]");

while(<FHIN>)
 {
   if(/^\s/)
    {
      next;
    }
   chomp;my @array=split;
   $result2{$array[0]}=$array[1];

 }
close FHIN;



foreach my $value (keys %snp)
 {
   my $diff1;
   my $diff2;
   if(exists($result1{$value}))
    {
      my @arr1=split(/\s+/,$snp{$value});
      my @arr2=split(/\s+/,$result1{$value});
      $diff1=abs($arr1[2]-$arr2[0]);
    } else
    {
      my @arr1=split(/\s+/,$snp{$value});
      $diff1=abs($arr1[2]);
    }
   if(exists($result2{$value}))
    {
      my @arr1=split(/\s+/,$snp{$value});
      my @arr2=split(/\s+/,$result2{$value});
      $diff2=abs($arr1[2]-$arr2[0]);
    } else
    {
      my @arr1=split(/\s+/,$snp{$value});
      $diff2=abs($arr1[2]);
    }
   print "$snp{$value}\t$diff1\t$diff2\n";
 }

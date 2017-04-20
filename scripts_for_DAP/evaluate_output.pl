#!/usr/bin/perl


my $causal_file=$ARGV[0];

my $output_DAP=$ARGV[1];

my %summary;

my %causal;
my %result_DAP;

open(FHIN,"<$causal_file");
<FHIN>;
while(<FHIN>)
 {
   chomp;my @array=split;
   $causal{$array[0]}=1;
 }

close FHIN;

my $result_num=0;

my $candidate=0;

$detected_causal=0;

open(FHIN,"<$output_DAP");

#foreach my $value (keys %causal)
# {
#   print "causal variants are: $value\n";
# }

while(<FHIN>)
 {
   chomp;my @array=split;
   if(/^Posterior expected model size/)
    {
      $result_num=$array[4];
      next;
    }
   if($array[0]=~/[1-9]*/ && @array==4)
    {
      $result_DAP{$array[1]}=$_;
      $candidate++;
 #     print "SNP is: $array[1]\n";
      if(exists($causal{$array[1]}))
       {
     #    print "SNP is: $array[1]\n";
         $detected_causal++;
       }
    }
 }

close FHIN;


print "Number of SNPs in confidence interval is: $candidate, detected causal variants is: $detected_causal, and the estimated causal variant number is $result_num\n";

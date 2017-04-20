#!/usr/bin/sh


file=$1

arrIN=($(echo ${file} | tr "/" "\n"))

#echo ${arrIN[-1]}



output=${arrIN[${#arrIN[@]}-1]}


perl ~/data/my_script/convert_expr-file_to_plink-pheno.pl  ${file} |awk 'NR>1' |perl -e 'my $num=0;while(<>){$num++;if($num==1){chomp;my @array=split;print "pheno phe CAGE";my $string=join(" ",@array[1..$#array]);print " $string\n";} else{chomp;my @array=split;print "geno $array[0] CAGE";my $string=join(" ",@array[1..$#array]);print " $string\n";}}'   >${output}

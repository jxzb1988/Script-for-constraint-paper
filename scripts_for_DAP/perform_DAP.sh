

file=$1

arrIN=($(echo ${file} | tr "/" "\n"))

#echo ${arrIN[-1]}



output=${arrIN[${#arrIN[@]}-1]}




sh ./convert_to_sbam_format.sh   ${file}

/nv/hp10/bzeng30/data2/software/dap-master/dap_src/dap   -d  ${output}  -t 2 -msize 4 -g /nv/hp10/bzeng30/data2/software/dap-master/dap_src/sample_data/sample.fixed_effect.grid   >${output}_DAP

#/nv/hp10/bzeng30/data2/software/dap-master/dap_src/dap   -d  ${output} -msize  2  -g /nv/hp10/bzeng30/data2/software/dap-master/dap_src/sample_data/sample.fixed_effect.grid   >${output}_DAP



snp_file=($(echo ${output} | tr "_" "\n"))

prefix=""

if [[ $file =~ "no_enrichment" ]];then
 prefix=${snp_file[10]}"_"${snp_file[11]}"_"${snp_file[12]}"*_no_enrichment*"
else
 prefix=${snp_file[10]}"_"${snp_file[11]}"_"${snp_file[12]}
fi



#prefix=${snp_file[10]}"_"${snp_file[11]}"_"${snp_file[12]}


dir=$(dirname "${file}")

causal_file=${dir}"/""SNPs_*"${prefix}"*"

final_causal_file=$(echo ${causal_file})
echo "Causal file is:"
echo ${final_causal_file}

perl ./evaluate_output.pl  ${final_causal_file}  ${output}_DAP   >${output}_DAP_summary_result



#rm ${output}  ${output}_DAP  

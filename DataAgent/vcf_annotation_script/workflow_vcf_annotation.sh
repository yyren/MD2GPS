#!/bin/bash
input_vcf=$1
out_file=$2
Docker_Image=$3
REF_PATH=$4


bed_name_with_path='/home/ubuntu/0'
save_bed_with_path='/home/ubuntu/0'
fileType='0'
panelId='29'
pipeline='Mendelian-Disease'
psid='bmap'
race='eas'
single_or_compare='single'

DNA_panel_path=`dirname $0`

echo "single path:$DNA_panel_path"
source $DNA_panel_path/DNA_panel_all_pid_parameter.sh $Docker_Image
source $DNA_panel_path/main_scripts/annotation.sh

file_with_path=`realpath $input_vcf`
work_path=${file_with_path%/*}

cd $work_path

if [[ ! -d $work_path/log ]]; then
    mkdir $work_path/log
else
    rm -r $work_path/log
    mkdir $work_path/log
fi
switch_vcf_to_standard_snpEff single $input_vcf
singularity exec $Docker_Image perl $DNA_panel_path/main_scripts/11_remove_duplicates.pl tumor_for_snpEff_adjust.vcf ${psid}_merged_for_snpEff.vcf


#left or right most occurance
mv ${psid}_merged_for_snpEff.vcf ${psid}_merged_for_snpEff_need_changed.vcf 
left_or_right_align ${psid}_merged_for_snpEff_need_changed.vcf ${psid}_merged_for_snpEff.vcf

#variants annotation
snpEff_annovar_annotation ${psid}_merged_for_snpEff.vcf $psid $save_bed_with_path $race $max_insertion_length ${psid}_annotated_results.txt $fileType $only_exon_splice $single_or_compare $panelId $pipeline
## search chemotherapy and prognosis related snp
#Add_snp_annotation $fileType $psid
filter_by_MAF=no
filter_by_snp=no
classification ${psid}_annotated_results.txt results_itd.txt $psid $bed_name_with_path $panelId $pipeline $race $out_file


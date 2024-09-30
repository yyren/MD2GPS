#!/bin/bash

check_error()  {
  if [ $1 -ne 0 ]; then
    echo "ERROR CODE: $2"
    exit 1
  fi
}
check_success()  {
	if [ -f sucess.txt ]; then
		rm sucess.txt
	fi
	if [ -f results_all.txt ]; then
		results_line_numb=`cat results_all.txt |wc -l`
		if [[ $results_line_numb -eq 0 || $results_line_numb -eq 1 ]]; then
			echo "Empty result! Please check the selected analysis options (eg. Select NGS Panel option)" >> success.txt
		fi
	fi
}

# $1: filetype; $2: customized filter(yes/no); $3: input mpileup file; $4: output the setting coverage file with input mpileup file name; $5: output mpileup file of target region; $6: bed name with path;
set_total_cov(){
	if [ $2 = "yes" ]; then
		total_cov_setted=$total_cov
		echo "$3	$total_cov_setted	$total_cov_setted" > $4
	else
		if [[ $1 != "0" && $1 != "4" ]]; then
			if [[ -f $6 && $2 = "no" ]]; then
				$docker_run_default perl $DNA_panel_path/main_scripts/21_target_off_target_readdepth.pl $6 $3 $5
				$docker_run_default Rscript $DNA_panel_path/main_scripts/22_statistic_readdepth_target.r $5 $4
			else
				$docker_run_default Rscript $DNA_panel_path/main_scripts/22_statistic_readdepth_all.r $total_cov $3 $4
			fi
			while read line
			do
				rd_cf=${line#*\t}
				break
			done < $4
			c_arr=($rd_cf)
			total_cov_setted=${c_arr[1]}
		else
			total_cov_setted=$total_cov_vcf
			echo "$3	$total_cov_setted	$total_cov_setted" > $4
		fi
	fi

}

set_var_cov(){
	if [[ $2 != "0" && $2 != "4" ]]; then
		if [ $1 = "Illumina" ]; then
			var_cov_snp_setted=$varcov_snp_illumina
			var_cov_indel_setted=$varcov_indel_illumina
			var_freq_setted=$varfreq_illumina
		elif [ $1 = "Ion_torrent" ]; then
			var_cov_snp_setted=$varcov_snp_Ion_torrent
			var_cov_indel_setted=$varcov_indel_Ion_torrent
			var_freq_setted=$varfreq_Ion_torrent
		fi
	else
		var_cov_snp_setted=$varcov_snp_vcf
		var_cov_indel_setted=$varcov_indel_vcf
		var_freq_setted=$varfreq_vcf
	fi

}

# $1: fileType; $2: single, multi tumor sample or tumor-normal compare sample; $3: input setting coverage file; $4: input setting support variant reads for snp
# $5: input setting support variant reads for indel; $6: normal or cancer(normal/cancer)
# output: filter_clinical_cn1.txt, filter_clinical_en1.txt
generate_filters(){
	echo "generate filter"
	$docker_run_default perl $DNA_panel_path/main_scripts/35_filter_auto.pl $1 $2 $3 $4 $5 $6
	check_error $? generate_filter
}

#prefilter for raw vcf file
#$1:single or compare; $2:input vcf file; $3: tumor psid; $4: normal psid
switch_vcf_to_standard_snpEff(){
	if [ $1 = "single" ]; then
		$docker_run_default perl $DNA_panel_path/main_scripts/11_vcf_switch.pl $2
        $docker_run_default perl $DNA_panel_path/main_scripts/11_vcf_standard_mnp.pl tumor_for_snpEff.vcf tumor_for_snpEff_mnp.vcf
		$docker_run_default perl $DNA_panel_path/main_scripts/12_change_merge_format.pl tumor_for_snpEff_mnp.vcf tumor_for_snpEff_adjust.vcf $REF_PATH/b37/human_g1k_v37_decoy.fasta
	else
		$docker_run_default perl $DNA_panel_path/main_scripts/11_vcf_switch.pl $2 $3 $4
        $docker_run_default perl $DNA_panel_path/main_scripts/11_vcf_standard_mnp.pl tumor_for_snpEff.vcf tumor_for_snpEff_mnp.vcf
        $docker_run_default perl $DNA_panel_path/main_scripts/11_vcf_standard_mnp.pl normal_for_snpEff.vcf normal_for_snpEff_mnp.vcf
        $docker_run_default perl $DNA_panel_path/main_scripts/11_vcf_standard_mnp.pl compare_for_snpEff.vcf compare_for_snpEff_mnp.vcf
        
		$docker_run_default perl $DNA_panel_path/main_scripts/12_change_merge_format.pl tumor_for_snpEff_mnp.vcf tumor_for_snpEff_adjust.vcf $REF_PATH/b37/human_g1k_v37_decoy.fasta
		$docker_run_default perl $DNA_panel_path/main_scripts/12_change_merge_format.pl normal_for_snpEff_mnp.vcf normal_for_snpEff_adjust.vcf $REF_PATH/b37/human_g1k_v37_decoy.fasta
		$docker_run_default perl $DNA_panel_path/main_scripts/12_change_merge_format.pl compare_for_snpEff_mnp.vcf compare_for_snpEff_adjust.vcf $REF_PATH/b37/human_g1k_v37_decoy.fasta
	fi


}

#extract somatic from normal and tumor vcf
#$1: normal for snpeff file; $2: tumor for snpeff file; $3: out somatic file
combine_somatic_for_snpEff(){
	$docker_run_default perl $DNA_panel_path/main_scripts/13_combine_normal_tumor_for_snpEff.pl $1 $2 $3
}

# $1: the source of file1(e.g varscan); $2: file1 with full path; $3: the source of file2(e.g mutect); $4: file2 with full path ; $5: merged file according to the value(chromsome, start position,ref,alt)
merge_same_data(){
	if [ -f $4 ]; then
		$docker_run_default perl $DNA_panel_path/main_scripts/11_merge_diff_source_for_snpEff.pl $1 $2 $3 $4 $5
		check_error $? merge_diff_source
	fi
}

#right or left align(Ion torrent(bwa,tmap): -:pass; +:occour at most right coordinate. Illumina(stampy): +: pass; -:occour at most left coordinate)
# old:$1:bed file with the tag of strand; $2: input file need to be changed; $3: platform(Illumina or Ion_torrent); $4: reference with path; $5:output file
#using now: $1: input file; $2: bed with strand information; $3: reference; $4: output file
left_or_right_align(){
	##perl $DNA_panel_path/main_scripts/12_standardize_according_to_HGVS.pl $1 $DNA_panel_path/panel_with_strand/GRCh37_all_gene_region.bed $REF_PATH/b37/human_g1k_v37_decoy.fasta snpEff_inputfile_hgvs_format.vcf
	##perl $DNA_panel_path/main_scripts/05_remove_duplicates.pl snpEff_inputfile_hgvs_format.vcf $2
	$docker_run_default perl $DNA_panel_path/main_scripts/05_remove_duplicates.pl $1 $2
	#perl $DNA_panel_path/main_scripts/12_add_strand_refseq.pl $1 $2 $3 $4 snpEff_inputfile_add_strand_seq.vcf
	#perl $DNA_panel_path/main_scripts/12_right_or_left_align.pl snpEff_inputfile_add_strand_seq.vcf $3 $5
}

# $1: software_index;
change_index_to_name(){
	if [ $1 = "1" ]; then
		software_name="tvc"
		file_numb=2
		input_file1="tvc_for_snpEff.vcf"
		input_file2="tvc_indel_for_snpEff.vcf"
	elif [ $1 = "2" ]; then
		software_name="mutect"
		file_numb=1
		input_file1="mutect_for_snpEff.vcf"
	elif [ $1 = "4" ]; then
		software_name="pindel"
		file_numb=2
		input_file1="pindel_SI_for_snpEff.vcf"
		input_file2="pindel_del_for_snpEff.vcf"
	elif [ $1 = "8" ]; then
		software_name="platypus"
		file_numb=1
		input_file1="platypus_for_snpEff.vcf"
	elif [ $1 = "16" ]; then
		software_name="varscan"
		file_numb=2
		input_file1="varscan_snp_for_snpEff.vcf"
		input_file2="varscan_indel_for_snpEff.vcf"
	elif [ $1 = "32" ]; then
		software_name="gatk"
		file_numb=3
		input_file1="gatk_snp_for_snpEff.vcf"
		input_file2="gatk_indel_for_snpEff.vcf"
		input_file3="gatk_mixed_for_snpEff.vcf"
	fi


}

# $1:software_strategy; $2: psid; $3: all_software index in array
merge_all(){
software_strategy=$1
psid=$2
all_software=$3
merged_file=$4
array_idx=0
select_numb=0
echo $3
for software_index in ${all_software[*]}
do
	change_index_to_name $software_index
	software_name_temp=$software_name
	echo "s1:$software_index|$software_name_temp|$file_numb"
	if [[ $[software_strategy&software_index] = "$software_index" ]]; then #traverse all selected software in pipeline
		echo "$select_numb"
		if [ $select_numb = "0" ]; then #the first appearance in the selected variant callers
			if [ $file_numb = "1" ]; then #if there is only one output file from the selected software
				echo "s2 merged:$software_name_temp|$file_numb"
				touch temp_for_merge.vcf
				merge_same_data $software_name_temp ${psid}_${input_file1} temp temp_for_merge.vcf ${psid}_merged_${software_name_temp}_for_snpEff.vcf
				cp ${psid}_merged_${software_name_temp}_for_snpEff.vcf ${psid}_merged_for_snpEff.vcf
			elif [ $file_numb = "2" ]; then # if there are two output files from the selected software
				echo "s3:$software_name_temp|$file_numb"
				merge_same_data $software_name_temp ${psid}_${input_file1} $software_name_temp ${psid}_${input_file2} ${psid}_merged_${software_name_temp}_for_snpEff.vcf
				cp ${psid}_merged_${software_name_temp}_for_snpEff.vcf $merged_file
			fi
		elif [[ $select_numb -gt 0 && $((software_strategy&software_index)) = "$software_index" ]]; then #not the first appearance in the selected variant callers
			echo "s4:$select_numb|$software_index|"
			for software_index_less in ${all_software[*]:$array_idx}
			do
				if [[ $[software_strategy&software_index_less] = "$software_index_less" ]]; then
					change_index_to_name $software_index_less
					software_name_temp2=$software_name
					echo "s5:$software_name_temp2|$file_numb"
					if [ $file_numb = "1" ]; then
						merge_same_data temp ${psid}_merged_for_snpEff.vcf $software_name_temp2 ${psid}_${input_file1} ${psid}_merged_${software_name_temp2}_for_snpEff.vcf
						cp ${psid}_merged_${software_name_temp2}_for_snpEff.vcf $merged_file
					elif [ $file_numb = "2" ]; then
						merge_same_data $software_name_temp2 ${psid}_${input_file1} $software_name_temp2 ${psid}_${input_file2} ${psid}_merged_${software_name_temp2}_for_snpEff.vcf
						merge_same_data temp ${psid}_merged_for_snpEff.vcf $software_name_temp2 ${psid}_merged_${software_name_temp2}_for_snpEff.vcf ${psid}_merged_${software_name_temp}_${software_name_temp2}_for_snpEff.vcf
						cp ${psid}_merged_${software_name_temp}_${software_name_temp2}_for_snpEff.vcf $merged_file
					fi
				fi
			done
			break 1
		fi
		select_numb=$[ select_numb + 1 ]
	fi
	array_idx=$[ numb + 1 ]
done

}


#$1: input file; $2: output file; $3: panelID; $4: idx
snpeff_annotation(){
    echo "$1|$2|$3|$4"
    idx=
    if [[ -z $4 ]]; then
        idx=1
    else
        idx=$4
    fi
	if [[ $3 -eq 38 ]]; then
		$docker_run_default java -jar /home/bmap/software/snpeff/snpEff/snpEff.jar -noNextProt -v GRCh37.75 -formatEff -onlyTr /home/ubuntu/ngs_script/DNA_panel/qc/panel_bed/38/Transcript.txt -onlyProtein $1 > $2
	else
        $docker_run_default perl $DNA_panel_path/main_scripts/built_liftover_input_bed.pl $1 ${idx}_hg19_liftover_input_all.bed
        $docker_run_default /home/bmap/software/liftover/liftOver ${idx}_hg19_liftover_input_all.bed /home/bmap/software/liftover/hg19ToHg38.over.chain.gz ${idx}_hg38_input_part.bed ${idx}_unMapped.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/divide_hg19_hg38_records.pl $1 ${idx}_hg38_input_part.bed ${idx}_merged_for_snpEff_hg19.vcf ${idx}_merged_for_snpEff_hg38.vcf
        #$docker_run_default java -Xmx10g -jar /home/bmap/software/snpeff/snpEff/snpEff.jar -v GRCh37.75 -noNextProt -formatEff -canon -onlyProtein ${idx}_merged_for_snpEff_hg19.vcf > ${idx}_snpEff_annotated_hg19.vcf
		$docker_run_default java -Xmx10g -jar /home/bmap/software/snpeff/snpEff/snpEff.jar -v GRCh37.75 -noNextProt -formatEff -canon -onlyProtein $1 > ${idx}_snpEff_annotated_hg19.vcf
        $docker_run_default java -Xmx10g -jar /home/bmap/software/snpeff/snpEff/snpEff.jar -v GRCh38.86 -noNextProt -formatEff -canon -onlyProtein ${idx}_merged_for_snpEff_hg38.vcf > ${idx}_snpEff_annotated_hg38.vcf
        
		$docker_run_default perl $DNA_panel_path/main_scripts/16_clean_by_bed_for_vcf.pl $DNA_panel_path/qc/panel_bed/snpeff_cano/hg38_pos_bed.txt ${idx}_merged_for_snpEff_hg38.vcf ${idx}_hg38_in_given_transcript_region.vcf
		$docker_run_default java -Xmx10g -jar /home/bmap/software/snpeff/snpEff/snpEff.jar -noNextProt -v GRCh38.86 -formatEff -onlyTr $DNA_panel_path/qc/panel_bed/snpeff_cano/Transcript.txt -onlyProtein ${idx}_hg38_in_given_transcript_region.vcf > ${idx}_snpEff_annotated_hg38_in_given_transcript_region.vcf
        $docker_run_default perl $DNA_panel_path/main_scripts/merge_hg19_hg38_snpeff_anno.pl ${idx}_snpEff_annotated_hg19.vcf ${idx}_snpEff_annotated_hg38_in_given_transcript_region.vcf ${idx}_snpEff_annotated_hg38.vcf $2
	fi
}

#$1: input file; $2: patientID_sampleID($psid); $3 save bed with path; $4 race(e.g eas); $5: maximum length of insertion; $6: output file (merge snpeff, annovar and add ci and p-value); $7: filetype; $8: only exon or splice site (yes/no); $9: single or compare or multi; ${10}:panleID; ${11}:pipeline
#function: snpEff and annovar annotation; add p-value and CI; remove the insertion variants whose length is bigger than $4; change the ref characters to the number of the length if the length of ref is bigger than 1000
snpEff_annovar_annotation(){
	echo "$1|$2|$3|$4|$5|$6"
	annovar_database_path=$REF_PATH/annovar-database
	snpeff_annotation $1 "$2"_snpEff_annotated.vcf ${10} 1 >> ./log/log_only_snpEff_1.txt 2>&1
	if [[ -f $3 && -f "$2"_snpEff_annotated.vcf ]]; then
		$docker_run_default perl $DNA_panel_path/main_scripts/12_change_snpEff_format.pl "$2"_snpEff_annotated.vcf "$2"_snpEff_annotated_change_format.vcf $9
		$docker_run_default perl $DNA_panel_path/main_scripts/12_extract_mutation_in_save_bed.pl $3 "$2"_snpEff_annotated_change_format.vcf "$2"_snpEff_clean_in_save_bed.vcf
		if [[ $8 = "no" || ${11} = "Pharma-Genom" ]]; then
			cat "$2"_snpEff_annotated_change_format.vcf "$2"_snpEff_clean_in_save_bed.vcf > "$2"_snpEff_annotated_clean.vcf
		elif [[ $8 = "yes" ]]; then
			$docker_run_default perl $DNA_panel_path/main_scripts/12_clean_filter_snpEff_annotated.pl "$2"_snpEff_annotated.vcf "$2"_snpEff_annotated_without_save_bed.vcf $9
			cat "$2"_snpEff_annotated_without_save_bed.vcf "$2"_snpEff_clean_in_save_bed.vcf > "$2"_snpEff_annotated_clean.vcf
		fi
	elif [ -f "$2"_snpEff_annotated.vcf ]; then
		echo "clean results and not exist save_bed_name" >> step_log.txt
		if [[ $8 = "no" || ${11} = "Pharma-Genom" ]]; then
			$docker_run_default perl $DNA_panel_path/main_scripts/12_change_snpEff_format.pl "$2"_snpEff_annotated.vcf "$2"_snpEff_annotated_clean.vcf $9
		elif [[ $8 = "yes" ]]; then
			$docker_run_default perl $DNA_panel_path/main_scripts/12_clean_filter_snpEff_annotated.pl "$2"_snpEff_annotated.vcf "$2"_snpEff_annotated_clean.vcf $9
		fi
	else
		echo "not exist ${2}_snpEff_annotated.vcf" >> error_log.txt
	fi
	$docker_run_default perl $DNA_panel_path/main_scripts/12_extract_null_exon.pl "$2"_snpEff_annotated_clean.vcf "$2"_snpEff_null_exon.vcf "$2"_snpEff_null_exon_for_snpEff.vcf
	snpeff_annotation "$2"_snpEff_null_exon_for_snpEff.vcf "$2"_snpEff_null_exon_annotated.vcf ${10} 2 >> ./log/log_only_snpEff_2.txt 2>&1
	$docker_run_default perl $DNA_panel_path/main_scripts/12_combine_null_exon.pl "$2"_snpEff_null_exon.vcf "$2"_snpEff_null_exon_annotated.vcf "$2"_snpEff_annotated_clean.vcf "$2"_snpEff_annotated_clean_combine_exon.vcf
	$docker_run_default perl $DNA_panel_path/main_scripts/13_change_format_for_annovar_with_ancor.pl "$2"_snpEff_annotated_clean_combine_exon.vcf "$2"_snpEff_result_for_annovar.txt
	$docker_run_default perl /home/bmap/software/annovar/annovar/table_annovar.pl "$2"_snpEff_result_for_annovar.txt $annovar_database_path -buildver hg19 -out "$2"_annovar_annotated -remove -protocol gnomad_exome,ljb26_all,cosmic74,clinvar_20150330,avsnp142 -operation f,f,f,f,f -nastring NA
    #perl $ngs_software_path/annovar/table_annovar.pl "$2"_snpEff_result_for_annovar.txt $annovar_database_path -buildver hg38 -out "$2"_annovar_annotated -remove -protocol gnomad_exome,ljb26_all,cosmic70,clinvar_20190305,avsnp150 -operation f,f,f,f,f -nastring NA
	$docker_run_default perl $DNA_panel_path/main_scripts/13_merge_snpEff_annovar_with_ancor.pl "$2"_snpEff_annotated_clean_combine_exon.vcf "$2"_annovar_annotated.hg19_multianno.txt "$2"_merge_snpEff_annovar.txt
	$docker_run_default perl $DNA_panel_path/main_scripts/14_rm_long_indel.pl "$2"_merge_snpEff_annovar.txt "$2"_merge_snpEff_annovar_rm_long_indel.txt $5
	$docker_run_default Rscript $DNA_panel_path/main_scripts/15_add_p_value_ci.r "$2"_merge_snpEff_annovar_rm_long_indel.txt $6 $7

}


#$1: input file; $2: input itd result(results_itd.txt); $3: patientID_sampleID($psid); $4: bed name; $5: panel id; $6:pipeline; $7:race; $8: out file
classification(){
	if [ $filter_by_MAF = "yes" ]; then
		MAF_index=8
	else
		MAF_index=0
	fi
	if [[ $run_itd_analysis = "yes" && -f $2 ]]; then
		itd_index=4
	else
		itd_index=0
	fi
	if [ $filter_by_snp = "yes" ]; then
		snp_index=2
	else
		snp_index=0
	fi
	if [[ $filter_by_bed = "yes" && -f $4 ]]; then
		bed_index=1
	else
		bed_index=0
	fi
	class_idx=$[ itd_index + snp_index + bed_index ]
	if [ $class_idx -eq 0 ]; then
		echo "no filter" >> step_log.txt
		cp $1 "$3"_result_filtered.txt
	elif [ $class_idx -eq 1 ]; then
		echo "clean by bed" >> step_log.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/16_clean_off_target.pl $4 $1 "$3"_result_filtered.txt
		check_error $? classification_by_bed
	elif [ $class_idx -eq 2 ]; then
		echo "filtered common snp" >> step_log.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/16_filter_snp.pl $1 "$3"_result_filtered.txt
		check_error $? classification_snp
	elif [ $class_idx -eq 4 ];then
		echo "remove duplicate from ITD" >> step_log.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/16_rm_dup_from_itd.pl $1 $2 "$3"_result_filtered.txt
		check_error $? classification_ITD
	elif [ $class_idx -eq 3 ]; then
		echo "clean by bed and filtered common snp " >> step_log.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/16_clean_off_target.pl $4 $1 "$3"_clean_by_bed.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/16_filter_snp.pl "$3"_clean_by_bed.txt "$3"_result_filtered.txt
		check_error $? classification_bed_snp
	elif [ $class_idx -eq 5 ]; then
		echo "clean by bed and remove duplicate from ITD" >> step_log.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/16_clean_off_target.pl $4 $1 "$3"_clean_by_bed.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/16_rm_dup_from_itd.pl "$3"_clean_by_bed.txt $2 "$3"_result_filtered.txt
		check_error $? classification_bed_ITD
	elif [ $class_idx -eq 6 ]; then
		echo "filtered common snp and remove duplicate from ITD" >> step_log.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/16_filter_snp.pl $1 "$3"_clean_by_snp.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/16_rm_dup_from_itd.pl "$3"_clean_by_snp.txt $2 "$3"_result_filtered.txt
		check_error $? classification_snp_ITD
	elif [ $class_idx -eq 7 ]; then
		echo "clean by bed and filtered common snp and remove duplicate from ITD" >> step_log.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/16_clean_off_target.pl $4 $1 "$3"_clean_by_bed.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/16_filter_snp.pl "$3"_clean_by_bed.txt "$3"_clean_by_bed_snp.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/16_rm_dup_from_itd.pl "$3"_clean_by_bed_snp.txt $2 "$3"_result_filtered.txt
		check_error $? classification_bed_snp_itd
	fi
	class_idx=$[ bed_index + MAF_index ]
	if [ $class_idx -eq 0 ]; then
		echo "do not filter by bed and MAF" >> step_log.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/20_only_add_tag.pl "$3"_result_filtered.txt results_all.txt
	elif [ $class_idx -eq 1 ]; then
		echo "filtered by bed but do not filter by MAF" >> step_log.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/20_only_add_tag.pl "$3"_result_filtered.txt results_all.txt
	elif [ $class_idx -eq 8 ]; then
		echo "filter by MAF but do not filter by bed" >> step_log.txt
		$docker_run_default perl $DNA_panel_path/main_scripts/18_classification_by_MAF_panel.pl "$3"_result_filtered.txt results_all.txt $7
	elif [ $class_idx -eq 9 ]; then
		echo "filter by bed and MAF" >> step_log.txt
		if [[ $5 -gt 0 && $5 -lt 10 ]]; then
			echo "filter by bed and hotsport:$5" >> step_log.txt
			$docker_run_default perl $DNA_panel_path/main_scripts/19_classification_by_hotsport_panel.pl $DNA_panel_path/main_scripts/"$5"_panel_target.txt "$3"_result_filtered.txt results_all.txt $7
		else
			echo "filter by bed and MAF" >> step_log.txt
			$docker_run_default perl $DNA_panel_path/main_scripts/18_classification_by_MAF_panel.pl "$3"_result_filtered.txt results_all.txt $7
		fi
	fi
	mv results_all.txt results_all_temp.txt
    #$docker_run_default perl $DNA_panel_path/main_scripts/filter_genetic_variants.pl results_all_temp.txt results_all_filtered_temp.txt
    #$docker_run_default perl $DNA_panel_path/main_scripts/24_add_ancorbase_change_pos_MNP_for_MutationTaster.pl $REF_PATH/b37/human_g1k_v37_decoy.fasta results_all_filtered_temp.txt $8
	$docker_run_default perl $DNA_panel_path/main_scripts/24_add_ancorbase_change_pos_MNP_for_MutationTaster.pl $REF_PATH/b37/human_g1k_v37_decoy.fasta results_all_temp.txt $8
}

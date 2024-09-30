#!/bin/bash
Docker_Image=$1

config_path=`dirname $0`
echo "path:$config_path"
docker_run_default="singularity exec $Docker_Image"

DNA_panel_config_path=$config_path
echo "2customized config:$DNA_panel_config_path"
for((i=1;i<4;i++));
do
	if [ ! -f $DNA_panel_config_path/DNA_panel_config.json ]; then
		DNA_panel_config_path=${DNA_panel_config_path%/*}
	else
		DNA_panel_config_path=$DNA_panel_config_path
		break
	fi
done
customized_total_cov=`$docker_run_default jq '.customized_total_cov' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
total_cov=`$docker_run_default jq '.total_cov' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
varcov_snp_illumina=`$docker_run_default jq '.varcov_snp_illumina' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
varcov_indel_illumina=`$docker_run_default jq '.varcov_indel_illumina' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
varfreq_illumina=`$docker_run_default jq '.varfreq_illumina' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
varcov_snp_Ion_torrent=`$docker_run_default jq '.varcov_snp_Ion_torrent' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
varcov_indel_Ion_torrent=`$docker_run_default jq '.varcov_indel_Ion_torrent' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
varfreq_Ion_torrent=`$docker_run_default jq '.varfreq_Ion_torrent' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
varcov_germline=`$docker_run_default jq '.varcov_germline' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
varcov_itd=`$docker_run_default jq '.varcov_itd' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
similar_precent_itd=`$docker_run_default jq '.similar_precent_itd' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
min_itd_length=`$docker_run_default jq '.min_itd_length' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
max_itd_length=`$docker_run_default jq '.max_itd_length' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
vd_germline_cf=`$docker_run_default jq '.vd_germline_cf' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
max_insertion_length=`$docker_run_default jq '.max_insertion_length' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
run_itd_analysis=`$docker_run_default jq '.run_itd_analysis' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
filter_by_bed=`$docker_run_default jq '.filter_by_bed' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
filter_by_snp=`$docker_run_default jq '.filter_by_snp' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
filter_by_MAF=`$docker_run_default jq '.filter_by_MAF' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
rm_duplicates=`$docker_run_default jq '.rm_duplicates' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
only_exon_splice=`$docker_run_default jq '.only_exon_splice' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
TVC_filter_cancer=`$docker_run_default jq '.TVC_filter_cancer' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
TVC_filter_germline=`$docker_run_default jq '.TVC_filter_germline' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
varcov_snp_vcf=`$docker_run_default jq '.varcov_snp_vcf' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
varcov_indel_vcf=`$docker_run_default jq '.varcov_indel_vcf' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
varfreq_vcf=`$docker_run_default jq '.varfreq_vcf' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
total_cov_vcf=`$docker_run_default jq '.total_cov_vcf' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
Illumina_software=`$docker_run_default jq '.Illumina_software' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
Ion_torrent_software=`$docker_run_default jq '.Ion_torrent_software' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
illumina_mapping_software=`$docker_run_default jq '.illumina_mapping_software' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
RNA_fusion_supporting_reads=`$docker_run_default jq '.RNA_fusion_supporting_reads' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
amplicon_lowest_depth_for_qc=`$docker_run_default jq '.amplicon_lowest_depth_for_qc' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
FFPE_target_seq_depth=`$docker_run_default jq '.FFPE_target_seq_depth' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
Whole_Blood_target_seq_depth=`$docker_run_default jq '.Whole_Blood_target_seq_depth' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`
Uniformity_arg_for_qc=`$docker_run_default jq '.Uniformity_arg_for_qc' $DNA_panel_config_path/DNA_panel_config.json | sed 's/\"//g'`





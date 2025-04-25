#!/bin/bash
vcf_file=$1
hpo_file_patient=$2
analysis_path=$3
dataAgent_output_file=$4
vcf_docker_sif_image=$5
mt_docker_sif_image=$6
ranking_docker_sif_image=$7
MD2GPS_docker_sif_image=$8
database_path=$9
NUM_THREADS=${10}
OPENAI_API_KEY=${11}
MODEL=${12}


# get path of this script
SCRIPT_DIR=$(cd $(dirname ${BASH_SOURCE:-$0});pwd)
DataAgent_Path=$SCRIPT_DIR
hgmd_db_file=$database_path/HGMD_Pro_2024.1_hg19.vcf
clinvar_db_file=$database_path/clinvar_20240416_hg19_add_time.vcf

# VCF annotation. input:
bash $DataAgent_Path/vcf_annotation_script/workflow_vcf_annotation.sh $vcf_file $analysis_path/results_all.txt $vcf_docker_sif_image $database_path $analysis_path

# filter variants
## 1 filter by MAF and SnpEff
singularity exec --cleanenv $vcf_docker_sif_image perl $DataAgent_Path/variant_filter/filter_by_MAF_and_SnpEff.pl $analysis_path/results_all.txt $analysis_path/results_all_filtered_maf_snpeff.txt

## 2 prediction of variant pathogenicity using MutationTaster2
singularity exec --cleanenv $mt_docker_sif_image python $DataAgent_Path/variant_filter/mt_annotation.py \
	--input_file $analysis_path/results_all_filtered_maf_snpeff.txt \
	--output_file $analysis_path/results_all_filtered_pop_mt_annotation.txt \
	--refseq_file $database_path/refseq/hg19.fa

## 3 filter by clinvar, HGMD, SIFT, and Polyphen2
singularity exec --cleanenv $vcf_docker_sif_image perl $DataAgent_Path/variant_filter/filter_variants_by_database_and_sif_polyphen2.pl $hgmd_db_file $clinvar_db_file $analysis_path/results_all_filtered_pop_mt_annotation.txt $analysis_path/results_all_filtered_pop_mt_annotation_database_sif_poly2.txt $analysis_path/results_all_notpass.txt

## 4 filter variants by DataAgent
singularity exec --cleanenv $MD2GPS_docker_sif_image python $DataAgent_Path/variant_filter/DataAgent_LLM_filter.py \
	--input_file $analysis_path/results_all_filtered_pop_mt_annotation_database_sif_poly2.txt \
	--output_file $analysis_path/DataAgent_candidate_genetic_variants.txt \
	--hpo_file_patient $hpo_file_patient \
	--hpo_file_obo $database_path/hp.obo \
	--OPENAI_API_KEY $OPENAI_API_KEY \
	--model $MODEL

## 5 Ranking of variant pathogenicity
singularity exec --cleanenv $ranking_docker_sif_image python $DataAgent_Path/variant_filter/enrichment_by_hpo_tree_multi_threads_v2.py \
	--hpo_file_patient $hpo_file_patient \
	--disease_hpo $database_path/Disease_HPO_tree.txt --disease_gene $database_path/Disease_Gene_2col.txt \
	--gene_relation $database_path/KEGG_Gene_Relation_pvalue_symbolID.txt --gene_hpo $database_path/Gene_HPO_symbol_2col.txt \
	--hpo_file_obo $database_path/hp.obo --inheritance $database_path/Disease_inheritance_v2.txt \
	--threads $NUM_THREADS --sid MD2GPS \
	--freq $database_path/Disease_HPO_frequency.txt \
	--input $analysis_path/DataAgent_candidate_genetic_variants.txt \
	--out_all $analysis_path/ranking_results_all.txt --out_filter $analysis_path/ranking_results_temp.txt --out_filter2 $analysis_path/DataAgent_gene_results.txt


## 6 explain results with natural language
singularity exec --cleanenv $MD2GPS_docker_sif_image python $DataAgent_Path/DataAgent.py \
	--input_file $analysis_path/DataAgent_gene_results.txt \
	--output_file $dataAgent_output_file \
	--hpo_file_patient $hpo_file_patient \
	--hpo_file_obo $database_path/hp.obo \
	--prompt_file $DataAgent_Path/DataAgent_prompt.json \
	--OPENAI_API_KEY $OPENAI_API_KEY \
	--model $MODEL



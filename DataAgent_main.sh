#!/bin/bash

script_root_path=/export/home/pmap/xinyang/genetic_disease/code/2_PDS
pds_image=/export/home/pmap/xinyang/genetic_disease/database/env/pathway_image.sif
out_path=/export/home/zhouxinyang/Genetic_Diagnosis/Step2/results_for_xinyang_debug

freq_file=/export/home/pmap/xinyang/genetic_disease/database/Disease_HPO_frequency.txt
disease_hpo_file=/export/home/pmap/xinyang/genetic_disease/database/Disease_HPO_tree.txt
disease_gene_file=/export/home/pmap/xinyang/genetic_disease/database/Disease_Gene_2col.txt
gene_relation_file=/export/home/pmap/xinyang/genetic_disease/database/KEGG_Gene_Relation_pvalue_symbolID.txt
gene_hpo_file=/export/home/pmap/xinyang/genetic_disease/database/Gene_HPO_symbol_2col.txt
HPO_obo_file/export/home/pmap/xinyang/genetic_disease/database/hp.obo
inheritance_file=/export/home/pmap/xinyang/genetic_disease/database/Disease_inheritance.txt


raw_data_path=/export/home/renyongyong/data/docker_share/Data/src_file/121_disease
sample_path=/export/home/renyongyong/project/PHD/rare_disease/data/121_disease/vcf_data
hgmd_db_file=/export/home/renyongyong/project/PHD/rare_disease/database/HGMD/HGMD_Pro_2024.1_hg19.vcf
clinvar_db_file=/export/home/renyongyong/data/docker_share/Database/clinvar/clinvar_20240416_hg19_add_time.vcf
sif_image=/export/home/renyongyong/software/singularity_image/ubuntu2004_MEI.sif
ref_seq=/export/home/renyongyong/project/PHD/rare_disease/database/refseq/hg19.fa

result_path=/export/home/renyongyong/project/PHD/rare_disease/data/121_disease/results/annotation
file_names=($(ls -p "$raw_data_path" | grep -v /))


samples=($(ls -F $sample_path | grep '/$'|sed 's/\///g'))

######################################################################
#####################################################################

#!/bin/bash
vcf_file=$1
hpo_file=$2
analysis_path=$3
dockers_path=$4
database_path=$5
out_file_name=$6
NUM_THREADS=$7
OPENAI_API_KEY=${8}
MODEL=gpt-4

mt_docker_sif_image=$dockers_path/ubuntu2004_MT.sif
ranking_docker_sif_image=$dockers_path/ubuntu2004_Rank.sif
vcf_docker_sif_image=$dockers_path/ubuntu1604_py3_VCF.sif
MD2GPS_docker_sif_image=$dockers_path/ubuntu2004_MD2GPS.sif

# get path of this script
SCRIPT_DIR=$(echo -n "$SCRIPT_DIR" | tr -d '\n' | sed 's:/*$::')
DataAgent_Path=$SCRIPT_DIR/DataAgent

# VCF annotation. input:
singularity exec --cleanenv $vcf_docker_sif_image bash $DataAgent_Path/workflow_vcf_annotation.sh $vcf_file $analysis_path/results_all.txt $vcf_docker_sif_image $database_path

# filter variants
## 1 filter by MAF and SnpEff
singularity exec --cleanenv $vcf_docker_sif_image perl $DataAgent_Path/variant_filter/filter_by_MAF_and_SnpEff.pl $analysis_path/results_all.txt $analysis_path/results_all_filtered_maf_snpeff.txt

## 2 prediction of variant pathogenicity using MutationTaster2
singularity exec --cleanenv $mt_docker_sif_image python $DataAgent_Path/variant_filter/mt_annotation.py \
	--input_file $analysis_path/results_all_filtered_maf_snpeff.txt \
	--output_file $folder_name_with_path/results_all_filtered_pop_mt_annotation.txt \
	--refseq_file $ref_seq

## 3 filter by clinvar, HGMD, SIFT, and Polyphen2
singularity exec --cleanenv $vcf_docker_sif_image perl $DataAgent_Path/variant_filter/filter_variants_by_database_and_sif_polyphen2.pl $hgmd_db_file $clinvar_db_file $analysis_path/results_all_filtered_pop_mt_annotation.txt $analysis_path/results_all_filtered_pop_mt_annotation_database_sif_poly2.txt $analysis_path/results_all_notpass.txt

## 4 filter variants by knowledge in GPT(等新阳修改)
singularity exec --cleanenv $MD2GPS_docker_sif_image python $DataAgent_Path/variant_filter/BioAgent_1_LLM.py \
	--input_path $analysis_path \
	--log_path $analysis_path \
	--vcf_file_name results_all_filtered.txt \
	--hpo_file_name $hpo_file \
	--hpo_dataset_path $database_path/hp.obo \
	--OPENAI_API_KEY $OPENAI_API_KEY \
	--model $MODEL

## 5 Ranking of variant pathogenicity
singularity exec --cleanenv $ranking_docker_sif_image python $DataAgent_Path/variant_filter/enrichment_by_hpo_tree_multi_threads.py \
	--hpo $hpo_file \
	--disease_hpo $database_path/Disease_HPO_tree.txt --disease_gene $database_path/Disease_Gene_2col.txt \
	--gene_relation $database_path/KEGG_Gene_Relation_pvalue_symbolID.txt --gene_hpo $database_path/Gene_HPO_symbol_2col.txt \
	--hpo_obo $database_path/hp.obo --inheritance $database_path/Disease_inheritance.txt \
	--threads $NUM_THREADS --sid MD2GPS \
	--freq $database_path/Disease_HPO_frequency.txt \
	--input $analysis_path/results_all_filtered.txt 
	--out_all $analysis_path/ranking_results_all.txt --out_filter $analysis_path/ranking_results_temp.txt --out_filter2 $analysis_path/gene_ranking_results.txt


## 6 explain results with natural language
singularity exec --cleanenv $MD2GPS_docker_sif_image python $DataAgent_Path/DataAgent.py \
    --input $analysis_path/gene_ranking_results.txt \
    --hpo_file $hpo_file \
    --prompt_file $DataAgent_Path/DataAgent_prompt.json \
    --hpo_dataset_path $database_path/hp.obo \
    --OPENAI_API_KEY $OPENAI_API_KEY \
    --model $MODEL \
    --out DataAgent_diagnosis_result.txt


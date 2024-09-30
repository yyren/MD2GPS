#!/bin/bash
in_vcf=$1
in_hpo=$2
analysis_work_path=$3
out_file_with_path=$4
NUM_THREADS=$5
OPENAI_API_KEY=$6


# get path of this script
SCRIPT_DIR=$(echo -n "$SCRIPT_DIR" | tr -d '\n' | sed 's:/*$::')
database_path=$SCRIPT_DIR/Database
dockers_path=$SCRIPT_DIR/Docker_image
mt_docker_sif_image=$docker_path/ubuntu2004_MT.sif
ranking_docker_sif_image=$docker_path/ubuntu2004_Rank.sif
vcf_docker_sif_image=$docker_path/ubuntu1604_py3_VCF.sif
MD2GPS_docker_sif_image=$docker_path/ubuntu2004_MD2GPS.sif

###1. DataAgent
bash DataAgent_main.sh $in_vcf $in_hpo $analysis_work_path $dockers_path $database_path $out_file_with_path $NUM_THREADS $OPENAI_API_KEY

### 2. KnowledgeAgent
singularity exec --cleanenv $MD2GPS_docker_sif_image python KnowledgeAgent.py \
    --PDS_result_file $PDS_RESULT_FILE \
    --hpo_file $HPO_FILE \
    --prompt_file $KNOWLEDGE_AGENT_PROMPT_FILE \
    --hpo_dataset_path $HPO_DATASET_PATH \
    --OPENAI_API_KEY $OPENAI_API_KEY \
    --model $MODEL \
    --out KnowledgeAgent_response.txt

### 3. DebateAgent
singularity exec --cleanenv $MD2GPS_docker_sif_image python DebateAgent.py \
    --PDS_result_file $PDS_RESULT_FILE \
    --hpo_file $HPO_FILE \
    --prompt_file $DEBATE_AGENT_PROMPT_FILE \
    --hpo_dataset_path $HPO_DATASET_PATH \
    --OPENAI_API_KEY $OPENAI_API_KEY \
    --model $MODEL \
    --out diagnosis_results.txt

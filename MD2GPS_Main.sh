#!/bin/bash
Project_config_with_path=$1
output_file=$2
Docker_Image_Default=$3


### input: fileName and HPO_fileName in $Project_config_with_path
### output: $output_file, e.g. Debate_diagnosis_result.txt

echo "Runing: MD2GPS_Main.sh"
default_docker_image_run="singularity exec $Docker_Image_Default"

#-----------------------------1. get the arraries of the public parameters from the $Project_config_with_path -----------------------------------------
threads=`$default_docker_image_run jq '.code | .[] | .coreNumber' $Project_config_with_path | sed 's/\"//g'`
sampleNames=`$default_docker_image_run jq '.data | .[] | .sampleName' $Project_config_with_path | sed 's/\"//g'`
fileNames=`$default_docker_image_run jq '.data | .[] | .fileName' $Project_config_with_path | sed 's/\"//g'`
HPO_fileNames=`$default_docker_image_run jq '.data | .[] | .HPO_fileName' $Project_config_with_path | sed 's/\"//g'`


imageFile_dataAgent_annotation=`$default_docker_image_run jq '.software| .[] | .imageName' $Project_config_with_path | sed 's/\"//g'`
imageFile_dataAgent_MT=`$default_docker_image_run jq '.software| .[] | .imageName2' $Project_config_with_path | sed 's/\"//g'`
imageFile_dataAgent_Rank=`$default_docker_image_run jq '.software| .[] | .imageName3' $Project_config_with_path | sed 's/\"//g'`
imageFile_dataAgent_md2gps=`$default_docker_image_run jq '.software| .[] | .imageName4' $Project_config_with_path | sed 's/\"//g'`


Analysis_Path=`$default_docker_image_run jq '.Analysis_Path' $Project_config_with_path | sed 's/\"//g'`
GPT_MODEL=`$default_docker_image_run jq '.GPT_MODEL' $Project_config_with_path | sed 's/\"//g'`
OPENAI_API_KEY=`$default_docker_image_run jq '.OPENAI_API_KEY' $Project_config_with_path | sed 's/\"//g'`
MD2GPS_Database_Path=`$default_docker_image_run jq '.MD2GPS_Database_Path' $Project_config_with_path | sed 's/\"//g'`
Max_debate_rounds=`$default_docker_image_run jq '.Max_debate_rounds' $Project_config_with_path | sed 's/\"//g'`

#-----------------------------2. set the arguments used in the customized pipeline  -----------------------------------------
script_root_path=$(cd $(dirname ${BASH_SOURCE:-$0});pwd)
software_image_annotation=$imageFile_dataAgent_annotation
software_image_MT=$imageFile_dataAgent_MT
software_image_Rank=$imageFile_dataAgent_Rank
software_image_md2gps=$imageFile_dataAgent_md2gps

infile=${fileNames[0]}
hpofile=${HPO_fileNames[0]}

####################### 3. Customized pipeline ###########################################################
########################3.1 run customized scripts #######################################################
cd $Analysis_Path

###3.1.1 DataAgent
bash $script_root_path/DataAgent/DataAgent_workflow.sh \
$infile $hpofile \
$Analysis_Path $Analysis_Path/DataAgent_diagnosis_result_round_0.txt \
$software_image_annotation \
$software_image_MT \
$software_image_Rank \
$software_image_md2gps \
$MD2GPS_Database_Path \
$threads \
$OPENAI_API_KEY \
$GPT_MODEL

echo "start KG"
### 3.1.2 KnowledgeAgent
singularity exec --cleanenv $software_image_md2gps python $script_root_path/KnowledgeAgent.py \
    --input_file $Analysis_Path/DataAgent_gene_results.txt \
    --output_file $Analysis_Path/KnowledgeAgent_diagnosis_result_round_0.txt \
    --hpo_file_patient $hpofile \
    --hpo_file_obo $MD2GPS_Database_Path/hp.obo \
    --prompt_file $script_root_path/KnowledgeAgent_prompt.json \
    --OPENAI_API_KEY $OPENAI_API_KEY \
    --model $GPT_MODEL


### 3.1.3. DebateAgent
singularity exec --cleanenv $software_image_md2gps python $script_root_path/DebateAgent.py \
    --DataAgent_result $Analysis_Path/DataAgent_gene_results.txt \
    --DataAgent_ranking_result $Analysis_Path/DataAgent_diagnosis_result_round_0.txt \
    --KnowledgeAgent_result $Analysis_Path/KnowledgeAgent_diagnosis_result_round_0.txt \
    --output_file $output_file \
    --hpo_file_patient $hpofile \
    --hpo_file_obo $MD2GPS_Database_Path/hp.obo \
    --prompt_file $script_root_path/DebateAgent_prompt.json \
    --OPENAI_API_KEY $OPENAI_API_KEY \
    --model $GPT_MODEL \
    --Max_debate_rounds $Max_debate_rounds


############################ Finish the customized pipeline #########################




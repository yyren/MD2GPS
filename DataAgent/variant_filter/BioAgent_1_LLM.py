import argparse
import os
import re
import time
import pronto
import requests
import openai
import pandas as pd
from io import StringIO
from collections import defaultdict
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from Bio import SeqIO
from openai import OpenAI
from tqdm import tqdm

########## The parameters need to run the script ##########
parser = argparse.ArgumentParser(description='The script of BioAgent prediction.')
parser.add_argument('--input_path', type=str, required=True, help='Path to the input folder.')
parser.add_argument('--vcf_file_name', type=str, required=True, help='The vcf file name in the input folder(format: sep = \n.).')
parser.add_argument('--hpo_file_name', type=str, required=True, help='The file name includes the hpo id. (format: sep = \n.)')
parser.add_argument('--hpo_dataset_path', type=str, required=True, help='Path to the hpo dataset file. (file format: .obo)')
parser.add_argument('--log_path', type=str, required=True, help='Path to save processing log file (folder). (the folders have been done.)')
parser.add_argument('--OPENAI_API_KEY', type=str, required=True, help='OpenAI API key.')
parser.add_argument('--model', type=str, required=True, help='The openai model you choose to use.')
args = parser.parse_args()

input_path = args.input_path
vcf_file_name = args.vcf_file_name
hpo_file_name = args.hpo_file_name
hpo_dataset_path = args.hpo_dataset_path
log_path = args.log_path
OPENAI_API_KEY = args.OPENAI_API_KEY
chat_model = args.model
########## set the parameters here ##########
client = OpenAI(api_key=OPENAI_API_KEY)


########## def MT-predict and  LLM-predict here ##########
##### MT_predict use the html_api and multiprocess to reduce the run time #####
def get_LLM_prediction(filtered_records, prompt, client):
    ### construct the unique gene list here, return Numpy
    unique_gene_name = filtered_records['Gene_Name'].unique()

    LLM_predict_result_dict = {}
    for gene in tqdm(unique_gene_name, desc="LLM is predicting genetic and phenotypic correlations.", leave=False):
        full_filter_input = prompt + gene
        retries = 0
        success = False

        while not success and retries < 3:
            try:
                response = client.chat.completions.create(
                    model=chat_model,
                    messages=[
                        {"role": "system", "content": "Please analyze the gene's relation to the abnormal phenotypes and respond with a simple 'Yes' or 'No'."},
                        {"role": "user", "content": full_filter_input}
                    ]
                )
                if response.choices and response.choices[0].message.content.strip():
                    LLM_prediction = str(response.choices[0].message.content).replace('\n', '')
                    ### construct the LLM prediction result dict here
                    LLM_predict_result_dict[gene] = LLM_prediction
                    success = True
                    time.sleep(1)
                else:
                    print(f"\n{gene} has no response from LLM. This is {retries} attempts.")

            except openai.RateLimitError as e:
                print(f"\nAttempt {retries + 1}: Rate Limit Exceeded - {e}")
            except openai.APIError as e:
                print(f"\nAttempt {retries + 1}: API Error - {e}")
            except openai.APIConnectionError as e:
                print(f"\nAttempt {retries + 1}: API Connection Error - {e}")
            except Exception as e:
                print(f"\nAttempt {retries + 1}: Unexpected Error - {e}")
            finally:
                if not success:
                    time.sleep(2 ** retries)
                    retries += 1

        if not success and retries >= 3:
            LLM_predict_result_dict[gene] = "no prediction"
            print("\nFailed to get a response after 3 attempts.")

    filtered_records['LLM_prediction'] = filtered_records['Gene_Name'].map(LLM_predict_result_dict).fillna("no prediction")
    # updated_records = filtered_records[(filtered_records['LLM_prediction'].str.contains("Yes", na=False)) | (filtered_records['LLM_prediction'] == "no prediction")]
    LLM_filtered_records = filtered_records[raw_records['LLM_prediction'] != 'No']

    return LLM_filtered_records, filtered_records


########## create the HPO_dict, follow the structure: HPO_dict[ID] = Name ##########
ontology = pronto.Ontology(hpo_dataset_path)
HPO_dict = {term.id: term.name for term in ontology.terms()}
processed_cases = []

log_path = os.path.join(log_path, "BioAgent_1_processed_id.txt")
########## main ###########
for id_folder in tqdm(os.listdir(input_path), desc="Analyzing the patient's cases in folder."):
    id_folder_path = os.path.join(input_path, id_folder)
    if os.path.isdir(id_folder_path):
        input_variant_path = os.path.join(id_folder_path, vcf_file_name)
        HPO_path = os.path.join(id_folder_path, hpo_file_name)

        MT_and_LLM_result_path = os.path.join(id_folder_path, f"{id_folder}_BioAgent_filtered_result_GPT4.txt")
        MT_filter_and_LLM_result_path = os.path.join(id_folder_path, f"{id_folder}_BioAgent_GPT4_annotation.txt")

        ##### construct the abnormal phenotypes list here #####
        symptom_list = []
        with open(HPO_path, 'r', encoding='utf-8') as file:
            for line in file:
                try:
                    symptom_list.append(HPO_dict[line.strip()])
                except KeyError:
                    print(f"\n{line.strip()} is wrong. Please check the hpo database.")
        symptom_str = ', '.join(symptom_list)

        filter_gene_based_symptom_prompt = f"""
            Based on your existing biological and medical knowledge, please determine whether [gene] is related to at least one abnormal phenotype in the [symptom_list]. The result should be the most frequent answer out of five judgments, and you must only answer "Yes" or "No" without any explanation.
            [phenotypes_list]: {symptom_str}
            [gene]: """

        raw_records = pd.read_csv(input_variant_path, sep='\t', keep_default_na=False)

        ##### LLM is prediction now #####
        MT_and_LLM_filter_result, MT_filter_and_LLM_result = get_LLM_prediction(raw_records, filter_gene_based_symptom_prompt, client)
        MT_filter_and_LLM_result.to_csv(MT_filter_and_LLM_result_path, sep='\t', index=False)
        MT_and_LLM_filter_result.to_csv(MT_and_LLM_result_path, sep='\t', index=False)

        ### record the done sub-folder name ###
        processed_cases.append(id_folder)
        with open(log_path, 'w', encoding='utf-8') as log_file:
            for pid in processed_cases:
                log_file.write(f"{pid}\n")

print("Analysis of BioAgent-1 finish!")

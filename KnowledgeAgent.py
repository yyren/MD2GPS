import re
import os
import time
import json
import pronto
import openai
import argparse
import pandas as pd
from openai import OpenAI


########## The parameters need to run the script ##########
parser = argparse.ArgumentParser(description='The script of BioAgent prediction.')
parser.add_argument('--PDS_result_file', type=str, required=True,
                    help='The PDS result file in the input folder(format: sep = \n.).')
parser.add_argument('--hpo_file', type=str, required=True,
                    help='The file name includes the hpo id. (format: sep = \n.)')
parser.add_argument('--prompt_file', type=str, required=True,
                    help='The file includes prompt (.json format). (format: sep = \n.)')
parser.add_argument('--hpo_dataset_path', type=str, required=True,
                    help='Path to the hpo dataset file. (file format: .obo)')
parser.add_argument('--OPENAI_API_KEY', type=str, required=True, help='OpenAI API key.')
parser.add_argument('--model', type=str, required=True, help='The openai model you choose to use.')
parser.add_argument('--out', type=str, required=True, help='The output file.')
args = parser.parse_args()


PDS_results_file = args.PDS_result_file
hpo_file = args.hpo_file
prompt_file = args.prompt_file
hpo_dataset_path = args.hpo_dataset_path
OPENAI_API_KEY = args.OPENAI_API_KEY
chat_model = args.model
outFile = args.out

client = OpenAI(api_key=OPENAI_API_KEY)

ontology = pronto.Ontology(hpo_dataset_path)
HPO_dict = {term.id: term.name for term in ontology.terms()}


def extract_gene(response_str):
    json_content_pattern = re.compile(r'\{(.*)\}', re.DOTALL)
    json_match = json_content_pattern.search(response_str)

    if json_match:
        json_content = json_match.group(0).replace('\n', '')
    else:
        json_content = response_str.replace('\n', '')

    pattern = re.compile(r'"Genetic_mutation_order": \[(.*?)\]', re.DOTALL)
    match = pattern.search(json_content)
    if match:
        gene_sequence_str = match.group(1)
        gene_sequence = [gene.strip().strip('"') for gene in gene_sequence_str.split(',')]
        return ','.join(gene_sequence)
    else:
        return "Gene sequence order not found"


def read_target_hpo(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        id_list = file.read().splitlines()
    symptom_list = [ontology[id].name for id in id_list if id in ontology]
    return ','.join(symptom_list)


def read_pds_diagnosis_result(file_path):
    df = pd.read_csv(file_path, sep='\t', keep_default_na=False, nrows=20)

    if df.empty or 'Pathogenic_Gene' not in df.columns:
        return 'no gene list order', 'no gene list', 0

    gene_mutation_pre_order = df['Pathogenic_Gene'].drop_duplicates().tolist()
    gene_list_count = len(gene_mutation_pre_order)

    unique_gene_list_str = ','.join(set(df['Pathogenic_Gene'])) if gene_list_count > 0 else 'no gene list'
    gene_mutation_pre_order_str = ','.join(gene_mutation_pre_order) if gene_list_count > 0 else 'no gene list order'

    return gene_mutation_pre_order_str, unique_gene_list_str, gene_list_count


def LLM_diagnose_pathogenic_gene_order(prompt):
    LLM_prediction = None
    retries = 0
    success = False
    while retries < 3 and not success:
        try:
            response = client.chat.completions.create(
                model=chat_model,
                messages=[
                    {"role": "system", "content": "You are a scientist with expertise in both biology and medicine, facing a complex patient case."},
                    {"role": "user", "content": prompt}
                ]
            )
            if response.choices and response.choices[0].message.content.strip():
                LLM_prediction = response.choices[0].message.content
                success = True
                return LLM_prediction
            else:
                print(f"No response from LLM. This is {retries + 1} attempts.")
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
    return "No valid response received after multiple attempts."


max_attempts = 3
KnowledgeAgent_pathogenic_gene_order = None


# read the phenotypes
if os.path.exists(hpo_file):
    phenotype_list = read_target_hpo(hpo_file)
else:
    raise FileNotFoundError(f"The hpo list file is not exist: {hpo_file}")


# read the ranking result
if os.path.exists(PDS_results_file):
    gene_list_ordered, gene_list, gene_num = read_pds_diagnosis_result(PDS_results_file)
else:
    raise FileNotFoundError(f"The DataAgent ranking result file is not exist: {PDS_results_file}")


# If there is no gene in PDS result, diagnosis end.
if gene_num == 0:
    KnowledgeAgent_pathogenic_gene_order = "Not Found"
    with open(outFile, 'w', encoding='utf-8') as file:
        file.write(f"The patient has no gene list from DataAgent rank result, diagnosis end.")
else:
    # read the KnowledgeAgent prompt
    with open(prompt_file, 'r', encoding='utf-8') as json_file:
        prompt = json.load(json_file)
        KnowledgeAgent_prompt = prompt['KnowledgeAgent_prompt']

        KnowledgeAgent_prompt_fill = KnowledgeAgent_prompt.replace("{phenotype_list_placeholder}", phenotype_list)
        KnowledgeAgent_prompt_fill = KnowledgeAgent_prompt_fill.replace("{gene_list_placeholder}", gene_list)
    
    # get the KnowledgeAgent response
    for attempt in range(max_attempts):
        KnowledgeAgent_analysis = LLM_diagnose_pathogenic_gene_order(KnowledgeAgent_prompt_fill)  # KPAgent response
        KnowledgeAgent_pathogenic_gene_order = extract_gene(KnowledgeAgent_analysis)
        KnowledgeAgent_order_num = KnowledgeAgent_pathogenic_gene_order.count(',') + 1 if KnowledgeAgent_pathogenic_gene_order not in ['no gene list order', ''] else 0
        time.sleep(1)
        if KnowledgeAgent_order_num == gene_num:
            break
        else:
            tries = attempt + 1
            all_KnowledgeAgent_analysis = f"This is the {tries} answer: " + KnowledgeAgent_analysis + f"The Genetic_mutation_order of this answer does not include all the genes in the gene list：{gene_list}. Please Answer again."
    
    # store the KnowledgeAgent response
    with open(outFile, 'w', encoding='utf-8') as file:
        file.write(KnowledgeAgent_analysis)

print("Knowledge Agent analysis is finished!")

import re
import os
import time
import json
import pronto
import openai
import argparse
import pandas as pd
import warnings
from openai import OpenAI


########## The parameters need to run the script ##########
parser = argparse.ArgumentParser(description='The script of Debate.')
parser.add_argument('--DataAgent_ranking_result', type=str, required=True,
                    help='DataAgent analysis file with gene evidence and ranking result. (txt)')
parser.add_argument('--DataAgent_result', type=str, required=True,
                    help='DataAgent results without explanation. (txt)')
parser.add_argument('--KnowledgeAgent_result', type=str, required=True,
                    help='KnowledgeAgent analysis file. (txt)')
parser.add_argument('--output_file', type=str, required=True,
                    help='Output file path and name.(txt)')
parser.add_argument('--hpo_file_patient', type=str, required=True,
                    help='The file includes the hpo id. (format: sep = \n.)')
parser.add_argument('--hpo_file_obo', type=str, required=True,
                    help='Path to the hpo dataset file. (file format: .obo)')
parser.add_argument('--prompt_file', type=str, required=True,
                    help='The file includes prompt (.json format). (format: sep = \n.)')
parser.add_argument('--OPENAI_API_KEY', type=str, required=True,
                    help='OpenAI API key.')
parser.add_argument('--model', type=str, required=True,
                    help='The openai model you choose to use.')
parser.add_argument('--Max_debate_rounds', type=int, required=True,
                    help='The max rounds of debating, which the debate will end when reaches the max rounds. ')


args = parser.parse_args()
DataAgent_result = args.DataAgent_result
DataAgent_ranking_result = args.DataAgent_ranking_result
KnowledgeAgent_ranking_result = args.KnowledgeAgent_result
output_file = args.output_file
hpo_file_patient = args.hpo_file_patient
hpo_dataset_path = args.hpo_file_obo
prompt_file_path = args.prompt_file
OPENAI_API_KEY = args.OPENAI_API_KEY
chat_model = args.model
max_rounds = args.Max_debate_rounds

warnings.filterwarnings("ignore", category=UnicodeWarning)

client = OpenAI(api_key=OPENAI_API_KEY)

ontology = pronto.Ontology(hpo_dataset_path)
HPO_dict = {term.id: term.name for term in ontology.terms()}

max_attempts = 3


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


def read_DataAgent_diagnosis_result(file_path):
    df = pd.read_csv(file_path, sep='\t', keep_default_na=False, nrows=20)

    if df.empty or 'Pathogenic_Gene' not in df.columns:
        return 'no gene list order', 'no gene list', 0

    gene_mutation_pre_order = df['Pathogenic_Gene'].drop_duplicates().tolist()
    gene_list_count = len(gene_mutation_pre_order)

    unique_gene_list_str = ','.join(set(df['Pathogenic_Gene'])) if gene_list_count > 0 else 'no gene list'
    gene_mutation_pre_order_str = ','.join(gene_mutation_pre_order) if gene_list_count > 0 else 'no gene list order'

    return gene_mutation_pre_order_str, unique_gene_list_str, gene_list_count


def read_DataAgent_evidence(file_path):
    df = pd.read_csv(file_path, sep='\t',
                     usecols=['Pathogenic_Gene', 'Amino_Change', 'Effect', 'Clinvar_latest_evidence'], nrows=20)

    required_columns = ['Pathogenic_Gene', 'Amino_Change', 'Effect', 'Clinvar_latest_evidence']
    for col in required_columns:
        if col not in df.columns:
            return f"Error: Column '{col}' is missing from the data."

    result_str = df.to_string(index=False)
    return result_str


def LLM_diagnose_pathogenic_gene_order(prompt):
    LLM_prediction = None
    retries = 0
    success = False
    while retries < 3 and not success:
        try:
            response = client.chat.completions.create(
                model=chat_model,
                messages=[
                    {"role": "system",
                     "content": "You are a scientist with expertise in both biology and medicine, facing a complex patient case."},
                    {"role": "user", "content": prompt}
                ],
                seed=928
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


### read the prompt ###
with open(prompt_file_path, 'r', encoding='utf-8') as json_file:
    prompt = json.load(json_file)
    transfer_evidence_prompt = prompt['transfer_evidence_prompt']
    Debate_prompt = prompt['Debate_prompt']
    DebateAgent_summary_prompt = prompt['DebateAgent_summary_prompt']


# read the analysis before debate
DataAgent_debate = open(DataAgent_ranking_result, 'r', encoding='utf-8').read()
KnowledgeAgent_debate = open(KnowledgeAgent_ranking_result, 'r', encoding='utf-8').read()

# We name the analysis as round 0 and record them.
with open(output_file, 'a', encoding='utf-8') as file:
    file.write(f"***DataAgent_round_0***\n\n{DataAgent_debate}\n\n")
with open(output_file, 'a', encoding='utf-8') as file:
    file.write(f"***KnowledgeAgent_round_0***\n\n{KnowledgeAgent_debate}\n\n")

# read the phenotypes
if os.path.exists(hpo_file_patient):
    phenotype_list = read_target_hpo(hpo_file_patient)
else:
    raise FileNotFoundError(f"The hpo list file is not exist: {hpo_file_patient}")

# read the ranking result of DataAgent
if os.path.exists(DataAgent_result):
    gene_list_ordered, gene_list, gene_num = read_DataAgent_diagnosis_result(DataAgent_result)
    evidence_table = read_DataAgent_evidence(DataAgent_result)
    transfer_evidence_prompt_full = transfer_evidence_prompt + evidence_table
    evidence = LLM_diagnose_pathogenic_gene_order(transfer_evidence_prompt_full)
    # print(evidence)
else:
    raise FileNotFoundError(f"The DataAgent ranking result file is not exist: {DataAgent_result}")

if gene_num == 0:
    Debate_process = "DataAgent sorting algorithm calculation has no results."
    with open(output_file, 'a', encoding='utf-8') as file:
        file.write(f"The patient has no gene list from DataAgent rank result, Debate end.")

elif gene_num > 0:
    Debate_prompt_fill = Debate_prompt.replace("{phenotype_list_placeholder}", phenotype_list)
    Debate_prompt_fill = Debate_prompt_fill.replace("{gene_list_placeholder}", gene_list)
    Debate_prompt_fill = Debate_prompt_fill.replace("{evidence_placeholder}", evidence)

    DebateAgent_summary_prompt_fill = DebateAgent_summary_prompt.replace("{phenotype_list_placeholder}", phenotype_list)
    DebateAgent_summary_prompt_fill = DebateAgent_summary_prompt_fill.replace("{gene_list_placeholder}", gene_list)
    DebateAgent_summary_prompt_fill = DebateAgent_summary_prompt_fill.replace("{evidence_placeholder}", evidence)

    for rounds in range(max_rounds):
        print(f"Debate round {rounds + 1} starts.")

        # DataAgent_Debate
        DataAgent_prompt_fill = Debate_prompt_fill.replace("{Your_opinion_placeholder}", DataAgent_debate)
        DataAgent_prompt_fill = DataAgent_prompt_fill.replace("{Other_perspectives_placeholder}", KnowledgeAgent_debate)
        for attempt in range(max_attempts):
            DataAgent_debating = LLM_diagnose_pathogenic_gene_order(DataAgent_prompt_fill)  # DataAgent Debate
            DataAgent_ranking = extract_gene(DataAgent_debating)
            DataAgent_order_num = DataAgent_ranking.count(',') + 1 if DataAgent_ranking not in ['no gene list order',''] else 0
            time.sleep(1)
            if DataAgent_order_num == gene_num and DataAgent_ranking == gene_list_ordered:
                break
            else:
                tries = attempt + 1
                Debate_prompt_fill_1 = f"This is the {tries} answer: " + DataAgent_debating + f"The Genetic_mutation_order of this answer does not include all the genes in the gene list：{gene_list}. Please Answer again."

        # KnowledgeAgent_Debate
        KnowledgeAgent_prompt_fill = Debate_prompt_fill.replace("{Your_opinion_placeholder}", KnowledgeAgent_debate)
        KnowledgeAgent_prompt_fill = KnowledgeAgent_prompt_fill.replace("{Other_perspectives_placeholder}", DataAgent_debating)
        for attempt in range(max_attempts):
            KnowledgeAgent_debating = LLM_diagnose_pathogenic_gene_order(KnowledgeAgent_prompt_fill)  # KnowledgeAgent Debate
            KnowledgeAgent_ranking = extract_gene(KnowledgeAgent_debating)
            KnowledgeAgent_order_num = KnowledgeAgent_ranking.count(',') + 1 if KnowledgeAgent_ranking not in ['no gene list order', ''] else 0
            time.sleep(1)
            if KnowledgeAgent_order_num == gene_num and KnowledgeAgent_ranking == gene_list_ordered:
                break
            else:
                tries = attempt + 1
                Debate_prompt_fill_2 = f"This is the {tries} answer: " + KnowledgeAgent_debating + f"The Genetic_mutation_order of this answer does not include all the genes in the gene list：{gene_list}. Please Answer again."

        DataAgent_debate = DataAgent_debating  # replace the DataAgent opinion with last round response
        KnowledgeAgent_debate = KnowledgeAgent_debating  # replace the KnowledgeAgent opinion with last round response

        with open(output_file, 'a', encoding='utf-8') as file:
            file.write(f"***DataAgent_round_{rounds + 1}***\n\n{DataAgent_debate}\n\n")
            file.write(f"***KnowledgeAgent_round_{rounds + 1}***\n\n{KnowledgeAgent_debate}\n\n")

        if DataAgent_ranking != KnowledgeAgent_ranking:
            # if "inconsistent" or "Inconsistent" in DebateAgent_response:
            print(f"Debate_round_{rounds + 1}: Rankings differ, Next round.")
            continue
        else:
            # Debate summary
            DebateAgent_summary_prompt_fill = DebateAgent_summary_prompt_fill.replace("{Analysis_1_placeholder}", DataAgent_debate)
            DebateAgent_summary_prompt_fill = DebateAgent_summary_prompt_fill.replace("{Analysis_2_placeholder}", KnowledgeAgent_debate)
            DebateAgent_summary = LLM_diagnose_pathogenic_gene_order(DebateAgent_summary_prompt_fill)

            print(f"Debate round {rounds + 1}: Rankings are consistent, summarizing debate.")
            with open(output_file, 'a', encoding='utf-8') as file:
                file.write(f"***DebateAgent_round{rounds + 1}***\n{DebateAgent_summary}\n")
            print("Debate ends due to consistency reached!")
            break
    else:
        print(f"Debate round {rounds + 1}: reach max rounds, summarizing debate.")
        # Debate summary
        DebateAgent_summary_prompt_fill = DebateAgent_summary_prompt_fill.replace("{Analysis_1_placeholder}", DataAgent_debate)
        DebateAgent_summary_prompt_fill = DebateAgent_summary_prompt_fill.replace("{Analysis_2_placeholder}", KnowledgeAgent_debate)
        DebateAgent_summary = LLM_diagnose_pathogenic_gene_order(DebateAgent_summary_prompt_fill)

        with open(output_file, 'a', encoding='utf-8') as file:
            file.write(f"***DebateAgent_round{rounds + 1}***\n{DebateAgent_summary}\n")
            file.write("Debate completed all rounds without consistent conclusion.\n")

        print("Debate completed all rounds without consistent conclusion.")

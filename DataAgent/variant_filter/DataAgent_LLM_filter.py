import argparse
import os
import time
import pronto
import warnings
import openai
from openai import OpenAI
from tqdm import tqdm
import pandas as pd

########## The parameters need to run the script ##########
parser = argparse.ArgumentParser(description='The script of DataAgent Filter.')
parser.add_argument('--input_file', type=str, required=True,
                    help='Path to the input file. (txt)')
parser.add_argument('--output_file', type=str, required=True,
                    help='The output file name and path you want. (txt)')
parser.add_argument('--hpo_file_patient', type=str, required=True,
                    help='The file includes the hpo id. (format: sep = \n.)')
parser.add_argument('--hpo_file_obo', type=str, required=True,
                    help='Path to the hpo dataset file. (file format: .obo)')
parser.add_argument('--OPENAI_API_KEY', type=str, required=True,
                    help='OpenAI API key.')
parser.add_argument('--model', type=str, required=True,
                    help='The openai model you choose to use.')
args = parser.parse_args()


input_file = args.input_file
output_file = args.output_file
hpo_file_patient = args.hpo_file_patient
hpo_dataset_path = args.hpo_file_obo
OPENAI_API_KEY = args.OPENAI_API_KEY
chat_model = args.model


########## set the parameters here ##########
client = OpenAI(api_key=OPENAI_API_KEY)

warnings.filterwarnings("ignore", category=UnicodeWarning)


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
                        {"role": "system",
                         "content": "Please analyze the gene's relation to the abnormal phenotypes and respond with a simple 'Yes' or 'No'."},
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

    filtered_records['LLM_prediction'] = filtered_records['Gene_Name'].map(LLM_predict_result_dict).fillna(
        "no prediction")  # NOW, filtered_records is the LLM annotation records file
    LLM_filtered_records = filtered_records[raw_records['LLM_prediction'] != 'No']

    return LLM_filtered_records, filtered_records


########## create the HPO_dict, follow the structure: HPO_dict[ID] = Name ##########
ontology = pronto.Ontology(hpo_dataset_path)
HPO_dict = {term.id: term.name for term in ontology.terms()}


filter_gene_based_symptom_prompt = """
        Based on your existing biological and medical knowledge, please determine whether [gene] is related to at least one abnormal phenotype in the [symptom_list]. You must only answer "Yes" or "No" without any explanation.
        [phenotypes_list]: {symptom_str_placeholder}
        [gene]: """


########## main ###########
if os.path.isfile(input_file):
    ##### construct the abnormal phenotypes list here #####
    if os.path.isfile(hpo_file_patient):
        symptom_list = []
        with open(hpo_file_patient, 'r', encoding='utf-8') as file:
            for line in file:
                try:
                    symptom_list.append(HPO_dict[line.strip()])
                except KeyError:
                    print(f"\n{line.strip()} is wrong. Please check the hpo database.")
        symptom_str = ', '.join(symptom_list)
    else:
        print(f"Error: The file '{hpo_file_patient}' does not exist. Please check the file path and try again.")

    filter_gene_based_symptom_prompt_fill = filter_gene_based_symptom_prompt.replace("{symptom_str_placeholder}", symptom_str)

    raw_records = pd.read_csv(input_file, sep='\t', keep_default_na=False)

    ##### LLM is prediction now #####
    LLM_filter_result, LLM_annotation_result = get_LLM_prediction(raw_records, filter_gene_based_symptom_prompt, client)
    # LLM_annotation_result.to_csv(LLM_annotation_result_path, sep='\t', index=False)
    LLM_filter_result.to_csv(output_file, sep='\t', index=False)

else:
    print(f"Error: The file '{input_file}' does not exist. Please check the file path and try again.")

print("DataAgent filtering completed!")

import argparse
import os
import re
import time
import pronto
import requests
import openai
import pandas as pd

from io import StringIO
from openai import OpenAI
from tqdm import tqdm
from collections import defaultdict
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from Bio import SeqIO


########## The parameters need to run the script ##########
parser = argparse.ArgumentParser(description='The script of BioAgent prediction.')
parser.add_argument('--input_file', type=str, required=True,
                    help='The input file name with path.')
parser.add_argument('--output_file', type=str, required=True,
                    help='The output file name with path.')
parser.add_argument('--refseq_file', type=str, required=True,
                    help='hg19.fa with path')
args = parser.parse_args()

input_file_name = args.input_file
output_file_name = args.output_file
refseq_file = args.refseq_file


##### MT_predict use the html_api and multiprocess to reduce the run time #####
def get_MT_mutation_prediction(input_file, refseq_file):
    # 定义目标URL
    target_url = 'https://genecascade.org/MT2021/MT_API102.cgi'
    
    # 生成更新后的CSV文件路径
    # updated_csv_file_path = csv_file_path.replace('.csv', '_updated.csv')
    
    # 初始化
    variant_predictions = defaultdict(set)
    max_variants_per_batch = 50
    session = requests.Session()
    retries = Retry(total=10, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
    session.mount('http://', HTTPAdapter(max_retries=retries))
    session.mount('https://', HTTPAdapter(max_retries=retries))
    fasta = SeqIO.to_dict(SeqIO.parse(refseq_file, "fasta"))
    
    raw_records = pd.read_csv(input_file, sep='\t', keep_default_na=False)
    
    def unify_prediction(prediction):
        if 'disease causing' in str(prediction).lower():
            return 'disease causing'
        elif 'polymorphism' in str(prediction).lower():
            return 'polymorphism'
        elif "nan" in str(prediction).lower():
            return 'no prediction'
        else:
            return prediction
    
    def correct_variant(chrom, ref, alt, position):
        chrom_with_prefix = f"chr{chrom}" if not str(chrom).startswith('chr') else chrom
        position_int = int(position)
        try:
            preceding_base = fasta[chrom_with_prefix][position_int - 2:position_int - 1].seq.upper()
        except KeyError:
            try:
                chrom_without_prefix = chrom.replace('chr', '')
                preceding_base = fasta[chrom_without_prefix][position_int - 2:position_int - 1].seq.upper()
            except KeyError:
                print(f"Cannot find chrom {chrom_with_prefix} or {chrom_without_prefix} in fasta.")
                return None
        new_ref = preceding_base + ref
        new_alt = preceding_base + alt
        new_position = str(position_int - 1)
        return f"{chrom}:{new_position}{new_ref}>{new_alt}"
    
    variants_transfer_dict = {}
    for start in tqdm(range(0, len(raw_records), max_variants_per_batch), desc='MT working', leave=False):
        raw_records_batch = raw_records.iloc[start:start + max_variants_per_batch]
        variants_set = {f"{row['CHROM']}:{row['ID']}{row['REF']}>{row['ALT']}" for _, row in
                        raw_records_batch.iterrows()}
        original_variants_str = ','.join(variants_set)
        
        all_errors_resolved = False
        while not all_errors_resolved:
            data = {'variants': original_variants_str, 'format': 'tsv'}
            response = session.post(target_url, data=data, timeout=None)
            
            if response.status_code == 200:
                error_lines = [line for line in response.text.split('\n') if "Variant not in VCF style:" in line]
                if not error_lines:
                    all_errors_resolved = True
                else:
                    corrected_variants_dict = {}
                    for error_line in error_lines:
                        match = re.search(r"Variant not in VCF style: '(.+)' '(.+)' '(.+)'", error_line)
                        if match:
                            ref, alt, variant_info = match.groups()
                            variant_match = re.search(r"([XYMT]|\d+):(\d+)([\dA-Z]+)>([\dA-Z]+)", variant_info)
                            if variant_match:
                                chrom, position, ref_bases, alt_bases = variant_match.groups()
                                original_variant = variant_info
                                corrected_variant = correct_variant(chrom, ref_bases, alt_bases, position)
                                corrected_variants_dict[original_variant] = corrected_variant
                                variants_transfer_dict[original_variant] = corrected_variant
                    for orig, corr in corrected_variants_dict.items():
                        original_variants_str = original_variants_str.replace(orig, corr)
            else:
                print(f"\n error! Request failed for batch starting at {start}: {response.status_code}")
                all_errors_resolved = True
        
        if all_errors_resolved:
            cleaned_lines = [line for line in response.text.split('\n') if not line.startswith('ERROR:')]
            cleaned_text = '\n'.join(cleaned_lines)
            try:
                raw_records_response = pd.read_csv(StringIO(cleaned_text), sep='\t')
                for _, row in raw_records_response.iterrows():
                    variant_key = f"{row['chr']}:{row['pos']}{row['ref']}>{row['alt']}"
                    prediction = row['prediction']
                    variant_predictions[variant_key].add(prediction)
            except KeyError as e:
                missing_column = e.args[0]
                print(f"\n error! KeyError encountered in batch starting at {start}: Column {missing_column} is missing.")
    
    for index, row in raw_records.iterrows():
        chrom = str(row['CHROM']).replace('X', '23')
        ori_variant_key = f"{chrom}:{row['ID']}{row['REF']}>{row['ALT']}"
        MT_predictions_set = variant_predictions.get(ori_variant_key, set())

        if not MT_predictions_set and ori_variant_key in variants_transfer_dict:
            MT_predictions_set = variant_predictions.get(variants_transfer_dict[ori_variant_key], {"Not Found"})

        if not MT_predictions_set:
            MT_predictions_set = {"Not Found"}
        
        MT_unified_predictions_set = {unify_prediction(pred) for pred in MT_predictions_set}
        prediction = "|".join([str(item) for item in MT_unified_predictions_set])
        raw_records.at[index, 'MutationTaster_predict'] = prediction
    
    raw_records['MutationTaster_predict'] = raw_records['MutationTaster_predict'].replace("", "no prediction")
    return raw_records
    #MT_filtered_records = raw_records[raw_records['MutationTaster_predict'] != 'polymorphism']
    
    # MT_filtered_records.to_csv(updated_csv_file_path, sep='\t', index=False)
    
    #return MT_filtered_records

MT_result = get_MT_mutation_prediction(input_file_name, refseq_file)
MT_result.to_csv(output_file_name, sep='\t', index=False)
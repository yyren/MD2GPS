from __future__ import unicode_literals
import json
import sys
import argparse
import pandas as pd
import os
import chardet
import codecs
import re

parser = argparse.ArgumentParser(
  description='generate the json file for the results' )

parser.add_argument('--default_json', type=str, default="default_json.json", 
            help="the prebuild json file with utf-8-sig format(utf-8无bom格式), default:default_json.json")

parser.add_argument('--folder_name', type=str, default="188_201908087234", 
            help="the name of the analysis folder , default:188_201908087234")
parser.add_argument('--result_json', type=str, default="file_summary.json", 
            help="the description of the results file , default:file_summary.json")
args = parser.parse_args()

default_json=args.default_json
folder_name=args.folder_name
result_json=args.result_json

#change to utf-8 format
content = open(default_json,"rb")# str类型
raw_encoding = chardet.detect(content.read())['encoding']
target_encoding='utf'
if not re.search(target_encoding,raw_encoding,re.IGNORECASE):
    print(default_json+" is "+raw_encoding+" format not in utf-8 format")
    os._exit(0)


#load_dict={}
json_file=open(default_json,'r')
load_dict = json.load(json_file,encoding=raw_encoding)
load_dict['folderName'] = folder_name

with open(result_json, 'w') as summary_f:
    json_data=json.dumps(load_dict, indent=4, ensure_ascii=False)
    summary_f.write(json_data)
summary_f.close()
json_file.close()

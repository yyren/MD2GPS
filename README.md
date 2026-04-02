# 🧬 MD2GPS

<p align="center">
  <b>MD2GPS: An LLM-Driven Multi-Agent Debate System for Mendelian Disease Diagnosis</b>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Status-Active-success">
  <img src="https://img.shields.io/badge/License-Academic-blue">
  <img src="https://img.shields.io/badge/Python-3.x-yellow">
  <img src="https://img.shields.io/badge/Docker-Supported-informational">
</p>

---

## 📌 Table of Contents

- [🧬 MD2GPS](#-md2gps)
  - [📌 Table of Contents](#-table-of-contents)
  - [🧠 Introduction](#-introduction)
  - [🚀 Key Features](#-key-features)
  - [🧪 Workflow](#-workflow)
  - [📊 Results](#-results)
  - [🖥️ Environment Requirements](#️-environment-requirements)
  - [🛠️ Usage](#️-usage)
    - [OpenAI API Key Requirements](#openai-api-key-requirements)
    - [💻 Local Version](#-local-version)
      - [🐳 Docker Images](#-docker-images)
  - [📂 Data \& Requirements](#-data--requirements)
  - [📥 Database](#-database)
    - [📌 Required Files (User-Provided)](#-required-files-user-provided)
      - [📥 ANNOVAR Database](#-annovar-database)
  - [Configuration File (Project\_Config.json)](#configuration-file-project_configjson)
    - [📌 Configuration Overview](#-configuration-overview)
    - [1️⃣ Computational Resources](#1️⃣-computational-resources)
    - [2️⃣ Input Data](#2️⃣-input-data)
    - [3️⃣ Software \& Docker Images](#3️⃣-software--docker-images)
    - [4️⃣ Paths \& Output](#4️⃣-paths--output)
    - [5️⃣ LLM Configuration](#5️⃣-llm-configuration)
  - [🚀 Command Line Usage](#-command-line-usage)
  - [🌐 Web Version](#-web-version)
      - [Required Files](#required-files)
      - [⏱️ Runtime](#️-runtime)

---

## 🧠 Introduction

Accurate diagnosis of Mendelian diseases is crucial for precision therapy and preimplantation genetic diagnosis. However, existing methods often fail to meet clinical standards or rely heavily on large-scale pretrained datasets.

We propose **MD2GPS**, an **LLM-driven multi-agent debate system** that integrates:

- 🧬 Data-driven reasoning
- 📚 Knowledge-driven inference
- 🤖 Natural language explanation via LLMs

The system enables **collaborative reasoning through agent debate**, improving diagnostic accuracy and interpretability.

---

## 🚀 Key Features

- 🔥 Multi-agent debate framework
- 🧠 LLM-based explanation generation
- ⚡ Plug-and-play modular design
- 🧬 Designed for Mendelian disease diagnosis
- 🔄 Easily extendable to other complex diseases

---

## 🧪 Workflow

![AUC](Pictures/worflow_of_MD2GPS.jpg)

---

## 📊 Results

- ✅ Evaluated on **1,185 samples (4 datasets)**
- 📈 TOP1 accuracy improved:
  - **42.9% → 66%**
- 🧬 In 72 challenging cases:
  - Identified pathogenic genes in **12 patients**
  - ⏱️ Reduced diagnostic time by **90%**

![](Pictures/performance_on_different_incidence_rates.jpg)

---
## 🖥️ Environment Requirements

The local version requires:

- Linux operating system
- Singularity / Apptainer (for running `.sif` containers)
- Python 3.x

⚠️ Note:
- The pipeline relies on Singularity containers rather than Docker directly
- Users must ensure Singularity is correctly installed before running

---

## 🛠️ Usage


### OpenAI API Key Requirements

- How to obtain: Sign up at platform.openai.com

- Environment variable: OPENAI_API_KEY


### 💻 Local Version

#### 🐳 Docker Images

Download and place in `Docker_image/`:

Windows download via web:
- `ubuntu1604_py3_VCF.sif`: https://bmap.sjtu.edu.cn/softstorage/details/54
- `ubuntu2004_Rank.sif`: https://bmap.sjtu.edu.cn/softstorage/details/56
- `ubuntu2004_MT.sif`: https://bmap.sjtu.edu.cn/softstorage/details/57
- `ubuntu2004_MD2GPS.sif`: https://bmap.sjtu.edu.cn/softstorage/details/55

linux (can be download directly by wget):
- `ubuntu1604_py3_VCF.sif`:https://bmap.sjtu.edu.cn/softstorage/download?id=54
- `ubuntu2004_Rank.sif`: https://bmap.sjtu.edu.cn/softstorage/download?id=56
- `ubuntu2004_MT.sif`: https://bmap.sjtu.edu.cn/softstorage/download?id=57
- `ubuntu2004_MD2GPS.sif`: https://bmap.sjtu.edu.cn/softstorage/download?id=55

After downloading, rename each file to the corresponding name above (e.g. `ubuntu1604_py3_VCF.sif`).

---

## 📂 Data & Requirements

## 📥 Database

We provide a set of required databases, which can be downloaded from:  
🔗 https://bmap.sjtu.edu.cn/datastorage/main/63

Linux (wget):
🔗 https://bmap.sjtu.edu.cn/datastorage/download?id=291
After downloading, rename the file "Database.zip" and unzip to "MD2GPS/Database".

⚠️ However, this repository **does not include** certain databases such as **HGMD** and **ANNOVAR**, due to their licensing restrictions. Users must obtain these resources separately.

To assist users, we provide an example file:

- `HGMD_database_example.vcf`  
  → An example illustrating the required format of the HGMD database  

---

### 📌 Required Files (User-Provided)

The following files are **not included** in this repository and must be obtained by users:

| File | Source | Location |
|------|--------|---------|
| HGMD_Pro_2024.1_hg19.vcf | HGMD website (license required) | `MD2GPS/Database/` |
| hg19.fa | UCSC Genome Browser | `MD2GPS/Database/refseq/` |

---

#### 📥 ANNOVAR Database

To download ANNOVAR-related databases:

1. Register at:  
   https://www.openbioinformatics.org/annovar/annovar_download_form.php  

2. After registration, you will receive an email with download instructions  

3. Download and extract `annovar.latest.tar.gz`

4. Place `download_annovar_databases.sh` into the extracted ANNOVAR directory  

5. (Recommended) Set the output directory in the script:

   **OUT_DIR** = MD2GPS/Database/annovar-database

6. Run the script to download the required databases:
  ```
  bash download_annovar_databases.sh
  ```

---

## Configuration File (Project_Config.json)

The `Project_Config.json` file defines all parameters required for running MD2GPS.  
Users must configure this file **before execution**.

### 📌 Configuration Overview

The configuration includes the following components:

---

### 1️⃣ Computational Resources

- **coreNumber**  
  Number of CPU cores used for analysis.

---

### 2️⃣ Input Data

- **sampleName**  
  Identifier of the sample of vcf.

- **fileName**  
  Path to the input VCF file.

- **HPO_fileName**  
  Path to the file containing HPO IDs (one id per line).

---

### 3️⃣ Software & Docker Images

- **imageName / imageName2 / imageName3 / imageName4**  
  Paths to required Docker images for different modules of MD2GPS.

---

### 4️⃣ Paths & Output

- **Analysis_Path**  
  Directory for storing analysis results.

- **MD2GPS_Database_Path**  
  Path to the MD2GPS database.

---

### 5️⃣ LLM Configuration

- **GPT_MODEL**  
  Model used for generating explanations (e.g., `gpt-4`).

- **OPENAI_API_KEY**  
  API key for accessing the LLM service.

- **Max_debate_rounds**  
  Number of debate rounds between agents.

---

## 🚀 Command Line Usage

```bash
bash ${absolute_path}/MD2GPS_Main.sh \
    Project_Config.json \
    MD2GPS_diagnosis_result.txt \
    /Docker_image/ubuntu1604_py3_VCF.sif
```

---
MD2GPS_diagnosis_result.txt: the analysis result file contains the diagnosis result and explanations.

## 🌐 Web Version

We provide a web application for **academic, personal, and non-commercial use**:  
🔗 https://bmap.sjtu.edu.cn/customanalysis/analysisdatas/107/1

#### Required Files

1. **VCF file**  
   - WGS sequencing data  
   - Example: `sample.vcf`

2. **HPO file (TXT)**  
   - One HPO ID per line  
   - Example: `sample_hpo_id.txt`

3. **Compressed package**
   - `.zip` or `.rar`
   - ❗ No subfolders

Test data can be downloaded from:  
https://bmap.sjtu.edu.cn/datastorage/details/292  

⚠️ The file names of the VCF and HPO files inside the compressed package must match those in the test data.

#### ⏱️ Runtime

- ~30 minutes per sample  
- Using 20 CPU cores (server-side)
[![Build Status](https://travis-ci.com/Cargill-Capstone2021/ADME-modeling-for-Cargill.svg?branch=main)](https://travis-ci.com/github/Cargill-Capstone2021/ADME-modeling-for-Cargill)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/licenses/MIT)

# ADME Modeling and Prediction for Non-Human Pharmaceutical Lead Screening 

**ADME Modeling and Prediction for Non-Human Pharmaceutical Lead Screening** is a capstone project that was sponsored by Cargill, Inc. and was developed by graduate students from the DIRECT data science program in the University of Washington.  

## Description
Do you know that developing a new drug takes at least 10 years, with the lead generation taking 3 to 4 years and clinical trails taking six to seven years on average? Any method that can shorten this prosses will be life saving and money saving. The goal of this project is to build a machine learning model that takes the SMILES strings of the drug candidates that are targetted for farm animals and predicts their ADME and toxicity properties in the early drug discovery process. ADME is an abbreviation in pharmacokinetics and pharmacology for "absorption, distribution, metabolism, and excretion", and describes the disposition of a pharmaceutical compound within an organism. The four criteria all influence the drug levels and kinetics of drug exposure to the tissues and hence influence the performance and pharmacological activity of the compound as a drug. With this *in silico* approach, the drug discovery process can be accelerated by 2-3 years. In addition, unnecessary cost and some of the animal experiments can be avoided.


## Get started
### Clone the repository
```
https://github.com/Cargill-Capstone2021/ADME-modeling-for-Cargill.git
```
### Install rdkit
```
! wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh
! chmod +x Miniconda3-py37_4.8.2-Linux-x86_64.sh
! bash ./Miniconda3-py37_4.8.2-Linux-x86_64.sh -b -f -p /usr/local
! conda install -c rdkit rdkit -y
import sys
sys.path.append('/usr/local/lib/python3.7/site-packages/')
```
### Install the environment

```
conda env create -f environment.yml
```
### Load the environment
```
conda activate adme
```

### Using the code
1. Data extraction
- use data_extract.py to extract data from ChemBL database

2. Data prep
- use data_prep.py to process the data obtained from the last step

3. Calculate the molecular descriptors using Descriptors.py

4. Train the default model for each ADME feature (under adme/modules)

5. Perform hyperparameter tuning


### Examples

See [example_notebook.ipynb](https://github.com/Cargill-Capstone2021/ADME-modeling-for-Cargill.git) for more details on how to use the code.

## Authors
Sara Aalinezhad - Chemical Engineering, MS

Salek Segid - Chemical Engineering, MS

Liwen Xing - Molecular Engineering & Sciences Institute, PhD


## License

This project is licensed under the MIT License - see `LICENSE` for details.




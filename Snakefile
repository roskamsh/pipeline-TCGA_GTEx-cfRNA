import datetime
import sys
import os
import pandas as pd
import json
import numpy as np

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"omic_config.yaml"

CONTRASTS = ['Liver_Tumor_TCGA-vs-Liver_Adj_Normal_TCGA','Liver_Tumor_TCGA-vs-Liver_Normal_GTEx','Liver_Adj_Normal_TCGA-vs-Liver_Normal_GTEx','Lung_Adenocarcinoma_Tumor_TCGA-vs-Lung_Adenocarcinoma_Adj_Normal_TCGA','Lung_Adenocarcinoma_Tumor_TCGA-vs-Lung_Normal_GTEx','Lung_Adenocarcinoma_Adj_Normal_TCGA-vs-Lung_Normal_GTEx','Lung_Squamous_Tumor_TCGA-vs-Lung_Squamous_Adj_Normal_TCGA','Lung_Squamous_Tumor_TCGA-vs-Lung_Normal_GTEx','Lung_Squamous_Adj_Normal_TCGA-vs-Lung_Normal_GTEx']

CONTRASTS_RUVSEQ = ['Liver_Tumor_TCGA-vs-Liver_Adj_Normal_TCGA-Liver_Normal_GTEx','Lung_Adenocarcinoma_Tumor_TCGA-vs-Lung_Adenocarcinoma_Adj_Normal_TCGA-Lung_Normal_GTEx','Lung_Squamous_Tumor_TCGA-vs-Lung_Squamous_Adj_Normal_TCGA-Lung_Normal_GTEx']

with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())

for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)


result_dirs = ['diffexp','tables']
for rule in result_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'results',rule)):
        log_out = os.path.join(os.getcwd(), 'results', rule)
        os.makedirs(log_out)
        print(log_out)

k = np.arange(1,12)

def message(mes):
    sys.stderr.write("|--- " + mes + "\n")

def get_contrast(wildcards):
    """Return each contrast provided in the configuration file"""
    return config["diffexp"]["contrasts"][wildcards.contrast]

rule all:
    input:
        expand("results/{contrast}/{contrast}-DEG.txt",contrast = CONTRASTS),
        expand("results/{contrast}_RUVseq/{contrast}-dds-RUVseqNorm-k{kvalue}.rds", contrast = CONTRASTS_RUVSEQ, kvalue = k)

rule deseq:
    input:
        "compiled_counts.txt"
    output:
        "results/{contrast}/{contrast}_dds.rds",
        "results/{contrast}/{contrast}_rlog.rds",
    params:
        metadata = config["metadata"]
    conda:
        "envs/deseq2.yaml"
    shell:
        """Rscript scripts/runDESeq2.R --counts={input} --metadata={params.metadata} --contrast={wildcards.contrast}"""

rule deseq_plots:
    input:
        dds = "results/{contrast}/{contrast}_dds.rds",
        rld = "results/{contrast}/{contrast}_rlog.rds"
    output:
        "results/{contrast}/{contrast}-PCA.pdf",
        "results/{contrast}/{contrast}-ma-plot.pdf",
        "results/{contrast}/{contrast}-heatmap-all-genes.pdf",
        "results/{contrast}/{contrast}-DEG.txt"
    conda:
        "envs/deseq2.yaml"
    shell:
        """Rscript scripts/DESeq2_plots.R --dds={input.dds} --rld={input.rld} --contrast={wildcards.contrast}"""

rule deseq_normd_init:
    input:
        "compiled_counts.txt"
    output:
        "results/{contrast}_RUVseq/{contrast}-dds.rds",
        "results/{contrast}_RUVseq/{contrast}-DEG.txt"
    params:
        metadata = config["metadata"]
    conda:
        "envs/RUVseq.yaml"
    shell:
        """Rscript scripts/RUVseq_init.R --counts={input} --metadata={params.metadata} --contrast={wildcards.contrast}"""

rule deseq_normd:
    input:
        dds = "results/{contrast}_RUVseq/{contrast}-dds.rds",
        res = "results/{contrast}_RUVseq/{contrast}-DEG.txt"
    output:
        "results/{contrast}_RUVseq/{contrast}-dds-RUVseqNorm-k{kvalue}.rds"
    conda:
        "envs/RUVseq.yaml"
    shell:
        """Rscript scripts/RUVseq_k_iteration.R --dds={input.dds} --res={input.res} --kvalue={wildcards.kvalue} --contrast={wildcards.contrast}"""


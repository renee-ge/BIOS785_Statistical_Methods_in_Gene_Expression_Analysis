# BIOS735_Methods_in_Gene_Expression_Analysis

This repository contains the data, analysis files, and final report for Renee's final project for the class "Statistical Methods in Gene Expression Analysis"

## Identifying Potential Therapeutic Targets for Nonalcoholic Steatohepatitis (NASH) with Fibrosis

With a global prevalence of 30%, Nonalcoholic Fatty Liver Disease (NAFLD) is the most common liver disease worldwide (Younossi et al., 20123). NAFLD is a condition where fat builds up in your liver. There are two forms of this disease: Nonalcoholic Fatty Liver (NAFL) and Nonalcoholic steatohepatitis (NASH). NAFL is simply fat retention in the liver, which is relatively benign compared to NASH. In NASH, however, one typically observes liver inflammation and damage along with fat in the liver. This inflammation and liver damage can lead to liver fibrosis (scarring). Liver fibrosis is the initial stage of irreversible liver damage known as cirrhosis which may also lead to cancer. (Suppli et al., 2019; U.S. Department of Health and Human Services, 2021)

For prevention and treatment of this prevalent and dangerous disease, better understanding of the disease mechanism and its progression is needed. As such, we’d like to use RNA-seq analysis techniques to examine how cell composition changes through NASH and liver fibrosis progression and identify potential gene targets for future research. In particular, at the cellular level, we are interested in hepatic stellate cells (HSCs) which are the major cell type involved in liver fibrosis. Although their current role is unclear, under healthy conditions we typically see quiescent HSCs (qHSCs). These cells can be triggered to differentiate into activated HSCs which are responsible for the production of collagen which can lead to cirrhosis. Thus, we will focus our analysis on HSCs in general. We replicated some of the results from the paper titled “Single-cell and bulk transcriptomics of the liver reveals potential targets of NASH with fibrosis” by Wang et al (2021).

## Data

| Dataset Name | Database | Accession Number|
|---|---|---|
|Single-cell RNA-seq data of human liver (17,810 cells from six healthy individuals profiled by Wang, ZY et al. and 157,619 cells from four healthy and three cirrhotic human livers of a public dataset)| ArrayExpress | E-MTAB-10553|
| Bulk RNA-seq data of liver fibrosis mouse models (Mouse hepatic stellate cells in vivo and in vitro) | GEO Database | GSE116987 |
| Human bulk RNA-seq data containing healthy normal weight (n=14) and obese (n=12) individuals, NAFL (n=15) and NASH (n=16) patients | GEO Databse | GSE126848 |

## Methods

Each of the datasets was first analyzed separately using appropriate preprocessing measures. In each dataset differential expression analysis was performed and common upregulated and downregulated genes between fibrotic and healthy cells among the four datasets were obtained. Gene ontology enrichment analysis was performed and the translatability of liver fibrosis associated genes between human and mouse models was discussed. 

## Project Organization

* Report.pdf contains the final written report detailing the full data, methods, analyses, and conclusions
* 

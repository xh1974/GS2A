GS2A
====

Gene Signature Association Analysis

INTRODUCTION

Gene Signature Association Analysis (GS2A) is a tool for screening genes that are associated with a set of "signature genes" in compendium of tumor cohorts of expression profiles. 

SYSTEM REQUIREMENT

Linux environment with gcc compiler

INSTALLATION

Under GS2A/ directory, run:

$./make

An executable file GS2A will be created under GS2A/bin/

RUNNING GS2A

1. Command line

GS2A -d <expression data file> -t <target gene id file> -c <candidate data file> -o <output file>

2. Format of Expression data file

A data matrix of expression data with header:

gene  <sample1> <sample2> ...
<gene1> <value11> <value12> ...
<gene2> <value21> <value22> ...
...

See GS2A/bin/example_expression.txt for example

3. Format of targe gene id file

List of target gene ids

<gene1>
<gene2>
...

See GS2A/bin/example_target.txt for example

4. Format of candidate data file

The same format as the expression data file. 

Note: The set of samples in the candidate data file is NOT necessary to be identical as the sample set in the expression file. However, there has to be significant overlap between two sample set. GS2A automatically determines the overlapping samples and uses them for analysis. Similarly, the set of genes in the candidate data file is NOT necessary to be identical to the gene set in expression data file. The type of genomic data in candidate data file is also NOT necessary to be the gene expression profiles. For example, user can generate copy number profiles in the candidate file, and input gene expression profiles in the expression data file. In such case, GS2A will screen for genes whose copy number abberations are associated with the expresssion of target genes.

See GS2A/bin/example_cand.txt for example.

5. Format of output file

<Gene ID> <is signature gene> <GS2A score>  <p-value>

Gene ID: the gene name shown in the candidate data file
is signature gene: If the gene is in the target gene set, flag 1; otherwise flag 0
score: the GS2A score
p-value: p-value computed from z-score.

Note: The p-value in the output file is for reference only. GS2A relies on Robust Rank Aggregation to determine statistical significance.

REFERENCE

Xu et al. "Gene signature association analysis identifies an EZH2-E2F1 collaborative pathway in cancer" Manuscript in Preparation.


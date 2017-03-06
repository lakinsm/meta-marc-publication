# Metagenomic Markov Models for Antimicrobial Resistance Characterization - Publication Repository

## Directories and Files

**analytic_data**

  * alignment\_cov.csv - this file contains the gene coverages from all SAM files used in the Alignment analyses
  * mmarc\_publication\_analytic\_data.csv - results from the Soil, Pediatric, and PRJNA292471 analyses
  * mmarc\_publication\_simulation\_results.csv - results from the Meta-MARC simulation cross-validation analysis
  * mmarc\_mismatch\_analytic\_data.csv - results from the genetic variation analysis of the PRJNA292471 dataset

**graphs** - exploratory and publication figures

**metadata**

  * mmarc\_model\_annotations.tsv - metadata describing the Meta-MARC annotation graph as it relates to the Meta-MARC models
  * pediatric\_soil\_metadata.csv - a mapping of the Pediatric and Soil sample names to the phenotypic truth label of each sample
  * resfams\_metadata.txt - metadata describing the resfams model annotations translated into MEGARes annotation terms
  * PRJNA292471\_metadata.csv - metadata describing the PRJNA292471 samples

**scripts**

  * alignment\_mismatch\_subroutine.py - Python script used to parse alignment SAM files for the genetic variation analysis of the PRJNA292471 dataset
  * fix\_formating.sh - Bash script to change short-hand sample names to publication sample names
  * mismatch\_regression.R - R script used to analyze the genetic variation data from the PRJNA292471 dataset
  * mmarc\_contig\_mismatch\_subroutine.py - Python script used to parse Meta-MARC assembled SAM files for the genetic variation analysis of the PRJNA292471 dataset
  * mmarc\_raw\_mismatch\_subroutine.py - Python script used to parse Meta-MARC HMMER tblout files for the genetic variation analysis of the PRJNA292471 dataset
  * mmarc\_reads\_aligned.py - Python script used to parse Meta-MARC assembled SAM files for the Soil, Pediatric, and PRJNA292471 dataset analyses
  * publication\_figures.R - R script used to generate publication figures
  * resfams\_contig\_mismatch\_subroutine.py - Python script used to parse Resfams SAM files for the genetic variation analysis of the PRJNA292471 dataset
  * resfams\_reads\_aligned.py - Python script used to parse Resfams SAM files for the Soil, Pediatric, and PRJNA292471 dataset analyses

All data is under an MIT license.  If you use this data, please cite the associated publication (TBD).


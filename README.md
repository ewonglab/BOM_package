# BOM: predictive models to study motif composition of cis-regulatory sequences 
BAG-of-Motifs (BOM) is an R package that uses XGBoost to construct predictive models of cis-regulatory sequences between cell states. It then interprets the importance of each motif in the model using SHAP scores.

## Installation

For installing and loading BOM, run:
```
install.packages("~/capsule/data/BagOfMotifs_0.0.2.tar.gz", repos = NULL, type = "source")
install.packages("cowplot")
install.packages("cvAUC")
```
## Quick start

## Tutorial

Please follow <a href="tutorial.md"> this vignette </a> instructions for a demonstration on how to use this softare. Alternative you can follow same instructions via this <a href="tutorial.ipynb">jupyter notebook</a> link.


Alternatively the tutorial can be run via <a href="https://codeocean.com/capsule/4079053/tree">this code ocean capsule</a> where all the software is preinstalled and ready to go. 
Clone the capsule and then run the jupyternotebook called tutorial_model_EW_20231201.ipynb

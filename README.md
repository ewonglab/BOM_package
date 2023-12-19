# BOM: predictive models to study motif composition of cis-regulatory sequences 
BAG-of-Motifs (BOM) is an R package that uses XGBoost to construct predictive models of cis-regulatory sequences between cell states. It then interprets the importance of each motif in the model using SHAP scores.

## Installation

For installing and loading BOM, run:
```
install.packages("~/capsule/data/BagOfMotifs_0.0.2.tar.gz", repos = NULL, type = "source")
install.packages("cowplot")
install.packages("cvAUC")
```
Install FIMO 

See <a href="https://meme-suite.org/meme/doc/install.html"> here </a> 

## Quick start

```
motifs_path <- "./data/gimme.vertebrate.v5.0.meme"
chr_sizes <- './data/mm10.chrom.sizes'
annot <- './data/mm10.knownGene.gtf.gz'

                        
# Generate FASTA 
generateAllFasta(bedDir = "./bed/", genome = "Mmusculus")

# Annotate  
BagOfMotifs::runFIMO(motifs_path = motifs_path, FIMO_path = '/path/to/fimo')

# Motif count and model training
BagOfMotifs::binModel()

# Prediction and performance
BagOfMotifs::predict_binary_multi()

# Estimate SHAP 
save_shap_multi()

# Plots (bar/beeswarm/waterfall)
shapPlots_multi()
shapPlots_multi(plotType = "beeswarm")
shapPlots_multi(plotType = "watterfall")

```


## Tutorial

<a href="tutorial.html"> How to use BOM to analyze single cell ATAC data </a>  


Alternatively, an example can be run via <a href="https://codeocean.com/capsule/4079053/tree">this code ocean capsule</a> where all the software is preinstalled and ready to go. 
Clone the capsule and then run the jupyternotebook called tutorial_model_EW_20231201.ipynb

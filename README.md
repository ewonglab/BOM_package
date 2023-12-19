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
See https://meme-suite.org/meme/doc/install.html

## Quick start

```


#motifs_path <- "./data/gimme.vertebrate.v5.0.meme"
#chr_sizes <- './data/mm10.chrom.sizes'
#annot <- './data/mm10.knownGene.gtf.gz'

# Generates BED file for each cell type
BagOfMotifs::filterCREs(inputBedFile = input_bed, annotFile = annot, chrSizes = chr_sizes)
                        
# Generate FASTA 
generateAllFasta(bedDir = "./bed/", genome = "Mmusculus")

# Annotate motifs 
BagOfMotifs::runFIMO(input_path = ./fasta/, motifs_path = motifs_path, 
                    FIMO_path = '/path/to/fimo')

# Motif count and model training
BagOfMotifs::binModel(target_ct = NULL,
                      data_path = motif_out,
                      qval_thresh = 0.5, 
                      outDir = "/results/", nthreads = 2)

# Prediction and plots
BagOfMotifs::predict_binary_multi(inputMotif_dir = "./results/"
                                  , inputXGB_dir = "./results/"
                                  , outputTrain_dir = "./results/"
                                  , pred_dir = "./results/"
                                  , outputFile = out_file
                                  )

# Estimate SHAP 
save_shap_multi(dataDir = "./results/", outDir = "./results/")

# Bar plots
shapPlots_multi(dataDir = "./results/")

# Beeswarm plots
shapPlots_multi(dataDir = "./results/", plotType = "beeswarm")

# Waterfall plots
shapPlots_multi(dataDir = "./results/", plotType = "watterfall"
                , = CRE_ids = )

```


## Tutorial

<a href="tutorial.html"> How to use BOM to analyze single cell ATAC data </a>  


Alternatively, an example can be run via <a href="https://codeocean.com/capsule/4079053/tree">this code ocean capsule</a> where all the software is preinstalled and ready to go. 
Clone the capsule and then run the jupyternotebook called tutorial_model_EW_20231201.ipynb

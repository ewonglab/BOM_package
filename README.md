# BOM: predictive models to find important motifs at cis-regulatory regions 
BAG-of-Motifs (BOM) is an R package that uses XGBoost to construct predictive models of cis-regulatory sequences between cell states. It then interprets the importance of each motif in the model using SHAP scores.

## Installation

For installing and loading BOM, run:
```
devtools::install_github("ewonglab/BOM_package")
library(BagOfMotifs)
```

If GenomicRanges and GenomicFeatures are not installed:
```
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicFeatures")
```

Our tutorial requires the mouse mm10 genome. If you don't have this already:
```
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
```

FIMO is required:
Installation information can be found <a href="https://meme-suite.org/meme/doc/install.html"> here </a> 

## Quick start

Input is a text file with four columns (chromosome, start, stop, condition)

```
library(BOM_package)
library("BSgenome.Mmusculus.UCSC.mm10")

motifs_path <- "./extdata/gimme.vertebrate.v5.0.meme"
chr_sizes <- './extdata/mm10.chrom.sizes'
annot <- './extdata/mm10.knownGene.gtf.gz'
FIMO_path <- '/path/to/fimo'
                    
# Generate FASTA and annotate motifs
generateAllFasta(bedDir = "./bed/", genome = "Mmusculus")
BagOfMotifs::runFIMO(motifs_path = motifs_path, FIMO_path = FIMO_path)

# Motif count and model training
BagOfMotifs::binModel()

# Prediction and performance
BagOfMotifs::predict_binary_multi()

# Estimate SHAP 
save_shap_multi()

# Plots (bar/beeswarm/waterfall)
shapPlots_multi()
shapPlots_multi(plotType = "beeswarm")
shapPlots_multi(plotType = "waterfall")

```


## Tutorial

<a href="tutorial.html"> How to use BOM to interrogate distal cis-regulatory elements from single cell ATAC data </a>  
<br>
<br>
A <a href="https://codeocean.com/capsule/4079053/tree"> CodeOcean example</a> is available where all the software is preinstalled and ready to go.
 
To use, clone the capsule and run the jupyternotebook "tutorial_model_EW_20231201.ipynb"

# BOM: predictive models to find important motifs at cis-regulatory regions 
BAG-of-Motifs (BOM) is an R package that uses XGBoost to construct predictive models of cis-regulatory sequences between cell states. It then interprets the importance of each motif in the model using SHAP scores.

It is based on the principle that the activity of these regions relies on the binding of TFs to specific motifs. Leveraging available TF binding motif profiles and the Extreme Gradient Boosting (XGBoost) algorithm, BOM has achieved high performance in classifying context-specific cis-regulatory elements. 

Through the use of SHapley Additive exPlanations (SHAP), BOM helps identify the important TF binding motifs that contribute to the classification. BOM provides several visualization options for motif counts and motif importance scores. They allow users to explore and interpret the most influential motifs learned by BOM, providing insights into the regulatory landscape of the analyzed cis-regulatory regions.


## Installation

For installing and loading BOM, run:
```
devtools::install_github("ewonglab/BOM_package")
library(BagOfMotifs)
```


The tutorial requires the mouse mm10 genome. If you don't have this already:
```
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
```

FIMO is required: 
Installation information can be found <a href="https://meme-suite.org/meme/doc/install.html"> here </a> 

## Quick start

Input is a text file with four columns (chromosome, start, stop, condition)

```
library(BOM_package)
library("BSgenome.Mmusculus.UCSC.mm10") # or genome of choice

extdata_path <- system.file("extdata",package = "BagOfMotifs")
motifs_path <- paste0(extdata_path,"/gimme.vertebrate.v5.0.meme")
chr_sizes <- paste0(extdata_path,"/mouse.chrom.sizes.txt")
annot <- paste0(extdata_path, "/Mus_musculus.GRCm38.92.gtf.gz")
FIMO_path <- "/path/to/fimo"
                    
# Generate FASTA and annotate motifs
generateAllFasta(bedDir = "path/to/input", genome = "Mmusculus")

# Run FIMO to scan motifs
# For FIMO ≥ 5.5.4 (default):
BagOfMotifs::runFIMO(
	motifs_path = motifs_path,
	FIMO_path   = FIMO_path
)

# For FIMO < 5.5.4, disable the updated parser:
BagOfMotifs::runFIMO(
	motifs_path   = motifs_path,
	FIMO_path     = FIMO_path,
	FIMO_updated  = FALSE
)

# Motif count and model training
BagOfMotifs::binModel()

# Prediction and performance
BagOfMotifs::predict_binary_multi()

# Estimate SHAP 
save_shap_multi()

# Plots (bar, beeswarm, waterfall)
shapPlots_multi()
shapPlots_multi(plotType = "beeswarm")
shapPlots_multi(plotType = "waterfall")
```

## Tutorial

<a href="tutorial.md"> How to use BOM to interrogate distal cis-regulatory elements from single cell ATAC data </a>  
<br>
<br>
A <a href="https://codeocean.com/capsule/4079053/tree"> CodeOcean example</a> is also available where all the software is preinstalled and ready to go.
 
To use, clone the capsule and run the jupyternotebook "tutorial_model_EW_20231201.ipynb"

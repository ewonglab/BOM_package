##############################################################################
## shapPlots
#' Makes SHAP plots for a single XGB model
#'
#' @param xgb_model  	    XGB model.
#' @param ts training set for the model provided to xgb_model.
#' @param plotType  Plottype, can be "bar" (default), "beeswarm", or "waterfall".
#' @param annotData Dataframe with columns "Motif" and "Factor". The values in "Motif" should match what was used 
#' @param max_display  Number of transcription factors to plot. The remaining transcriptions factors are plotted under other 
#' @param order  Ïf set to "decreasing" then plots the highest value first through to lowest value.
#' @param show_numbers  Boolean. If set to TRUE will display numbers on plot.
#' @param average_shap Boolean. Whether to average shap values across the specified CREs in waterfall plots (default TRUE).
#' 
#' @example
#' 
#' xgb_model <- readRDS("Cardiomyocytes.rds")
#' train_set <- "cardiom_trainSet.txt"
#' ts <- read.table(train_set, header = TRUE)
#' 
#' p <- shapPlots_test(xgb_model = xgb_model, ts = ts, CRE_id = "12:98725135-98725635", plotType = "waterfall", annotDat = NULL, annotLength = 30, order = "decreasing",show_numbers = FALSE)
#'  
#'
#' @export
shapPlots <- function(xgb_model, ts, plotType = "bar", max_display = 15, CRE_ids = NULL, annotDat = NULL, annotLength = 30
                           , order = "decreasing", show_numbers = FALSE, average_shap = TRUE,  ...)
{	
  require(ggplot2)
  require(shapviz)
  require(gggenes)
  require(shades)
  
  # removing labels
  ts$binary_celltype <- NULL
  shp <- shapviz::shapviz(object = xgb_model, data.matrix(ts))
  
  decreasing=TRUE
  if (order=="decreasing")
  {	decreasing = FALSE } 
  
  p <- {}
  if (plotType == "bar")
  {
    p <- shapviz::sv_importance(shp, kind="bar",fill="red", max_display = max_display, ...) + 
      theme(panel.background = element_blank(), axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black")) +
      geom_vline(xintercept = 0, linetype = "solid", color = "grey")
    if (! is.null(annotDat))
    {  p$data$feature <- annonTATE(p$data, annotDat, annotLength=annotLength, decreasing=decreasing)  }
    return(p)
  }
  else if (plotType == "beeswarm")
  {  	p <- shapviz::sv_importance(shp, kind="beeswarm", max_display = max_display, ...) +
    scale_color_gradient(low = "blue", high = "red") +
    theme(panel.background = element_blank()
          , axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"))
  
  if (! is.null(annotDat))
  { p$data$feature <- annonTATE(p$data, annotDat, annotLength=annotLength,decreasing=decreasing)   }
  return(p)
  }
  else if (plotType =="waterfall")
  {	
    if(length(CRE_ids) > 1){
      
      if(length(CRE_ids[!CRE_ids %in% rownames(ts)]) > 0 )
      { warning("Some CRE ids are not present in training set.")}
      if(length(CRE_ids[duplicated(CRE_ids)]) > 0 )
      {
        warning("Duplicated values present in CRE_ids. Unique CRE IDs will be used.")
        CRE_ids <- unique(CRE_ids)
        }
      CRE_idx <- which(rownames(ts) %in% CRE_ids)
      if(average_shap){
        p <- shapviz::sv_waterfall(shp, row_id = CRE_idx, fill_colors = c("blue", "red"), max_display = max_display, ...) + theme(panel.background = element_blank(),panel.grid.major.y = element_blank()) + 
          theme(panel.background = element_blank(),panel.grid.major.y = element_blank())
        if (! is.null(annotDat))
        {	p$data
          df <- data.frame(feature = row.names(p$data), value=0)
          p$data$label <- annonTATE(df, gimme_annot, waterfallOther = p$data["other","label"]
                                    , annotLength=annotLength,decreasing=decreasing)
        }
        return(p)
      }
      else
      {
          p_watterfall <- list()
          for(i in 1:length(CRE_idx)){
            print(CRE_idx[i])
            p <- shapviz::sv_waterfall(shp, row_id = CRE_idx[i], fill_colors = c("blue", "red"),...) + theme(panel.background = element_blank(),panel.grid.major.y = element_blank()) + 
              theme(panel.background = element_blank(),panel.grid.major.y = element_blank()) +
              ggtitle(rownames(ts)[CRE_idx[i]])
            if (! is.null(annotDat))
            {	p$data
              df <- data.frame(feature = row.names(p$data), value=0)
              p$data$label <- annonTATE(df, gimme_annot, waterfallOther = p$data["other","label"]
                                        , annotLength=annotLength,decreasing=decreasing)
            }
            p_watterfall[[i]] <- p
          }
      }
      return(p_watterfall)
    }
    else if(length(CRE_ids) == 1)
    {
      CRE_idx <- which(rownames(ts) == CRE_ids)
      print(CRE_idx)
      p <- shapviz::sv_waterfall(shp, row_id = CRE_idx, fill_colors = c("blue", "red"),...) + theme(panel.background = element_blank(),panel.grid.major.y = element_blank()) + 
        theme(panel.background = element_blank(),panel.grid.major.y = element_blank()) +
        ggtitle(CRE_id)       
      if (! is.null(annotDat))
      {	p$data
        df <- data.frame(feature = row.names(p$data), value=0)
        p$data$label <- annonTATE(df, gimme_annot, waterfallOther = p$data["other","label"]
                                  , annotLength=annotLength,decreasing=decreasing)
      }
      return(p)
      }
  }
}

##############################################################################
## shapPlots_multi
#' Makes SHAP plots for multiple models
#'
#' @param xgb_models  	List of file names for models
#' @param train_sets	List of file names for training sets (same order as models).
#' @param plotType  Plottype, can be "bar" (default), "beeswarm", or "waterfall"
#' @param annotData Dataframe with columns "Motif" and "Factor". The values in "Motif" should match what was used 
#' @param numTF  Number of transcription factors to plot. The remaining transcriptions factors are plotted under other 
#' @param order  Ïf set to "decreasing" then plots the highest value first through to lowest value.
#' @param show_numbers  Boolean. If set to TRUE will display numbers on plot.
#' @param average_shap Boolean. Whether to average shap values across the specified CREs in waterfall plots (default TRUE).
#' 
#' @example
#' 
#' train_list <- c("allantois_trainSet.txt", "cardiom_trainSet.txt", "endothelium_trainSet.txt", "erythroid_trainSet.txt" , "exe_endo_trainSet.txt")
#' model_list <- c("allantois.rds", "cardiom.rds", "endothelium.rds", "erythroid.rds", "exe_endo.rds")
#' 
#' plots <- shapPlots_multi(xgb_models = model_list, train_sets = train_list)
#' plots <- shapPlots_multi(xgb_models = model_list, train_sets = train_list, plotType = "beeswarm")
#' plots <- shapPlots_multi(xgb_models = model_list, train_sets = train_list, plotType = "waterfall", CRE_ids = c("19:10232248-10232748", "2:157923362-157923862"), average_shap = F)
#'  
#'
#' @export
shapPlots_multi <- function(xgb_models, train_sets, plotType = "bar", CRE_ids = NULL
                            , annotDat = NULL, annotLength = 30, order = "decreasing"
                            , show_numbers = FALSE, average_shap = TRUE, ...)
{
  xgb.models <- lapply(xgb_models, readRDS)
  train.sets <- lapply(train_sets, read.table, header = TRUE)
  
  if(length(xgb_models) != length(train_sets)){
    stop("The number of models and training sets is different.")
  }
  
  p_list <- list()
  for(i in 1:length(xgb_models)){
    
    p_list[[i]] <- shapPlots_test(xgb_model = xgb.models[[i]], ts = train.sets[[i]], plotType = plotType
                                  , CRE_ids = CRE_ids, annotDat = annotDat, annotLength = annotLength
                                  , order = order, show_numbers = show_numbers, average_shap = average_shap)
  }
  return(p_list)
}

##############################################################################
## annonTATE
#' Re-annotates a ggplot created from shapviz. 
#' Appends gene annotation onto exsiting annotation.
#'
#' @param plotDat              Data slot from a ggplot object
#' @param peakAnnotations      Motif annotation containing columns called "Motif" and "Factor"
#' @param annotLength       Number of characters to truncate annotation to. Default 30 characters.
#' @param waterfallOther      Annotation for waterfall plots "other" category
#' @param decreasing         Plot in decreasing order (default FALSE).
#'
annonTATE <- function(plotDat, peakAnnotations, annotLength = 30, waterfallOther=NULL, decreasing=FALSE)
{
    otherAnnotation <- 'other'
    if (! is.null(waterfallOther))
    {  otherAnnotation <- waterfallOther  }
    
    peakAnnotations$Motif <- gsub(pattern = "-", replacement = ".", x = peakAnnotations$Motif) # hack to ensure that matching occurs
    lookupIDs <- as.character(plotDat$feature)
    TF_IDs    <- {}
    for(i in 1:length(lookupIDs))
    {
        if (lookupIDs[i] == "other")
        {    
            TF_IDs <- c(TF_IDs,otherAnnotation)
        }
        else
        {     idx    <- which(peakAnnotations$Motif %in% lookupIDs[i]) 
            if(length(idx) > 0)
            {  TF_IDs <- c(TF_IDs, paste0(unlist(peakAnnotations$Factor[idx]), collapse=",")) }
            else  
            {     TF_IDs <- c(TF_IDs, lookupIDs[i]) }
        }
    }

    plotDat$feature <- substr(paste0(plotDat$feature,"_", TF_IDs), 1, annotLength)
    
    sum_by_feature <- aggregate(value ~ feature, data = plotDat, sum)
    sum_by_feature_ordered <- sum_by_feature[order(sum_by_feature$value, decreasing = decreasing),]
    sum_by_feature_ordered

    newFeatures <- factor(plotDat$feature, levels= unique(sum_by_feature_ordered$feature))
    

    return(newFeatures)

}

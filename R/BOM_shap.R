





##############################################################################
## plotShapWaterFall
#' Plots a waterfall plot of SHAP values for a selected locus
#'
#' @param scores  A data frame SHAP scores for a particular locus. Colnames must contain motif ID and
#'                    match entries inthe Motif column of peakAnnotations data frame.
#' @param peakAnnotations  Dataframe of transcription factor lookup ID information. Must have columns
#'                        names Motif and Factor. Values in Motif column should match column names in 
#'                        the scores dataframe input.
#' @param annotLength  Integer indicating length of annotation. Default is 30 characters.
#' @param numTF  Number of transcription factors to plot. The remaining transcriptions factors are plotted under other 
#' 
#' @example
#'
#' 
#'
#' @export
plotShapWaterFall <- function(scores,peakAnnotations,numTF=10, annotLength=30)
{
  
	scores <- t(scores)
	scores <- scores[order(abs(scores),decreasing=TRUE),]
	shapScoreSum <- sum(scores)

	waterFallDat <- data.frame(head(scores,numTF))
	colnames(waterFallDat ) <- "score"

	# Prepare annotation labels if requested
	if(! is.null(peakAnnotations))  # 
	{
		lookupIDs <- row.names(waterFallDat)
		TF_IDs <- {}
	  for(i in 1:length(lookupIDs))
	  {
		idx <- which(peakAnnotations$Motif %in% lookupIDs[i]) 
		TF_IDs <- c(TF_IDs, paste0(unlist(peakAnnotations$Factor[idx]), collapse=","))
	  }
	  row.names(waterFallDat) <- make.unique(TF_IDs, sep="__")
	}
	
	# Prepare waterfall plot
	waterFallDat$start <- 0
	waterFallDat$end <- 0	
    for(i in 1:nrow(waterFallDat))
    {
      waterFallDat$end[i]   <- shapScoreSum
      waterFallDat$start[i] <- waterFallDat$end[i] - waterFallDat$score[i]
      shapScoreSum <- waterFallDat$start[i]
    }  
	waterFallDat <- rbind(waterFallDat,data.frame(score=0, start=0, end= shapScoreSum))

	rownames( waterFallDat)[nrow(waterFallDat)] <- "Other" 
	rownames( waterFallDat)  <- substr(x=rownames( waterFallDat), start=1,stop=30)	
	waterFallDat$TF <- factor(rownames(waterFallDat),levels=rev(rownames(waterFallDat)))
	waterFallDat$id <- nrow(waterFallDat):1
	waterFallDat$barColour = "red"
	waterFallDat$barColour[which(waterFallDat$score < 0)] = "blue"

	p <- ggplot(waterFallDat, aes(TF,fill=barColour)) + 
		geom_rect(aes(TF,	xmin = id - 0.45, xmax = id + 0.45, ymin = end,	ymax = start),show.legend = FALSE) + 
		scale_fill_manual("legend", values = c("red" = "red", "blue" = "blue", "black" = "black")) +
		theme_classic() + coord_flip() #+ scale_x_reverse() +
    ggtitle(paste0("shap score ",colnames(scores)[1]))
	
	return(p)
}

##############################################################################
## ShapBarPlot
#' Plots a bar plot of SHAP values 
#'
#' @param scores  A data frame SHAP scores for a particular locus. Colnames must contain motif ID and
#'                    match entries inthe Motif column of peakAnnotations data frame.
#' @param peakAnnotations  Dataframe of transcription factor lookup ID information. Must have columns
#'                        names Motif and Factor. Values in Motif column should match column names in 
#'                        the scores dataframe input. If NULL no annotation is provided.
#' @param annotLength  Integer indicating length of annotation. Default is 30 characters.
#' @param numTF  Number of transcription factors to plot. The remaining transcriptions factors are plotted under other 
#' 
#' @example
#'
#' 
#'
#' @export
ShapBarPlot <- function(scores,peakAnnotations=NULL,numTF=10, annotLength=30)
{
	scores = (data.frame(colMeans(abs(scores))))
	rn <- row.names(scores)
	idx <- order(abs(scores),decreasing=TRUE)
	scores <- scores[idx,]
	names(scores) <- rn[idx]


    shapScoreSum <- sum(abs(scores))
    
    # Prepare data frame with requested number of entries
    plotDat <- data.frame(head(scores,numTF))
    colnames(plotDat ) <- "score"

	if(! is.null(peakAnnotations))  # 
	{
		# Lookup transcription factor ID and relabel rows
		lookupIDs <- row.names(plotDat)
		TF_IDs    <- {}
		for(i in 1:length(lookupIDs))
		{
			idx    <- which(peakAnnotations$Motif %in% lookupIDs[i]) 
			TF_IDs <- c(TF_IDs, paste0(unlist(peakAnnotations$Factor[idx]), collapse=","))
		}    
		row.names(plotDat) <- make.unique(TF_IDs, sep="__")
	}
    
    # Format dataframe
	rownames(plotDat) <- substr(x=rownames(plotDat), start=1,stop=30)
    plotDat$TF <- factor(rownames(plotDat),levels=rev(rownames(plotDat)))
    plotDat$id <- nrow(plotDat):1
    plotDat$barColour = "red"
    print(plotDat$score[1])
    
    p <- ggplot(plotDat, aes(score,fill=barColour)) + 
        geom_rect(aes(TF,	xmin = id - 0.45, xmax = id + 0.45, ymin = 0,	ymax = plotDat$score),show.legend = FALSE) + 
        scale_fill_manual("legend", values = c("red" = "red", "blue" = "blue", "black" = "black")) +
        theme_classic() + coord_flip() #+ scale_x_reverse() +
        ggtitle(paste0("shap score ",colnames(scores)[1]))
    return(p)
    
    
}
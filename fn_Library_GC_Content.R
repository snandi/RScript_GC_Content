########################################################################
## This function library is for analyzing intensity with gc content   ##
########################################################################

######################### Convert NAs to Zero ##########################
na.is.zero <- function(X)
{
  X1 <- X
  X1[is.na(X)] <- 0.0
  return(X1)
}
########################################################################

########################################################################
"%notin%" <- function(x, y){
  if(x %in% y){
    return(FALSE)
  } else{
    return(TRUE)
  }
}
########################################################################

########################################################################
"%w/o%" <- function(x, y){
  return(x[!x %in% y])
}
########################################################################

########################################################################
## This function converts the raw data to fragment level averages and ##
## sd of intensity. Each row of the output is a unique fragment. It   ##
## also includes the number of molecules aligned to those fragments.  ##
########################################################################
fn_aggregateByFrag <- function(Data=Data1){
  FragmentData1 <- aggregate(Data$moleculeID, by=list(Data$alignedChr, Data$alignedFragIndex), 
                             FUN=length)
  names(FragmentData1) <- c('alignedChr', 'alignedFragIndex', 'numMolecules')
  FragmentData.Median <- aggregate(Data[,c('fractionGC', 'intensityPerPixel')], 
                                   by=list(Data$alignedChr, Data$alignedFragIndex), 
                                   FUN=median)
  names(FragmentData.Median) <- c('alignedChr', 'alignedFragIndex', 'fractionGC_median', 
                                  'intensityPerPixel_median')
  FragmentData.Mean <- aggregate(Data[,c('fractionGC', 'intensityPerPixel')], 
                                 by=list(Data$alignedChr, Data$alignedFragIndex), 
                                 FUN=mean)
  names(FragmentData.Mean) <- c('alignedChr', 'alignedFragIndex', 'fractionGC_mean', 
                            'intensityPerPixel_mean')
  FragmentData.SD <- aggregate(Data[,c('fractionGC', 'intensityPerPixel')], 
                               by=list(Data$alignedChr, Data$alignedFragIndex), 
                               FUN=sd)
  names(FragmentData.SD) <- c('alignedChr', 'alignedFragIndex', 'fractionGC_sd', 
                            'intensityPerPixel_sd')
  FragmentData <- merge(x = merge(x = merge(x = FragmentData1, y = FragmentData.Mean), y = FragmentData.SD), y = FragmentData.Median)
  FragmentData[,c('fractionGC_sd', 'intensityPerPixel_sd')] <- sapply(X=FragmentData[,c('fractionGC_sd', 'intensityPerPixel_sd')], 
                                                                      FUN=na.is.zero)
  FragmentData <- FragmentData[order(FragmentData$numMolecules, decreasing=T),]  
  return(FragmentData)
}
########################################################################

########################################################################
## This function takes in the data created by the previous function & ##
## returns only the ones which have a certain number of coverage. The ##
## higher the coverage, the better the explanatory power of the model.##
########################################################################
fn_subsetByCoverage <- function(CoverageNum=15, filename='FragmentData_5000.RData'){
  Filename <- paste(FilePath, filename, sep='')
  load(Filename)
  
  FragmentData$fractionGC_sd <- NULL
  FragmentData$fractionGC_median <- NULL
  names(FragmentData)[names(FragmentData)=='fractionGC_mean'] <- 'fractionGC'
  
  FragmentData.Subset <- subset(FragmentData, numMolecules > CoverageNum)
  
  Fragment.Keep <- FragmentData.Subset[,c('alignedChr', 'alignedFragIndex', 'numMolecules', 
                                          'intensityPerPixel_sd')]
  Fragment.Keep$Keep <- 1
  
  Output <- list(Fragment.Keep=Fragment.Keep, FragmentData=FragmentData.Subset)
  return(Output)
}
########################################################################

########################################################################
## This function takes the summary() object of an lm output and retu- ##
## -rns a row of output in the desired format                         ##
########################################################################
fn_formatSummary <- function(Summary, Modelname){
  RSS <- round(Summary$sigma ^2, 4)
  RSquared <- round(Summary$r.squared, 4)
  Intercept_Beta <- round(Summary$coefficients['(Intercept)', 'Estimate'], 4)
  Intercept_TStat <- round(Summary$coefficients['(Intercept)', 't value'], 4)
  GC_Content_Beta <- round(Summary$coefficients['fractionGC', 'Estimate'], 4)
  GC_Content_StdError <- round(Summary$coefficients['fractionGC', 'Std. Error'], 4)
  GC_Content_TStat <- round(Summary$coefficients['fractionGC', 't value'], 4)

  Output <- cbind(RSS, RSquared, Intercept_Beta, Intercept_TStat, GC_Content_Beta, GC_Content_StdError, GC_Content_TStat)
  colnames(Output) <- paste(Modelname, colnames(Output), sep='_')
  rm(Summary)
  return(Output)
}
########################################################################

########################################################################
## This function takes identifies the molecule responsible for flagg- ##
## -ing the fragment location as an outlier and saves the subset of   ##
## the dataset with info on that fragment location, from all the mol- ##
## -ecules                                                            ##
########################################################################
fn_outlierMolecules <- function(RowIndex, OutlierFragments, Data, DataPath){
  Chr <- OutlierFragments[RowIndex, 'alignedChr']
  FragIndex <- OutlierFragments[RowIndex, 'alignedFragIndex']
  Subset <- subset(Data, alignedChr==Chr & alignedFragIndex==FragIndex)
  Subset$intensityPerPixel_median <- median(Subset$intensityPerPixel)
  #Subset$intensityPerPixel_median <- subset(OutlierFragments, alignedChr==Chr & alignedFragIndex==FragIndex)$intensityPerPixel_median
  Subset$Outlier <- Subset$intensityPerPixel > 2 * Subset$intensityPerPixel_median  ## First classification
  # Subset$Outlier[Subset$intensityPerPixel == max(Subset$intensityPerPixel)] <- TRUE    ## Second classification
  Subset$Outlier[Subset$intensityPerPixel > 75000] <- TRUE    ## Third classification (brute force)
  
  Colnames <- c('alignedChr', 'alignedFragIndex', 'moleculeID', 'intensityPerPixel', 'intensityPerPixel_sd', 'intensityPerPixel_median', 'Outlier')
  Subset <- Subset[,c(Colnames)]
  Filename <- paste(DataPath, 'Outliers/Subset_Chr', Chr, '_FragIndex', FragIndex, '.txt', sep='')
  write.table(Subset, file=Filename, sep='\t')
  if(sum(Subset$Outlier) > 0){
    return(subset(Subset, Outlier==TRUE))
  }
}
########################################################################

########################################################################
## This function eliminates the outliers identified in RScript06. It  ##
## takes in the big data and returns clean data.                      ##
########################################################################
fn_eliminateOutliers <- function(BigData, Outliers){
  BigData <- merge(x=BigData, y=Outliers[,c('alignedChr', 'alignedFragIndex', 'moleculeID', 'Outlier')], by=c('alignedChr', 'alignedFragIndex', 'moleculeID'), all.x=T, all.y=F)
  BigData <- subset(BigData, is.na(Outlier))
  BigData$Outlier <- NULL
  return(BigData)
}
########################################################################

########################################################################
## This function reads in Alignment chunk files created by Steve      ##
########################################################################
fn_readAlignmentChunks <- function(Filename, Header, Colnames, Colnames.Steve){
  AlChunk <- read.table(file=Filename, header=Header, sep='\t', 
                        skip=0, stringsAsFactors=FALSE)
  colnames(AlChunk) <- Colnames
  attributes(AlChunk)$comment <- Colnames.Steve
  return(AlChunk)
}

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
## This functions reads in the fragment data files produced by Steve' ##
## s code and extracts the intensity data and returns that in a data- ##
## -frame format. Each dataframe has a column of intensity values of  ##
## one molecule. That molecule name is in one of the other columns.   ##
########################################################################
fn_returnMoleculeIntensity <- function(MoleculeID, FragmentData){
  NColumns <- ncol(FragmentData)
  Subset <- subset(FragmentData, moleculeID==MoleculeID)
  Intensity <- as.vector(Subset[11:NColumns])
  Intensity <- Intensity[!is.na(Intensity)]
  IntensityData <- as.data.frame(Intensity)
  IntensityData$MoleculeID <- MoleculeID
  IntensityData$Intensity_Normalized <- IntensityData$Intensity/median(IntensityData$Intensity)
  IntensityData$PixelNum <- index(IntensityData)
  IntensityData <- within(data=IntensityData,{
    MoleculeID <- factor(MoleculeID)
  })
  return(IntensityData)
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
########################################################################

########################################################################
## This function takes in the Alignment chunk data and returns the n- ##
## -umber of molecules aligned to each fragments, by chromosome.      ##
########################################################################
fn_numMolAlignedperLoc <- function(AlChunk){
  Data <- unique(AlChunk[,c('refChr', 'refStartIndex', 'molID')])
  
  Chr_Frags_Mols.Table <- aggregate(x=Data$molID, by=list(Data$refChr, Data$refStartIndex), 
                                    FUN=length)
  names(Chr_Frags_Mols.Table) <- c('refChr', 'refStartIndex', 'numMolecules')
  return(Chr_Frags_Mols.Table)
}


########################################################################
## This function takes in a reference vector and a Test vector and    ##
## aligns them pixel by pixel, based on randomized pixel insertion/   ##
## deletion, a Monte-carlo type algorithm.                            ##
########################################################################
fn_alignPixels <- function(Test, Reference, Simulation_N=5000){
  if(length(Reference) < length(Test)){
    print("Reference Shorter")
    Diff <- length(Test) - length(Reference)
    SimulationLength <- min(Simulation_N, floor(choose(n=length(Test), k=Diff)/10))
    RSquared <- 0
    BestSample <- sample(x=1:length(Test), size=Diff)
    for(i in 1:SimulationLength){
      Sample <- sample(x=1:length(Test), size=Diff)
      Test1 <- Test[index(Test) %w/o% c(Sample)]
      Model <- lm(Test1 ~ Reference)
      if(summary(Model)$r.squared > RSquared){
        RSquared <- summary(Model)$r.squared
        BestSample <- Sample
      }
      #print(c(RSquared, summary(Model)$r.squared))
    }
    Test.Aligned <- Test[index(Test) %w/o% c(BestSample)]
    print(RSquared)
    return(Test.Aligned)
  } else if(length(Test) < length(Reference)){
    print("Reference Longer")
    Diff <- length(Reference) - length(Test)
    SimulationLength <- min(Simulation_N, floor(choose(n=length(Test), k=Diff)/10))
    RSquared <- 0
    Test.Aligned <- Test
    for(i in 1:SimulationLength){
#       print(paste(i, RSquared))
      Test1 <- Test
      for(PixelNum in 1:Diff){
        Insert_PixelNum <- sample(x=2:(length(Test1) -1), size=1)
        Insert_Intensity <- 0.5*(Test1[Insert_PixelNum] + Test1[Insert_PixelNum + 1])
        Test1 <- c(Test1[1:Insert_PixelNum], Insert_Intensity, 
                   Test1[(Insert_PixelNum + 1):length(Test1)])
      }
      Model <- lm(Test1 ~ Reference)
      if(summary(Model)$r.squared > RSquared){
        RSquared <- summary(Model)$r.squared
        Test.Aligned <- Test1
      }
    }
    print(RSquared)
    return(Test.Aligned)
  } else{
    print("Same Length")
    Test.Aligned <- Test
    return(Test.Aligned)
  }
}

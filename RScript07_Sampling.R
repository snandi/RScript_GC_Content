rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script is for conducting further analysis on intensity and gc ##
## content. This script includes a molecule level effect and randomly ##
## samples few molecules at a time, to avoid big data problems        ##
########################################################################

########################################################################
## Initialize Header file and function library                        ##
########################################################################
FilePath <- '~/RA_Genomics/Project_GC_Content/RScripts_GC_Content/'
DataPath <- '~/RA_Genomics/Project_GC_Content/'
DataPath.Steve <- DataPath
Filename.Header <- paste('~/RScripts/HeaderFile_Nandi.R', sep='')
source(Filename.Header)
source(paste(FilePath, 'fn_Library_GC_Content.R', sep=''))
########################################################################

########################################################################
## Data Input                                                         ##
########################################################################
Fragments <- fn_subsetByCoverage(CoverageNum=3, filename='FragmentData_7134.RData')
names(Fragments)

Fragment.Keep <- Fragments[['Fragment.Keep']]
FragmentData <- Fragments[['FragmentData']]

Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.7134Groups', sep='')
#Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.5000Groups', sep='')
#Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.1061Groups', sep='')
#Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.50Groups', sep='')
Data <- read.table(file=Filename, header=TRUE, sep='\t', 
                   skip=0, stringsAsFactors=FALSE, fill=TRUE)
Data$groupID <- substr(x=Data$moleculeID, start=1, stop=7)
str(Data)
dim(Data)

Data <- merge(x=Data, y=Fragment.Keep, all=T)
Data <- subset(Data, Keep==1)
Data$Weight1 <- Data$numMolecules
Data$Weight2 <- 1/Data$intensityPerPixel_sd

rm(Fragments)           ## Clean memory

Filename1 <- paste(DataPath, 'Outliers/OutlierFragments.RData', sep='')
load(Filename1)
OutlierFragments <- OutlierFragments[order(OutlierFragments$alignedChr, OutlierFragments$alignedFragIndex),]

Filename2 <- paste(DataPath, 'Outliers/OutlierMolecules.RData', sep='')
load(Filename2)

Data.Clean <- fn_eliminateOutliers(BigData=Data, Outliers=OutlierMolecules)
dim(Data.Clean)

Data.Clean <- subset(Data.Clean, intensityPerPixel < 75000)
Data.Clean$fragSize <- Data.Clean$numPixels * Data.Clean$fragConversionFactor

summary(Data.Clean$fractionGC)
summary(Data.Clean$intensityPerPixel)

#######################################################################
## Checking if intensityPerPixel_sd has any relationship with         ##
## numMolecules. This is for choosing the right weights.              ##
## DO NOT DELETE THIS SECTION!                                        ##
########################################################################
# Var_vs_Mol <- lm(log(intensityPerPixel_sd) ~ numMolecules, data=FragmentData)
# plot(log(intensityPerPixel_sd) ~ numMolecules, data=FragmentData, pch='.')
# lines(lowess(y=log(FragmentData$intensityPerPixel_sd), x=FragmentData$numMolecules), col="Red")
# 
# layout(matrix(1:4, ncol=2, byrow=TRUE))
# plot(Var_vs_Mol)
# layout(matrix(1:1))
## Conclusion: No pattern observed between variance of intensity and
## numMolecules.
########################################################################

Groups <- unique(Data.Clean$groupID)

## sampleGroups <- sample(x=Groups, size=15, replace=FALSE, prob=NULL)
## sampleData <-subset(Data.Clean, groupID %in% sampleGroups)

## length(unique(sampleData$moleculeID))

## Time1 <- Sys.time()
## ModelA <- lm( log(intensityPerPixel) ~ fractionGC + fragConversionFactor + as.factor(moleculeID), data=sampleData)

## ModelB <- lm( log(intensityPerPixel) ~ fractionGC + fragConversionFactor + as.factor(moleculeID), data=sampleData, weights=Weight1)

## ModelC <- lm( log(intensityPerPixel) ~ fractionGC + fragConversionFactor + as.factor(moleculeID), data=sampleData, weights=Weight2)

## summary(ModelA)$coefficients[1:3,]
## summary(ModelA)$r.squared

## summary(ModelB)$coefficients[1:4,]
## summary(ModelB)$r.squared

## summary(ModelC)$coefficients[1:4,]
## summary(ModelC)$r.squared
## Time2 <- Sys.time()
## print(Time2 - Time1)

# Colnames of output:
## Sample number
## Num of fragments
## Num of molecules
## ModelA
## ## RSS
## ## DF
## ## RSquared
## ## Intercept Estimate
## ## Intercept TStat
## ## fractionGC Estimate
## ## fractionGC TStat

## cbind(fn_formatSummary(Summary=summary(ModelA), Modelname='ModelA'),
##       fn_formatSummary(Summary=summary(ModelB), Modelname='ModelB'),
##       fn_formatSummary(Summary=summary(ModelC), Modelname='ModelC'),
##       numMolecules=length(unique(sampleData$moleculeID)))

fn_modelSample <- function(Data.Clean, Groups, SampleSize){
  sampleGroups <- sample(x=Groups, size=SampleSize, replace=FALSE, prob=NULL)
  sampleData <-subset(Data.Clean, groupID %in% sampleGroups)

  ModelA <- lm( log(intensityPerPixel) ~ fractionGC + fragConversionFactor + as.factor(moleculeID), data=sampleData)
  
  ModelB <- lm( log(intensityPerPixel) ~ fractionGC + fragConversionFactor + as.factor(moleculeID), data=sampleData, weights=Weight1)

  ModelC <- lm( log(intensityPerPixel) ~ fractionGC + fragConversionFactor + as.factor(moleculeID), data=sampleData, weights=Weight2)
  
  OutputRow <- cbind(fn_formatSummary(Summary=summary(ModelA), Modelname='ModelA'),
                     fn_formatSummary(Summary=summary(ModelB), Modelname='ModelB'),
                     fn_formatSummary(Summary=summary(ModelC), Modelname='ModelC'),
                     numMolecules=length(unique(sampleData$moleculeID)))
  return(OutputRow)
}

## Time1 <- Sys.time()
## MC_Output <- fn_modelSample(Data.Clean=Data.Clean, Groups=Groups, SampleSize=15)
## Time2 <- Sys.time()
## print(Time2 - Time1)

## for(Index in 1:4){
##   MC_Output <- rbind(MC_Output, fn_modelSample(Data.Clean=Data.Clean, Groups=Groups, SampleSize=15))
## }

########################################################################
## This part of the code is parallelized and should be run on bigmems ##
## only.                                                              ##
########################################################################
NumCores <- 8   ## Number of cores to be used
cl <- makeSOCKcluster(as.numeric(NumCores))
registerDoSNOW(cl)

MC_Output <- foreach(Indices=1:4500, .inorder=FALSE, .packages=MyAutoLoads, .combine=rbind) %dopar% fn_modelSample(Data.Clean=Data.Clean, Groups=Groups, SampleSize=15)

stopCluster(cl)
rm(cl)
########################################################################

Filename1 <- paste(DataPath, 'SamplingOutput.RData', sep='')
save(MC_Output, file=Filename1)

Filename2 <- paste(DataPath, 'SamplingOutput.txt', sep='')
write.table(MC_Output, file=Filename2, sep='\t')

MC_Output <- as.data.frame(MC_Output)
hist(MC_Output$ModelA_RSquared, breaks=40)
hist(MC_Output$ModelB_RSquared, breaks=40)

hist(MC_Output$ModelA_GC_Content_Beta, breaks=40)
hist(MC_Output$ModelC_GC_Content_Beta, breaks=40)

rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script is for conducting further analysis on intensity and gc ##
## content. This will use data produced by Steve on 07/10 and after.  ##
## This script will incorporate fragment level effect first, as per   ##
## discussion during the meeting with M. Newton and D. Schwartz on    ##
## 07/15                                                              ##
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
Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.7134Groups', sep='')
#Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.5000Groups', sep='')
#Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.1061Groups', sep='')
#Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.50Groups', sep='')
Data <- read.table(file=Filename, header=TRUE, sep='\t', 
                   skip=0, stringsAsFactors=FALSE, fill=TRUE)
str(Data)
dim(Data)
#View(head(Data))
#View(Data)
#Data <- Data[1:100000,]
########################################################################

########################################################################
## Some Molecule level summary                                        ##
########################################################################
length(unique(Data$moleculeID))
Data.Molecule <- unique(Data[,c('moleculeID', 'moleculeConversionFactor', 'numFrags', 
                                'numPixels', 'totalIntensity')])
length(unique(Data.Molecule$moleculeID))
Corr.Molecule <- corr.test(Data.Molecule[,-1])


#View(subset(Data, fractionGC==0.5834))

Data$groupID <- substr(x=Data$moleculeID, start=1, stop=7)
length(unique(Data$moleculeID))
length(unique(Data$groupID))

########################################################################
## Some Fragment level summary                                        ##
########################################################################
AlignedChr <- unique(Data$alignedChr)

for(Chr in AlignedChr){
  print(Chr)
  Data1 <- subset(Data, alignedChr==Chr)
  if('FragmentData' %in% ls()){
      FragmentData <- rbind(FragmentData, fn_aggregateByFrag(Data = Data1))
  } else{
    FragmentData <- fn_aggregateByFrag(Data = Data1)
  }
}
head(FragmentData)
dim(FragmentData)
summary(FragmentData)

Filename <- paste(FilePath, 'FragmentData_7134.RData', sep='')
#Filename <- paste(FilePath, 'FragmentData_5000.RData', sep='')
#Filename <- paste(FilePath, 'FragmentData_1061.RData', sep='')
save(FragmentData, file=Filename)

########################################################################
## Regression Models                                                  ##
########################################################################
Model1 <- lm(log(intensityPerPixel_mean) ~ fractionGC_mean, data=FragmentData)
summary(Model1)
#plot(Model1)

Model2 <- lm(log(intensityPerPixel_sd) ~ fractionGC_mean, data=FragmentData)
summary(Model2)

plot(log(intensityPerPixel_mean) ~ fractionGC_mean, data=FragmentData)


Model3 <- lm(log(intensityPerPixel_mean) ~ fractionGC_mean + as.factor(alignedChr), 
             data=FragmentData)
summary(Model3)
layout(matrix(1:4, ncol=2, byrow=TRUE))
plot(Model3)
layout(matrix(1:1))

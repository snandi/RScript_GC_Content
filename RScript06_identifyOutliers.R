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

#Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.7134Groups', sep='')
Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.5000Groups', sep='')
#Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.1061Groups', sep='')
#Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.50Groups', sep='')
Data <- read.table(file=Filename, header=TRUE, sep='\t', 
                   skip=0, stringsAsFactors=FALSE, fill=TRUE)
Data$groupID <- substr(x=Data$moleculeID, start=1, stop=7)
str(Data)
dim(Data)

Data <- merge(x=Data, y=Fragment.Keep, all=T)
Data <- subset(Data, Keep==1)

FragmentData$MeanMedian_Diff <- round(abs(FragmentData$intensityPerPixel_median - FragmentData$intensityPerPixel_mean)/FragmentData$intensityPerPixel_median, 3)
FragmentData$Outlier <- (FragmentData$MeanMedian_Diff >= 0.10)
OutlierFragments <- subset(FragmentData, Outlier==T)

rm(Fragments)           ## Clean memory

#View(subset(Data, moleculeID=='2388301_734_2060415'))

#hist(FragmentData$intensityPerPixel_sd, breaks=40)
#summary(FragmentData$intensityPerPixel_sd)
## Chromosome: 6, Location Index: 7753, NumMolecules: 8, 0.4128	237098.27	558403.9
#View(subset(Data, alignedChr==6 & alignedFragIndex==7753))

## fn_outlierMolecules(RowIndex=5, OutlierFragments=OutlierFragments, Data=Data, DataPath=DataPath)

NumCores <- 12   ## Number of cores to be used
cl <- makeSOCKcluster(as.numeric(NumCores))
registerDoSNOW(cl)

OutlierMolecules <- foreach(Rows=1:nrow(OutlierFragments), .inorder=FALSE, .packages=MyAutoLoads, .combine=rbind) %dopar% fn_outlierMolecules(RowIndex=Rows, OutlierFragments=OutlierFragments, Data=Data, DataPath=DataPath)

stopCluster(cl)
rm(cl)

Filename1 <- paste(DataPath, 'Outliers/OutlierFragments.RData', sep='')
save(OutlierFragments, file=Filename1)

Filename2 <- paste(DataPath, 'Outliers/OutlierMolecules.RData', sep='')
save(OutlierMolecules, file=Filename2)

Filename3 <- paste(DataPath, 'Outliers/OutlierMolecules.txt', sep='')
write.table(OutlierMolecules, file=Filename3, sep='\t')


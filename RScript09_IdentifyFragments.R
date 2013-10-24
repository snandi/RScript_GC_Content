rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script is just for identifying fragments for Steve, that have ##
## distinctive gc content patterns.                                   ##
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

HighGC <- subset(FragmentData, fractionGC > 0.75)
LowGC <- subset(FragmentData, fractionGC < 0.25)

Filename1 <- paste(DataPath, 'Fragments_HighGC.txt', sep='')
write.table(HighGC, file=Filename1, sep='\t')

Filename2 <- paste(DataPath, 'Fragments_LowGC.txt', sep='')
write.table(LowGC, file=Filename2, sep='\t')

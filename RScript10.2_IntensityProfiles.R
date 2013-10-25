rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script is for plotting intensity profiles of multiple molecu- ##
########################################################################

########################################################################
## Initialize Header file and function library                        ##
########################################################################
FilePath <- '~/Project_GC_Content/RScripts_GC_Content/'
DataPath <- '/exports/aspen/steveg/human_nMaps/GC_content/subdivideFragments/testDir/'
Filename.Header <- paste('~/RScripts/HeaderFile_lmcg.R', sep='')
source(Filename.Header)
source(paste(FilePath, 'fn_Library_GC_Content.R', sep=''))
########################################################################

setwd(DataPath)

AllFragmentNames <- list.files()[grep(pattern='intensities', x=list.files())]

AllFragmentNames <- AllFragmentNames %w/o% c('chr3.bp50713000.intensities', 
                                             'chr11.frag4854.intensities')
FragmentName <- 'chr10.bp107002000.intensities'
FragmentName.Sub <- substr(x=FragmentName, start=1, stop=(nchar(FragmentName)-12))
FragmentName.GC <- gsub(pattern='intensities', replacement='gcContent', x=FragmentName)

FragmentFilename <- paste(DataPath, FragmentName, sep='')
## These fragment files are created by Steve
FragmentData <- read.table(file=FragmentFilename, header=TRUE, 
                           stringsAsFactors=FALSE)
GCFilename <- paste(DataPath, FragmentName.GC, sep='')
GCFile <- read.table(file=GCFilename, header=FALSE, sep=',', skip=2,
                     stringsAsFactors=FALSE)
GCData <- as.data.frame(cbind(GC=GCFile[,5],
                              BP100=index(GCFile)))
#View(FragmentData)
#str(FragmentData)
Molecules <- FragmentData[,'moleculeID']
MoleculeID <- Molecules[1]
NumPages <- ceil(length(Molecules)/6) ## Depends on the number of molecules aligned

IntensityData <- do.call(what=rbind, lapply(X=Molecules, FUN=fn_returnMoleculeIntensity, 
                                            FragmentData=FragmentData))
View(IntensityData)

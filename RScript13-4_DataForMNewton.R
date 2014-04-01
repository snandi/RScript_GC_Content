rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script plots the intensity profiles of mFlorum intervals      ##
## along with the acgt content 
########################################################################

########################################################################
## Initialize Header file and function library                        ##
########################################################################
OutputDataPath <- '~/Project_GC_Content/RData/mflorum/'
RScriptPath <- '~/Project_GC_Content/RScripts_GC_Content/'
RDataPath <- '~/Project_GC_Content/RData/mflorum/'
DataPath.MFlorum <- '~/mflorum_nMaps/GC_Content/'
Filename.Header <- paste('~/RScripts/HeaderFile_lmcg.R', sep='')
source(Filename.Header)
source(paste(RScriptPath, 'fn_Library_GC_Content.R', sep=''))
########################################################################

########################################################################
## Load the binary RData file of the Alignment Chunks, previously cr- ##
## eated                                                              ##
########################################################################
Filename.Bin <- paste(OutputDataPath, 
                      'MF_cap348_inca34_cf209_minSize50_minFrag5_alignmentChunks.RData', sep='')
load(Filename.Bin)

Filename.bploc <- paste(RDataPath, 'mflorum.bploc', sep='')
bp.loc <- read.table(file=Filename.bploc, header=T, sep='\t')

########################################################################
## Enter Fragment Details                                             ##
########################################################################
Chr <- 'chr1'
FragIndex <- 25
BasePairInterval <- 209     ## Length of base pair interval to estimate gcat %      
## Starting base pair coordinate of fragment
FragBP_Start <- bp.loc[which(bp.loc$alignedFragIndex == FragIndex), 'refMapCoordStart'] 
## Ending base pair coordinate of fragment
FragBP_End <- bp.loc[which(bp.loc$alignedFragIndex == FragIndex), 'refMapCoordEnd']

NumBP_Frag <- FragBP_End - FragBP_Start ## Length of frag in BP
NumSubFrag <- round(NumBP_Frag/BasePairInterval, 0) ## Number of sub fragments
## Also, the number of pixels

#for(FragIndex in 1:39){


#########################################################################
## For All molecules aligned to the fragment                           ##
#########################################################################
FragmentName <- paste(Chr, '_frag', FragIndex, '_intensities', sep='')
FragmentName.Sub <- substr(x=FragmentName, start=1, stop=(nchar(FragmentName)-12))
FragmentName.GC <- gsub(pattern='intensities', replacement='gcContent', x=FragmentName)
FragmentFilename <- paste(RDataPath, FragmentName, sep='')
## These fragment files are created by Steve
FragmentData <- read.table(file=FragmentFilename, header=TRUE, 
                           stringsAsFactors=FALSE)
Molecules <- FragmentData[,'moleculeID']
MoleculeID <- Molecules[1]
NumPages <- ceil(length(Molecules)/6) ## Depends on the number of molecules aligned

IntensityData <- do.call(what=rbind, lapply(X=Molecules, FUN=fn_returnMoleculeIntensity, 
                                            FragmentData=FragmentData))

NumPixels <- aggregate(IntensityData$PixelNum, by=list(IntensityData$MoleculeID), FUN=max)
colnames(NumPixels) <- c('MoleculeID', 'Pixels')
table(NumPixels$Pixels)

#########################################################################
## Align all molecules by pixel numbers                                ##
#########################################################################
Pixel.Mode <- Mode(NumPixels$Pixels)
FragmentLength <- Pixel.Mode

IntensityData.199 <- subset(IntensityData, MoleculeID %in% 
                                  subset(NumPixels, Pixels==FragmentLength)[,'MoleculeID'])
IntensityData.199.Wide <- reshape(IntensityData.199[,-3], 
                                      timevar='MoleculeID', 
                                      idvar='PixelNum', 
                                      direction='wide')

MFlorumData_Interval25 <- list(NumPixels=NumPixels, 
                               IntensityData=IntensityData, 
                               IntensityData.199=IntensityData.199, 
                               IntensityData.199.Wide=IntensityData.199.Wide)

save(MFlorumData_Interval25, file=paste(RDataPath, 
                                        'MFlorumData_Interval25.RData', sep=''))


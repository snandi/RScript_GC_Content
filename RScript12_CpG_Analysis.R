rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script is for reading and extracting important information    ##
## out of any alignmentchunk file, produces after an alignment        ##
########################################################################

########################################################################
## Run Path definition file                                           ##
########################################################################
source('~/Project_GC_Content/RScripts_GC_Content/Paths_Header.R')
########################################################################

########################################################################
## Load the binary RData file of the Alignment Chunks, previously cr- ##
## eated                                                              ##
########################################################################
Filename.Bin <- paste(OutputDataPath, 
                      'alignmentChunks.withLength.all7134Groups.goldOnly.RData', sep='')
load(Filename.Bin)

########################################################################
## Chr3, reference Fragment Index: 13762                              ##
########################################################################
FragmentName <- 'chr3_frag13762_intensities'
FragmentName.Sub <- substr(x=FragmentName, start=1, stop=(nchar(FragmentName)-12))
FragmentName.GC <- gsub(pattern='intensities', replacement='gcContent', x=FragmentName)

FragmentFilename <- paste(DataPath.Nandi, FragmentName, sep='')
## These fragment files are created by Steve
FragmentData <- read.table(file=FragmentFilename, header=TRUE, 
                           stringsAsFactors=FALSE)
Molecules <- FragmentData[,'moleculeID']
MoleculeID <- Molecules[1]
NumPages <- ceil(length(Molecules)/6) ## Depends on the number of molecules aligned

#########################################################################
## For All molecules aligned to the fragment                           ##
#########################################################################
IntensityData <- do.call(what=rbind, lapply(X=Molecules, FUN=fn_returnMoleculeIntensity, 
                                            FragmentData=FragmentData))

NumPixels <- aggregate(IntensityData$PixelNum, by=list(IntensityData$MoleculeID), FUN=max)
NumPixels

Reference <- subset(IntensityData, MoleculeID=="2389364_734_2060253")[,'Intensity']
IntensityData.Aligned <- subset(IntensityData, MoleculeID=="2389364_734_2060253")
# Test <- subset(IntensityData, MoleculeID=="2398367_734_2060070")[,'Intensity']
# Test.Aligned <- fn_alignPixels(Test=Test, Reference=Reference)

for(Molecule in (levels(IntensityData$MoleculeID) %w/o% "2389364_734_2060253")){
  print(Molecule)
  Test <- subset(IntensityData, MoleculeID==Molecule)[,'Intensity']
  print(length(Test))
  Test.Aligned <- fn_alignPixels(Test=Test, Reference=Reference)
  Intensity <- Test.Aligned
  TestData <- as.data.frame(Intensity)
  TestData$MoleculeID <- Molecule
  TestData$Intensity_Normalized <- TestData$Intensity/median(TestData$Intensity)
  TestData$PixelNum <- index(TestData)
  TestData <- within(data=TestData,{
    MoleculeID <- factor(MoleculeID)
  })
  IntensityData.Aligned <- rbind(IntensityData.Aligned, TestData)
}

IntensityData.Aligned.Wide <- reshape(IntensityData.Aligned[,-1], 
                                      timevar='MoleculeID', 
                                      idvar='PixelNum', 
                                      direction='wide')
IntensityData.Aligned.Wide$Intensity_Mean <- rowMeans(IntensityData.Aligned.Wide[,2:21])
IntensityData.Aligned.Wide$Intensity_SD <- apply(X=IntensityData.Aligned.Wide[,2:21], 
                                                 MARGIN=1, FUN=sd)
IntensityData.Aligned.Wide$Intensity_Up <- IntensityData.Aligned.Wide$Intensity_Mean + 3*IntensityData.Aligned.Wide$Intensity_SD
IntensityData.Aligned.Wide$Intensity_Dn <- IntensityData.Aligned.Wide$Intensity_Mean - 3*IntensityData.Aligned.Wide$Intensity_SD




############################ PLOTTING ################################
Molecules <- levels(IntensityData$MoleculeID)
NumPages <- ceil(length(Molecules)/6) ## Depends on the number of molecules aligned

Today <- Sys.Date()
Filename.pdf <- paste('~/Project_GC_Content/RScripts_GC_Content/Plots/IntensityPlots_', 
                      FragmentName.Sub, '_', Today, '.pdf', sep='')

pdf(file=Filename.pdf, onefile=TRUE, pointsize=6)
trellis.par.set(fontsize=list(text=8,points=8))

Page <- 1
for(Page in 1:NumPages){
  StartMolecule <- 6*(Page - 1) + 1
  EndMolecule <- min((StartMolecule + 5), length(Molecules))
  NumPanels <- EndMolecule - StartMolecule + 1
  
#   PlotGC <- lattice::xyplot(GC ~ BP100, 
#                             data=GCData, 
#                             scales=list(x="free", y="free"), 
#                             type=c("l", "g"), lwd=1.25, col='gray87', 
#                             panel = function(...) { 
#                               panel.fill(col = 'royalblue4') 
#                               panel.xyplot(...) 
#                             }, 
#                             xlab='',
#                             main='GC Content of 100bp windows')
  
  PlotIntensity <- lattice::xyplot(Intensity_Normalized ~ PixelNum | MoleculeID, 
                                   data=subset(IntensityData, MoleculeID %in% levels(IntensityData$MoleculeID)[StartMolecule:EndMolecule]), 
                                   layout=c(1,6), 
                                   scales=list(x="free", y="free"), 
                                   type=c("l", "g"), lwd=1.25, col='darkblue',
                                   panel = function(...) { 
                                     panel.fill(col = 'gray78') 
                                     panel.xyplot(...) 
                                   }, 
                                   main=paste('Intensity plots of molecules aligned to', FragmentName.Sub))
#   grid.arrange(PlotGC, PlotIntensity, ncol=1, heights=c(1/5,4/5))
  print(PlotIntensity)
  ## This grid command prints the two charts into the pdf device
}
dev.off()

Filename.pdf <- paste('~/Project_GC_Content/RScripts_GC_Content/Plots/IntensityPlots_', 
                      FragmentName.Sub, '_Aligned_', Today, '.pdf', sep='')

pdf(file=Filename.pdf, onefile=TRUE, pointsize=6)
trellis.par.set(fontsize=list(text=8,points=8))

Page <- 1
for(Page in 1:NumPages){
  StartMolecule <- 6*(Page - 1) + 1
  EndMolecule <- min((StartMolecule + 5), length(Molecules))
  NumPanels <- EndMolecule - StartMolecule + 1
  
#   PlotGC <- lattice::xyplot(GC ~ BP100, 
#                             data=GCData, 
#                             scales=list(x="free", y="free"), 
#                             type=c("l", "g"), lwd=1.25, col='gray87', 
#                             panel = function(...) { 
#                               panel.fill(col = 'darkgreen') 
#                               panel.xyplot(...) 
#                             }, 
#                             xlab='',
#                             main='GC Content of 100bp windows')
  
  PlotIntensity <- lattice::xyplot(Intensity_Normalized ~ PixelNum | MoleculeID, 
                                   data=subset(IntensityData.Aligned, MoleculeID %in% levels(IntensityData.Aligned$MoleculeID)[StartMolecule:EndMolecule]), 
                                   layout=c(1,6), 
                                   scales=list(x="free", y="free"), 
                                   type=c("l", "g"), lwd=1.25, col='darkgreen',
                                   panel = function(...) { 
                                     panel.fill(col = 'gray78') 
                                     panel.xyplot(...) 
                                   }, 
                                   main=paste('After aligning all the pixels to a reference', FragmentName.Sub))
#   grid.arrange(PlotGC, PlotIntensity, ncol=1, heights=c(1/5,4/5))
  print(PlotIntensity)
  ## This grid command prints the two charts into the pdf device
}

PlotMeanIntensity <- lattice::xyplot(Intensity_Mean + Intensity_Up + Intensity_Dn ~ PixelNum, 
                                     data=subset(IntensityData.Aligned.Wide), 
                                     layout=c(1,1), 
                                     scales=list(x="free", y="free"), 
                                     type=c("l", "g"), col=c('darkgreen'),
                                     lty=c(1, 3, 3), lwd=c(1.5, 1.5, 1.5),
                                     panel = function(...) { 
                                       panel.fill(col = 'gray78') 
                                       panel.xyplot(...) 
                                     }, 
                                     main=paste('Mean Intensity plot after pixel alignment'), 
                                     ylab='Mean Intensity +/- 3 standard dev')
# grid.arrange(PlotGC, PlotMeanIntensity, ncol=1, heights=c(1/3,2/3))
print(PlotMeanIntensity)

dev.off()







AlChunk.chr7 <- subset(AlChunk, refChr=='chr7')
AlChunk.chr7$num_basepairs <- AlChunk.chr7$refEndCoord - AlChunk.chr7$refStartCoord

#head(AlChunk.chr7)
#tail(AlChunk.chr7)

1706711  1713935

StartBPCoord <- 1706711
EndBPCoord <- 1713935
View(subset(AlChunk.chr7, refStartCoord >= StartBPCoord & refEndCoord <= EndBPCoord))
View(subset(AlChunk.chr7, refStartCoord >= StartBPCoord))

str(AlChunk.chr7)
head(AlChunk.chr7)
tail(AlChunk.chr7)

rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script is for plotting intensity profiles of multiple molecu- ##
## -les aligned to the same fragment location on the genome. This file##
## reads in file created by Steve                                     ##
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

########################################################################
## Bigger for loop, one reference fragment location at a time         ##
########################################################################
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

#########################################################################
## For One molecule                                                    ##
#########################################################################
# IntensityData <- fn_returnMoleculeIntensity(MoleculeID=Molecules[1],FragmentData=FragmentData)
# str(IntensityData)

#########################################################################
## For All molecules aligned to the fragment                           ##
#########################################################################
IntensityData <- do.call(what=rbind, lapply(X=Molecules, FUN=fn_returnMoleculeIntensity, 
                                            FragmentData=FragmentData))
#dim(IntensityData)
#head(IntensityData)
IntensityData <- subset(IntensityData, MoleculeID != "2391387_734_2060372")

NumPixels <- aggregate(IntensityData$PixelNum, by=list(IntensityData$MoleculeID), FUN=max)
summary(NumPixels)

Reference <- subset(IntensityData, MoleculeID=="2393713_734_2060200")[,'Intensity']
IntensityData.Aligned <- subset(IntensityData, MoleculeID=="2393713_734_2060200")

for(MoleculeID in levels(IntensityData$MoleculeID) %w/o% "2393713_734_2060200"){
  print(MoleculeID)
  Test <- subset(IntensityData, MoleculeID==MoleculeID)[,'Intensity']
  Test.Aligned <- fn_alignPixels(Test=Test, Reference=Reference)
  Intensity <- Test.Aligned
  TestData <- as.data.frame(Intensity)
  TestData$MoleculeID <- MoleculeID
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
  
  PlotGC <- lattice::xyplot(GC ~ BP100, 
                            data=GCData, 
                            scales=list(x="free", y="free"), 
                            type=c("l", "g"), lwd=1.25, col='gray87', 
                            panel = function(...) { 
                              panel.fill(col = 'royalblue4') 
                              panel.xyplot(...) 
                            }, 
                            xlab='',
                            main='GC Content of 100bp windows')
  
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
  grid.arrange(PlotGC, PlotIntensity, ncol=1, heights=c(1/5,4/5))
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
  
  PlotGC <- lattice::xyplot(GC ~ BP100, 
                            data=GCData, 
                            scales=list(x="free", y="free"), 
                            type=c("l", "g"), lwd=1.25, col='gray87', 
                            panel = function(...) { 
                              panel.fill(col = 'darkgreen') 
                              panel.xyplot(...) 
                            }, 
                            xlab='',
                            main='GC Content of 100bp windows')
  
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
  grid.arrange(PlotGC, PlotIntensity, ncol=1, heights=c(1/5,4/5))
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
grid.arrange(PlotGC, PlotMeanIntensity, ncol=1, heights=c(1/3,2/3))

dev.off()


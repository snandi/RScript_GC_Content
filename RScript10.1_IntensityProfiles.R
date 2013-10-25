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
for(FragmentName in AllFragmentNames){
  #FragmentName <- 'chr10.bp107002000.intensities'
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
  
  Today <- Sys.Date()
  Filename.pdf <- paste('~/Project_GC_Content/Adhoc_2013-09-25/IntensityPlots_', FragmentName.Sub, 
                        '_', Today, '.pdf', sep='')
  
  pdf(file=Filename.pdf, onefile=TRUE, pointsize=6)
  trellis.par.set(fontsize=list(text=8,points=8))
  ## trellis is used to set the font size back to smaller size as pdf resets the fontsizes
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
  rm(IntensityData, Molecules)
  print(FragmentName.Sub)
}

rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script looks at the CpG islands in the Alpha-globin sites that##
## has a very large CpG island that is never methylated. This is      ##
## essentially the same script as RScript_12. Eventually, this will   ##
## be deprecated                                                      ##
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

########################################################################
## Enter Fragment Details                                             ##
########################################################################
Chr <- 'chr1'
FragIndex <- 35
BasePairInterval <- 209     ## Length of base pair interval to estimate gcat %      
# FragBP_Start <- 461077      ## Starting base pair coordinate of fragment
# FragBP_End <- 493266        ## Ending base pair coordinate of fragment
# NumBP_Frag <- FragBP_End - FragBP_Start ## Length of frag in BP
# NumSubFrag <- NumBP_Frag/BasePairInterval ## Number of sub fragments
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

IntensityData.Aligned <- subset(IntensityData, MoleculeID %in% 
                                  subset(NumPixels, Pixels==Pixel.Mode)[,'MoleculeID'])
IntensityData.Aligned.Wide <- reshape(IntensityData.Aligned[,-3], 
                                      timevar='MoleculeID', 
                                      idvar='PixelNum', 
                                      direction='wide')
NumCol <- ncol(IntensityData.Aligned.Wide)
IntensityData.Aligned.Wide$Intensity_Mean <- rowMeans(IntensityData.Aligned.Wide[,2:NumCol])
IntensityData.Aligned.Wide$Intensity_SD <- apply(X=IntensityData.Aligned.Wide[,2:NumCol], 
                                                 MARGIN=1, FUN=sd)
IntensityData.Aligned.Wide$Intensity_Up <- IntensityData.Aligned.Wide$Intensity_Mean + 1*IntensityData.Aligned.Wide$Intensity_SD
IntensityData.Aligned.Wide$Intensity_Dn <- IntensityData.Aligned.Wide$Intensity_Mean - 1*IntensityData.Aligned.Wide$Intensity_SD

#View(IntensityData.Aligned.Wide)

############################### ANOVA ################################
IntensityData.Aligned$PixelFactor <- as.factor(IntensityData.Aligned$PixelNum)
str(IntensityData.Aligned)
Model1 <- lm(Intensity ~ PixelFactor, data=IntensityData.Aligned)
anova(Model1)
summary(Model1)

Model2 <- lm(Intensity ~ MoleculeID, data=IntensityData.Aligned)
anova(Model2)
summary(Model2)

Model3 <- lm(Intensity ~ MoleculeID + PixelFactor, data=IntensityData.Aligned)
anova(Model3)
anova(Model2, Model3)
summary(Model3)

############################ PLOTTING ################################
IntensityData.Aligned$MoleculeID <- as.vector(IntensityData.Aligned$MoleculeID)
Molecules <- unique(IntensityData.Aligned$MoleculeID)
IntensityData.Aligned$MoleculeID <- as.factor(IntensityData.Aligned$MoleculeID)
Panels <- 5
NumPages <- ceil(length(Molecules)/Panels) ## Depends on the number of molecules aligned

Today <- Sys.Date()
Page <- 1
Filename.pdf <- paste('~/Project_GC_Content/RScripts_GC_Content/Plots/MF_', 
                      FragmentName.Sub, '_Pixels', nrow(IntensityData.Aligned.Wide), 
                      '_', Today, '.pdf', sep='')

pdf(file=Filename.pdf, onefile=TRUE, pointsize=6)
trellis.par.set(fontsize=list(text=8,points=8))

PlotHist <- lattice::barchart(x=table(NumPixels$Pixels), col='gray50',
                              xlab='Number of nMaps', ylab='Pixel length of nMaps',
                              main=paste('Histogram of nMap Pixels aligned to', FragmentName.Sub))
print(PlotHist)

for(Page in 1:NumPages){
  StartMolecule <- Panels*(Page - 1) + 1
  EndMolecule <- min((StartMolecule + Panels - 1), length(Molecules))
  NumPanels <- EndMolecule - StartMolecule + 1
  
  my.Colors=c('olivedrab4', 'olivedrab1', 'gray25', 'gray56')
#   PlotGC <- lattice::barchart(SplitSeq_GCAT, groups=TRUE, horizontal=FALSE,
#                               xlab='Nucleotide Position',
#                               scales=list(x=list(rot=90)),
#                               main='Sequence Composition', 
#                               par.settings = simpleTheme(col = my.Colors, border=my.Colors), 
#                               # For changing colors in barplot and legend and no black 
#                               # borders in the rectangles
#                               auto.key = list(adj = 1, columns = 4), 
#                               panel=function(...){
#                                 panel.barchart(...)
#                                 panel.abline(v=c(CpG_Start_Relative, CpG_End_Relative), col.line="red", lwd=2) 
#                               })
  
  PlotIntensity <- lattice::xyplot(Intensity_Normalized ~ PixelNum | MoleculeID, 
                                   data=subset(IntensityData.Aligned, MoleculeID %in% levels(IntensityData.Aligned$MoleculeID)[StartMolecule:EndMolecule]), 
                                   layout=c(1,NumPanels), 
                                   scales=list(x="free", y="free"), 
                                   type=c("l", "g"), lwd=1.25, col='darkgreen',
                                   panel = function(...) { 
                                     panel.fill(col = 'gray78') 
                                     panel.xyplot(...) 
                                   }, 
                                   main=paste('After aligning all the pixels to a reference', FragmentName.Sub))
  ##grid.arrange(PlotGC, PlotIntensity, ncol=1, heights=c(1/3,2/3))
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
                                     ylab='Mean Intensity +/- 1 standard dev')
#grid.arrange(PlotGC, PlotMeanIntensity, ncol=1, heights=c(1/2,1/2))
print(PlotMeanIntensity)

dev.off()

table(NumPixels$Pixels)
Pixel.Mode
#}


rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script looks at the CpG islands selected by Adi and plots     ##
## the aligned intensities for all the molecules aligned to each of   ##
## those fragments, along with the CpG islands.                       ##
########################################################################

########################################################################
## Initialize Header file and function library                        ##
########################################################################
OutputDataPath <- '~/Project_GC_Content/RData/'
RScriptPath <- '~/Project_GC_Content/RScripts_GC_Content/'
RDataPath <- '~/Project_GC_Content/RData/'
DataPath.Nandi <- '~/Project_GC_Content/Data/'
DataPath.CpG <- '~/Project_GC_Content/CpGIslands/'
DataPath.Steve <- '/aspen/steveg/human_nMaps/GC_content/'
Filename.Header <- paste('~/RScripts/HeaderFile_lmcg.R', sep='')
source(Filename.Header)
source(paste(RScriptPath, 'fn_Library_GC_Content.R', sep=''))
########################################################################

########################################################################
## Load the binary RData file of the Alignment Chunks, previously cr- ##
## eated                                                              ##
########################################################################
Filename.Bin <- paste(OutputDataPath, 
                      'alignmentChunks.withLength.all7134Groups.goldOnly.RData', sep='')
load(Filename.Bin)

Filename <- paste(DataPath.CpG, 'CpG_IslandsOfInterest.csv', sep='')
Fragments <- read.csv(file=Filename, header=TRUE, stringsAsFactors=F)
FragmentRow <- 5
## The first 6 rows are CpG islands provided by Adi and the last 2 are 
## Alpha-globin sites
#for(FragmentRow in 3:8){
#Chr <- 'chr19'
########################################################################
## Enter Fragment Details                                             ##
########################################################################
Chr <- Fragments[FragmentRow, 'Chr']
#FragIndex <- 65
FragIndex <- Fragments[FragmentRow, 'FragIndex']
#FragBP_Start <- 593752      ## Starting base pair coordinate of fragment
FragBP_Start <- Fragments[FragmentRow, 'FragBP_Start']
#FragBP_End <- 621533        ## Ending base pair coordinate of fragment
FragBP_End <- Fragments[FragmentRow, 'FragBP_End']
BasePairInterval <- 200     ## Length of base pair interval to estimate gcat %      
NumBP_Frag <- FragBP_End - FragBP_Start ## Length of frag in BP
NumSubFrag <- NumBP_Frag/BasePairInterval ## Number of sub fragments
#CpG_Start <- 615691         ## Starting base pair coordinate of CpG island
CpG_Start <- Fragments[FragmentRow, 'CpG_Start']
#CpG_End <- 623505           ## Ending base pair coordinate of CpG island
CpG_End <- Fragments[FragmentRow, 'CpG_End']
NumBP_CpG <- CpG_End - CpG_Start ## Length of CpG island in BP
CpG_Start_Relative <- round(max(1, (CpG_Start - FragBP_Start + 1)/BasePairInterval), 0)
                            ## max is used for those CpG regions that start before
                            ## the start of the fragment
CpG_End_Relative <- round(min(NumSubFrag, (CpG_End - FragBP_Start + 1)/BasePairInterval), 0)
                            ## min is used for those CpG regions that continue
                            ## beyond the length of this fragment

########################################################################
## DNA sequence data                                                  ##
########################################################################
SeqComp <- fn_returnSeqComp(Chr=Chr, FragIndex=FragIndex, Interval=BasePairInterval, 
                            DataPath.CpG=DataPath.CpG)
SplitSeq_GCAT <- SeqComp[['SplitSeq_GCAT']]
SeqData <- SeqComp[['SeqData']]
count(seq=SeqData, wordsize=2)
SeqData.CpGIsland <- SeqData[(CpG_Start_Relative*BasePairInterval):(CpG_End_Relative*BasePairInterval)]
count(seq=SeqData, wordsize=2)
count(seq=SeqData.CpGIsland, wordsize=2)
count(seq=SeqData.CpGIsland, wordsize=1)

SplitSeq_GCAT.DF <- as.data.frame(SplitSeq_GCAT)
colnames(SplitSeq_GCAT.DF) <- c('GG', 'CC', 'AA', 'TT')
SplitSeq_GCAT.DF$PixelNum <- index(SplitSeq_GCAT.DF)
SplitSeq_GCAT.DF$CpG <- 0
SplitSeq_GCAT.DF[CpG_Start_Relative:CpG_End_Relative, 'CpG'] <- 1
SplitSeq_GCAT.DF$CpG <- as.factor(SplitSeq_GCAT.DF$CpG)

#########################################################################
## For All molecules aligned to the fragment                           ##
#########################################################################
FragmentName <- paste(Chr, '_frag', FragIndex, '_intensities', sep='')
FragmentName.Sub <- substr(x=FragmentName, start=1, stop=(nchar(FragmentName)-12))
FragmentName.GC <- gsub(pattern='intensities', replacement='gcContent', x=FragmentName)
FragmentFilename <- paste(DataPath.CpG, FragmentName, sep='')
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
hist(NumPixels$Pixels, n=40)

#########################################################################
## Align all molecules by pixel numbers                                ##
#########################################################################
RefMoleculeID <- fn_getRefMolecule(NumPixels=NumPixels)
Reference <- subset(IntensityData, MoleculeID==RefMoleculeID)[,'Intensity']
IntensityData.Aligned <- subset(IntensityData, MoleculeID==RefMoleculeID)
# Test <- subset(IntensityData, MoleculeID=="2398367_734_2060070")[,'Intensity']
# Test.Aligned <- fn_alignPixels(Test=Test, Reference=Reference)

for(Molecule in (levels(IntensityData$MoleculeID) %w/o% RefMoleculeID)){
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
IntensityData.Reg <- merge(x=IntensityData.Aligned, y=SplitSeq_GCAT.DF, all=F)
IntensityData.Reg$GroupID <- as.factor(substr(x=IntensityData.Reg$MoleculeID, start=1, stop=7))
str(IntensityData.Aligned)
#View(IntensityData.Aligned)
str(IntensityData.Reg)
#View(IntensityData.Reg)
Filename.Out <- paste(DataPath.CpG, FragmentName.Sub, '_RegData.RData', sep='')
#save(IntensityData.Reg, file=Filename.Out)

Model1 <- lm(Intensity ~ PixelFactor, data=IntensityData.Reg)
anova(Model1)
summary(Model1)

Model2 <- lm(Intensity ~ GroupID, data=IntensityData.Reg)
anova(Model2)
summary(Model2)

Model3 <- lm(Intensity ~ GroupID + PixelFactor, data=IntensityData.Reg)
anova(Model3)
anova(Model1, Model3)
summary(Model3)

Model4 <- lm(Intensity ~ GroupID+GG+CC+TT, data=IntensityData.Reg)
summary(Model4)

Model5 <- lm(Intensity ~ GroupID + CpG, data=IntensityData.Reg)
summary(Model5)
anova(Model5)

############################ PLOTTING ################################
Molecules <- levels(IntensityData$MoleculeID)
Panels <- 5
NumPages <- ceil(length(Molecules)/Panels) ## Depends on the number of molecules aligned

Today <- Sys.Date()
Filename.pdf <- paste('~/Project_GC_Content/RScripts_GC_Content/Plots/IntensityPlots_', 
                      FragmentName.Sub, '_', Today, '.pdf', sep='')

pdf(file=Filename.pdf, onefile=TRUE, pointsize=6)
trellis.par.set(fontsize=list(text=8,points=8))

my.Colors=c('olivedrab4', 'olivedrab1', 'gray25', 'gray56')
PlotGC <- lattice::barchart(SplitSeq_GCAT, groups=TRUE, horizontal=FALSE,
                            xlab='Nucleotide Position',
                            scales=list(x=list(rot=90)),
                            main='Sequence Composition', 
                            par.settings = simpleTheme(col = my.Colors, border=my.Colors), 
                            # For changing colors in barplot and legend and no black 
                            # borders in the rectangles
                            auto.key = list(adj = 1, columns = 4), 
                            panel=function(...){
                              panel.barchart(...)
                              panel.abline(v=c(CpG_Start_Relative, CpG_End_Relative), col.line="red", lwd=2) 
                            })

Page <- 1
for(Page in 1:NumPages){
  StartMolecule <- Panels*(Page - 1) + 1
  EndMolecule <- min((StartMolecule + Panels - 1), length(Molecules))
  NumPanels <- EndMolecule - StartMolecule + 1
  
  PlotIntensity <- lattice::xyplot(Intensity_Normalized ~ PixelNum | MoleculeID, 
                                   data=subset(IntensityData, MoleculeID %in% levels(IntensityData$MoleculeID)[StartMolecule:EndMolecule]), 
                                   layout=c(1,NumPanels), 
                                   scales=list(x="free", y="free"), 
                                   type=c("l", "g"), lwd=1.25, col='darkblue',
                                   panel = function(...) { 
                                     panel.fill(col = 'gray78') 
                                     panel.xyplot(...) 
                                   }, 
                                   main=paste('Intensity plots of molecules aligned to', FragmentName.Sub))
  grid.arrange(PlotGC, PlotIntensity, ncol=1, heights=c(1/3,2/3))
  #print(PlotIntensity)
  ## This grid command prints the two charts into the pdf device
}
dev.off()

Filename.pdf <- paste('~/Project_GC_Content/RScripts_GC_Content/Plots/IntensityPlots_', 
                      FragmentName.Sub, '_Aligned_', Today, '.pdf', sep='')

pdf(file=Filename.pdf, onefile=TRUE, pointsize=6)
trellis.par.set(fontsize=list(text=8,points=8))

Page <- 1
PlotIntensity <- lattice::xyplot( Intensity ~ PixelNum, 
                                  data=IntensityData.Aligned, 
                                  groups=MoleculeID,
                                  layout=c(1,1), 
                                  scales=list(x="free", y="free"), 
                                  type=c("l", "g"), lwd=1.25,
                                  panel = function(...) { 
                                    panel.fill(col = 'gray78') 
                                    panel.superpose
                                    panel.xyplot(...) 
                                  }, 
                                  main=paste('Raw Intensities, after aligning all the pixels to a reference', FragmentName.Sub))
# my.Colors=c('olivedrab4', 'olivedrab1', 'gray25', 'gray56')
# PlotGC <- lattice::barchart(SplitSeq_GCAT, groups=TRUE, horizontal=FALSE,
#                             xlab='Nucleotide Position',
#                             scales=list(x=list(rot=90)),
#                             main='Sequence Composition', 
#                             par.settings = simpleTheme(col = my.Colors, border=my.Colors), 
#                             # For changing colors in barplot and legend and no black 
#                             # borders in the rectangles
#                             auto.key = list(adj = 1, columns = 4), 
#                             panel=function(...){
#                               panel.barchart(...)
#                               panel.abline(v=c(CpG_Start_Relative, CpG_End_Relative), col.line="red", lwd=2) 
#                             })

for(Page in 1:NumPages){
  StartMolecule <- Panels*(Page - 1) + 1
  EndMolecule <- min((StartMolecule + Panels - 1), length(Molecules))
  NumPanels <- EndMolecule - StartMolecule + 1
                              
  PlotIntensity.Norm <- lattice::xyplot(Intensity_Normalized ~ PixelNum | MoleculeID, 
                                        data=subset(IntensityData.Aligned, MoleculeID %in% levels(IntensityData.Aligned$MoleculeID)[StartMolecule:EndMolecule]), 
                                        layout=c(1,NumPanels), 
                                        scales=list(x="free", y="free"), 
                                        type=c("l", "g"), lwd=1.25, col='darkgreen',
                                        panel = function(...) { 
                                          panel.fill(col = 'gray78') 
                                          panel.xyplot(...) 
                                        }, 
                                        main=paste('After aligning all the pixels to a reference', FragmentName.Sub))
  grid.arrange(PlotGC, PlotIntensity.Norm, ncol=1, heights=c(1/3,2/3))
  #print(PlotIntensity)
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
grid.arrange(PlotGC, PlotIntensity, PlotMeanIntensity, ncol=1, 
             heights=c(1/3,1/3,1/3))
#print(PlotMeanIntensity)

dev.off()

#}

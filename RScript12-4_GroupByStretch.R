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
##########################################################################
## Up to this line its the same as RScript12_CpG_Analysis.R             ##
##########################################################################
MedianLength <- median(NumPixels$Pixels)
Length_Min <- 0.9*MedianLength
Length_Max <- 1.1*MedianLength
NumPixels$Outlier <- !(NumPixels$Pixels <= Length_Max & NumPixels$Pixels >= Length_Min)
NumPixels.NoOutlier <- subset(NumPixels, !(Outlier))
RangeofPixelLengths <- range(NumPixels.NoOutlier[,'Pixels'])
Length_Range <- RangeofPixelLengths[2] - RangeofPixelLengths[1]
NumPixels.NoOutlier$Pixels[which.min(NumPixels.NoOutlier$Pixels)]
NumGroups <- ifelse(Length_Range <=4, 1, 
                    ifelse(Length_Range <= 8, 2, 
                           ifelse(Length_Range <= 14, 3, 
                                  ifelse(Length_Range <= 24, 4, 
                                         5))))
VectorOfPixelLengths <- as.numeric(names(table(NumPixels.NoOutlier$Pixels)))

## Use fn_subseq function to split this vector into NumGroups groups
## Once divided up into groups, use fn_Mode to find the mode length in each groups
## Use the NMap with the mode length as the reference NMap and align the rest 
## of the NMaps in that group to this reference

PixelGroups <- fn_subseq(x=VectorOfPixelLengths, n=NumGroups)
NumPixels.NoOutlier$PixelGroup <- 0
for(PixelGroup in 1:length(PixelGroups)){
  Pixels.Selected <- PixelGroups[[PixelGroup]]
  NumPixels.NoOutlier[NumPixels.NoOutlier$Pixels %in% c(Pixels.Selected),'PixelGroup'] <- PixelGroup
}

sum(NumPixels$Outlier)
hist(NumPixels$Pixels[NumPixels$Outlier==FALSE], n=40)
hist(NumPixels$Pixels, n=40)


#########################################################################
## Align all molecules by pixel numbers                                ##
#########################################################################
PixelGroup <- 3
for(PixelGroup in 1:length(PixelGroups)){
  
  RefMoleculeID <- fn_getRefMolecule(NumPixels=NumPixels.NoOutlier[NumPixels.NoOutlier$PixelGroup==PixelGroup,])
  Reference <- subset(IntensityData, MoleculeID==RefMoleculeID)[,'Intensity']
  IntensityData.Aligned <- subset(IntensityData, MoleculeID==RefMoleculeID)
  # Test <- subset(IntensityData, MoleculeID=="2398367_734_2060070")[,'Intensity']
  # Test.Aligned <- fn_alignPixels(Test=Test, Reference=Reference)
  
  MoleculeIDs_toAlign <- NumPixels.NoOutlier[NumPixels.NoOutlier$PixelGroup==PixelGroup, 'MoleculeID'] %w/o% RefMoleculeID
  for(Molecule in MoleculeIDs_toAlign){
    print(Molecule)
    Test <- subset(IntensityData, MoleculeID==Molecule)[,'Intensity']
    print(length(Test))
    Test.Aligned <- fn_alignPixels(Test=Test, Reference=Reference, Simulation_N=10000, 
                                   verbose=TRUE)
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
  
  IntensityData.Norm.Aligned.Wide <- reshape(IntensityData.Aligned[,-1], 
                                             timevar='MoleculeID', 
                                             idvar='PixelNum', 
                                             direction='wide')
  
  NumCol <- ncol(IntensityData.Aligned.Wide)
  IntensityData.Aligned.Wide$Intensity_Mean <- rowMeans(IntensityData.Aligned.Wide[,2:NumCol])
  IntensityData.Aligned.Wide$Intensity_SD <- apply(X=IntensityData.Aligned.Wide[,2:NumCol], 
                                                   MARGIN=1, FUN=sd)
  IntensityData.Aligned.Wide$Intensity_Up <- IntensityData.Aligned.Wide$Intensity_Mean + 1*IntensityData.Aligned.Wide$Intensity_SD
  IntensityData.Aligned.Wide$Intensity_Dn <- IntensityData.Aligned.Wide$Intensity_Mean - 1*IntensityData.Aligned.Wide$Intensity_SD
  
  IntensityData.Norm.Aligned.Wide$Intensity_Mean <- rowMeans(IntensityData.Norm.Aligned.Wide[,2:NumCol])
  IntensityData.Norm.Aligned.Wide$Intensity_SD <- apply(X=IntensityData.Norm.Aligned.Wide[,2:NumCol], 
                                                        MARGIN=1, FUN=sd)
  IntensityData.Norm.Aligned.Wide$Intensity_Up <- IntensityData.Norm.Aligned.Wide$Intensity_Mean + 1*IntensityData.Norm.Aligned.Wide$Intensity_SD
  IntensityData.Norm.Aligned.Wide$Intensity_Dn <- IntensityData.Norm.Aligned.Wide$Intensity_Mean - 1*IntensityData.Norm.Aligned.Wide$Intensity_SD
  
  ########################################################################
  ## DNA sequence data                                                  ##
  ########################################################################
  BasePairInterval <- ceil(NumBP_Frag/nrow(IntensityData.Aligned.Wide))
  
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
  ########################################################################
  
  ############################### ANOVA ################################
  IntensityData.Aligned$PixelPosition <- as.factor(IntensityData.Aligned$PixelNum)
  IntensityData.Reg <- merge(x=IntensityData.Aligned, y=SplitSeq_GCAT.DF, all=F)
  IntensityData.Reg$GroupID <- as.factor(substr(x=IntensityData.Reg$MoleculeID, start=1, stop=7))
  str(IntensityData.Aligned)
  
  Model <- lm(Intensity_Normalized ~ PixelPosition, data=IntensityData.Reg)
  summary(Model)
  PValues <- summary(Model)$coefficients[,4]
  Model_pValue <- round(anova(Model)$'Pr(>F)'[1], 4)
  ########################################################################
  
  ################################ PLOTS #################################
  Today <- Sys.Date()
  Filename.pdf <- paste('~/Project_GC_Content/RScripts_GC_Content/Plots/IntensityPlots_', 
                        FragmentName.Sub, '_PixelGroup_', PixelGroup, '_', Today, '.pdf', sep='')
  pdf(file=Filename.pdf, onefile=TRUE, pointsize=6)
  trellis.par.set(fontsize=list(text=8,points=8))
  
  my.Colors=c('olivedrab4', 'olivedrab1', 'gray25', 'gray56')
  PlotGC <- lattice::barchart(SplitSeq_GCAT, groups=TRUE, horizontal=FALSE,
                              xlab='Nucleotide Position',
                              scales=list(x=list(rot=90)),
                              main=paste(FragmentName.Sub, 'PixelGroup:', PixelGroup), 
                              par.settings = simpleTheme(col = my.Colors, border=my.Colors), 
                              # For changing colors in barplot and legend and no black 
                              # borders in the rectangles
                              auto.key = list(adj = 1, columns = 4), 
                              panel=function(...){
                                panel.barchart(...)
                                panel.abline(v=c(CpG_Start_Relative, CpG_End_Relative), col.line="red", lwd=2) 
                              })

  PlotIntensity <- lattice::xyplot(Intensity_Normalized ~ PixelNum | MoleculeID, 
                                   #data=subset(IntensityData.Aligned, MoleculeID %in% levels(IntensityData.Aligned$MoleculeID)[StartMolecule:EndMolecule]), 
                                   data=subset(IntensityData.Aligned, MoleculeID %in% levels(IntensityData.Aligned$MoleculeID)[1:8]),
                                   layout=c(1,5), 
                                   scales=list(x="free", y="free"), 
                                   type=c("l", "g"), lwd=1.25, col='darkgreen',
                                   panel = function(...) { 
                                     panel.fill(col = 'gray78') 
                                     panel.xyplot(...) 
                                   }, 
                                   main=paste('After aligning all the pixels to a reference', FragmentName.Sub))
  
  PlotMeanIntensity <- lattice::xyplot(Intensity_Mean + Intensity_Up + Intensity_Dn ~ PixelNum, 
                                       data=subset(IntensityData.Norm.Aligned.Wide), 
                                       layout=c(1,1), 
                                       scales=list(x="free", y="free"), 
                                       type=c("l", "g"), col=c('purple4', 'purple', 'purple'),
                                       lty=c(1, 2, 2), lwd=c(1.5, 1.5, 1.5),
                                       panel = function(...) { 
                                         panel.fill(col = 'gray78') 
                                         panel.xyplot(...) 
                                       }, 
                                       main=paste('Mean Relative Intensity plot', 'p-value', 
                                                  Model_pValue, ', Num of NMaps', length(PixelGroups[[PixelGroup]])), 
                                       ylab='Mean Intensity +/- 1 standard dev', 
                                       xlab='Pixel Position')
  grid.arrange(PlotGC, PlotMeanIntensity, ncol=1, heights=c(1/2,1/2))
  grid.arrange(PlotGC, PlotIntensity, ncol=1, heights=c(1/3,2/3))
  
  textplot(capture.output(summary(Model)))
  dev.off()
}

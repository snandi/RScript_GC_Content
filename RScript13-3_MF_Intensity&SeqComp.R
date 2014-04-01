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
FragIndex <- 33
BasePairInterval <- 209     ## Length of base pair interval to estimate gcat %      
## Starting base pair coordinate of fragment
FragBP_Start <- bp.loc[which(bp.loc$alignedFragIndex == FragIndex), 'refMapCoordStart'] 
## Ending base pair coordinate of fragment
FragBP_End <- bp.loc[which(bp.loc$alignedFragIndex == FragIndex), 'refMapCoordEnd']

NumBP_Frag <- FragBP_End - FragBP_Start ## Length of frag in BP
NumSubFrag <- round(NumBP_Frag/BasePairInterval, 0) ## Number of sub fragments
                                          ## Also, the number of pixels

#for(FragIndex in 1:39){

########################################################################
## DNA sequence data                                                  ##
########################################################################
# Interval <- NumSubFrag
# Filename <- paste(DataPath.MFlorum, 'mesoplasma_florum_1.fasta', sep='')
# MF_Seq <- read.fasta(file = Filename)[[1]]
# SeqTable <- count(seq=MF_Seq, wordsize=1)                                                                                                           
# SeqGC <- GC(MF_Seq)
# BasePosition <- round(seq(from=FragBP_Start, to=FragBP_End, by=BasePairInterval), 0)
# SplitSeq <- fn_subseq(x=MF_Seq[FragBP_Start:FragBP_End], n=Interval, force.number.of.groups=TRUE)
# #SplitSeq <- fn_subseq(x=MF_Seq, n=ceil(max(BasePosition)/Interval), force.number.of.groups=TRUE)
# SplitSeq_ACGT <- do.call(what=rbind, lapply(X=SplitSeq, FUN=count, wordsize=1))
# rownames(SplitSeq_ACGT) <- BasePosition[1:nrow(SplitSeq_ACGT)]
# colnames(SplitSeq_ACGT) <- c('A', 'C', 'G', 'T')
# SplitSeq_ACGT <- SplitSeq_ACGT/rowSums(SplitSeq_ACGT)
# SplitSeq_GCAT <- SplitSeq_ACGT[,c('G', 'C', 'A', 'T')]
# head(SplitSeq_GCAT)
# 
# SeqComp <- fn_returnSeqComp(Chr=Chr, FragIndex=FragIndex, Interval=BasePairInterval, 
#                             DataPath=DataPath.MFlorum, numPixels=NumSubFrag, 
#                             paste(DataPath.MFlorum, 'mesoplasma_florum_1.fasta', sep=''))
# SplitSeq_GCAT1 <- SeqComp[['SplitSeq_GCAT']]
# SeqData <- MF_Seq
# SeqData1 <- SeqComp[['SeqData']]
# 
# count(seq=MF_Seq, wordsize=2)

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
                                            FragmentData=FragmentData, 
                                            Truncate=FALSE,
                                            TruncateLength=0))

NumPixels <- aggregate(IntensityData$PixelNum, by=list(IntensityData$MoleculeID), FUN=max)
colnames(NumPixels) <- c('MoleculeID', 'Pixels')
table(NumPixels$Pixels)

#########################################################################
## Output Table                                                        ##
#########################################################################
#1 
FragmentName.Sub
#2 Stretch
Stretch <- NULL
#3 overall pValue
Model_pValue <- NULL
#4 RSS (summary(Model1.Norm)$sigma^2)
Model_RSS <- NULL
#5 Model R^2 (summary(Model1.Norm)$r.squared)
Model_RSquared <- NULL
#6 Num of significant coeff pValues
Sig_pValues <- NULL
#7 Num of significant Benjamini Hockberg adjusted coeff pValues
Sig_pValues.BH <- NULL
#8 Number of nMaps (length(Molecules))
nMaps <- NULL

Index <- 0
#########################################################################
## Align all molecules by pixel numbers                                ##
#########################################################################
Pixel.Mode <- Mode(NumPixels$Pixels)
FragmentLength <- Pixel.Mode
FragmentLengths <- c(44:53)

for(FragmentLength in FragmentLengths){
  Index <- Index + 1
  Stretch[Index] <- FragmentLength
  IntensityData.Aligned <- subset(IntensityData, MoleculeID %in% 
                                    subset(NumPixels, Pixels==FragmentLength)[,'MoleculeID'])
  IntensityData.Aligned.Wide <- reshape(IntensityData.Aligned[,colnames(IntensityData.Aligned) %w/o% 'Intensity_Normalized'],
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
  IntensityData.Aligned$GroupID <- substr(IntensityData.Aligned$MoleculeID, 
                                          start=1, stop=7)
  IntensityData.Aligned <- within(data=IntensityData.Aligned, {
                                  PixelFactor <- as.factor(PixelNum)
                                  GroupID <- as.factor(GroupID)
                                  })
  #str(IntensityData.Aligned)
  Model1 <- lm(Intensity ~ PixelFactor, data=IntensityData.Aligned)
  pValues.1 <- coef(summary(Model1))[,4]
  #anova(Model1)
  #summary(Model1)
  pValues1.BH <- p.adjust(pValues.1, method='BH')
  
  Model1.Norm <- lm(Intensity_Normalized ~ PixelFactor, data=IntensityData.Aligned)
  pValues1.Norm <- coef(summary(Model1.Norm))[,4]
  
  #anova(Model1.Norm)
  #summary(Model1.Norm)
  pValues1.Norm.BH <- p.adjust(pValues1.Norm, method='BH')
  

  Model_pValue[Index] <- anova(Model1.Norm)$'Pr(>F)'[1]
  Model_RSS[Index] <- round(summary(Model1.Norm)$sigma^2, 4)
  Model_RSquared[Index] <- round(summary(Model1.Norm)$r.squared, 4)
  Sig_pValues[Index] <- sum(pValues1.Norm < 0.05)
  Sig_pValues.BH[Index] <- sum(pValues1.Norm.BH < 0.05)
  
  #   Model2 <- lm(Intensity ~ GroupID, data=IntensityData.Aligned)
  #   #anova(Model2)
  #   #summary(Model2)
  #   
  #   Model3 <- lm(Intensity ~ GroupID + PixelFactor, data=IntensityData.Aligned)
  #   #anova(Model3)
  #   #anova(Model2, Model3)
  #   #summary(Model3)
  #   pValues.3 <- coef(summary(Model3))[,4]  
  #   pValues.3.BH <- p.adjust(pValues.3, method='BH')
  # 
  #   Model3.Norm <- lm(Intensity_Normalized ~ GroupID + PixelFactor, data=IntensityData.Aligned)
  #   summary(Model3.Norm)
  #   anova(Model3.Norm)
  ############################ PLOTTING ################################
  IntensityData.Aligned$MoleculeID <- as.vector(IntensityData.Aligned$MoleculeID)
  Molecules <- unique(IntensityData.Aligned$MoleculeID)
  
  nMaps[Index] <- length(Molecules)
  
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
  
  SeqComp <- fn_returnSeqComp(Chr=Chr, FragIndex=FragIndex, Interval=BasePairInterval, 
                              DataPath=DataPath.MFlorum, numPixels=FragmentLength, 
                              paste(DataPath.MFlorum, 'mesoplasma_florum_1.fasta', sep=''), 
                              FragBP_Start=FragBP_Start, FragBP_End=FragBP_End)
  SplitSeq_GCAT <- SeqComp[['SplitSeq_GCAT']]
  
  my.Colors=c('olivedrab4', 'olivedrab1', 'gray25', 'gray56')
  PlotGC <- lattice::barchart(SplitSeq_GCAT, groups=TRUE, horizontal=FALSE,
                              xlab='Nucleotide Position',
                              # scales=list(x=list(rot=90)),
                              scales=list(x=list(draw=FALSE)),
                              main='Sequence Composition', 
                              par.settings = simpleTheme(col = my.Colors, border=my.Colors), 
                              # For changing colors in barplot and legend and no black 
                              # borders in the rectangles
                              auto.key = list(adj = 1, columns = 4)
  )
  
  for(Page in 1:NumPages){
    StartMolecule <- Panels*(Page - 1) + 1
    EndMolecule <- min((StartMolecule + Panels - 1), length(Molecules))
    NumPanels <- EndMolecule - StartMolecule + 1
    
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
    grid.arrange(PlotGC, PlotIntensity, ncol=1, heights=c(1/4,3/4))
    #  print(PlotIntensity)
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
                                       main=paste('Mean Intensity of', length(Molecules), 'nMaps, ', FragmentLength, 'Pixels'), 
                                       ylab='Mean Intensity +/- 1 standard dev')
  grid.arrange(PlotGC, PlotMeanIntensity, ncol=1, heights=c(1/2,1/2))
  #print(PlotMeanIntensity)
  
  dev.off()
  
  #table(NumPixels$Pixels)
  print(FragmentLength)
} #for(FragmentLength in FragmentLengths)
#} #for(FragIndex in 1:39)

Fragment <- rep(x=FragmentName.Sub, times=Index)

Output <- as.data.frame(cbind(Stretch, Model_pValue, Model_RSS, Model_RSquared, 
                              Sig_pValues, Sig_pValues.BH, 
                              nMaps), stringsAsFactors=FALSE)
Output <- within(data=Output, {
  Stretch <- as.numeric(Stretch)
  Model_pValue <- as.numeric(Model_pValue)
  Model_RSS <- as.numeric(Model_RSS)
  Model_RSquared <- as.numeric(Model_RSquared)
  Sig_pValues <- round(as.numeric(Sig_pValues), 0)
  Sig_pValues.BH <- round(as.numeric(Sig_pValues.BH), 0)
  nMaps <- round(as.numeric(nMaps), 0)
})

Output.xtable <- xtable(Output)
digits(Output.xtable) <- c(0,0, 4, 4, 4, 0, 0, 0)
caption(Output.xtable) <- FragmentName.Sub
print(Output.xtable, include.rownames=FALSE)

Cor.Table <- round(cor(Output), 4)
print(xtable(Cor.Table))

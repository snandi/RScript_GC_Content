rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script conducts both parametric and non-parametric anova, as  ##
## well as transformations of the data, trying to establish repro-    ##
## -ducibility of the intensity profiles across nMaps.                ##
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
bp.loc$BasePairLength <- bp.loc$refMapCoordEnd - bp.loc$refMapCoordStart
bp.loc$PixelLength <- round(bp.loc$BasePairLength/209, 0)
write.table(bp.loc, file=paste(RDataPath, 'mflorum.bploc.txt', sep=''), 
            row.names=FALSE, sep='\t')
table(bp.loc$PixelLength)
########################################################################
## Enter Fragment Details                                             ##
########################################################################
Chr <- 'chr1'
FragIndex <- 2
BasePairInterval <- 209     ## Length of base pair interval to estimate gcat %      
## Starting base pair coordinate of fragment
FragBP_Start <- bp.loc[which(bp.loc$alignedFragIndex == FragIndex), 'refMapCoordStart'] 
## Ending base pair coordinate of fragment
FragBP_End <- bp.loc[which(bp.loc$alignedFragIndex == FragIndex), 'refMapCoordEnd']

NumBP_Frag <- FragBP_End - FragBP_Start ## Length of frag in BP
NumSubFrag <- round(NumBP_Frag/BasePairInterval, 0) ## Number of sub fragments
## Also, the number of pixels

#for(FragIndex in 1:38){

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

Truncate <- TRUE
TruncateLength=10

IntensityData <- do.call(what=rbind, lapply(X=Molecules, FUN=fn_returnMoleculeIntensity, 
                                            FragmentData=FragmentData, 
                                            Truncate=Truncate,
                                            TruncateLength=TruncateLength))

NumPixels <- aggregate((IntensityData$PixelNum - TruncateLength), 
                       by=list(IntensityData$MoleculeID), FUN=max)
colnames(NumPixels) <- c('MoleculeID', 'Pixels')
NumPixels.Table <- table(NumPixels$Pixels)

#########################################################################
## Outlier Molecule IDs                                                ##
#########################################################################
# Outlier_MoleculeIDs <- c('2428070_734_2090032', '2428134_734_2090169',
#                          '2427995_734_2090121', '2428086_734_2090170',
#                          '2428000_734_2090104', '2428056_734_2090121',
#                          '2428065_734_2090074', '2428138_734_2090068', 
#                          '2428082_734_2090091', '2428155_734_2090009',
#                          '2428009_734_2090050')
# IntensityData <- subset(IntensityData, !(MoleculeID %in% Outlier_MoleculeIDs))

#########################################################################
## Output Table                                                        ##
#########################################################################
#1 
FragmentName.Sub
#2 Stretch
Stretch <- NULL
#3 overall pValue
pValue.Raw <- NULL
#3a overall pValue after log
pValue.Log <- NULL
#3b overall pValue after Kruskal Wallis rank sum test
pValue.KW <- NULL
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
FragmentLengths <- FragmentLength
if(length(NumPixels.Table[NumPixels.Table > 20]) > 0){
  FragmentLengths <- as.numeric(names(NumPixels.Table[NumPixels.Table >= 20]))
}

#FragmentLengths <- c(76:84) ## For Fragment 33
#FragmentLengths <- c(72:85) ## For Fragment 11
#FragmentLength <- 68

if(length(FragmentLengths) > 0){
  for(FragmentLength in FragmentLengths){
    Index <- Index + 1
    Stretch[Index] <- FragmentLength
    IntensityData.Aligned <- subset(IntensityData, MoleculeID %in% 
                                      subset(NumPixels, Pixels==FragmentLength)[,'MoleculeID'])
    ### Take out Outliers (more than +/- 20% different) ###
    Outlier_MoleculeIDs <- c(as.vector(unique(IntensityData.Aligned$MoleculeID[IntensityData.Aligned$Intensity_Normalized >= 1.2])), 
                             as.vector(unique(IntensityData.Aligned$MoleculeID[IntensityData.Aligned$Intensity_Normalized <= 0.8])))
    IntensityData.Aligned <- subset(IntensityData.Aligned, !(MoleculeID %in% Outlier_MoleculeIDs))
    #######################################################
    
    IntensityData.Aligned.Wide <- reshape(IntensityData.Aligned[,colnames(IntensityData.Aligned)%w/o%'Intensity_Normalized'], 
                                          timevar='MoleculeID', 
                                          idvar='PixelNum', 
                                          direction='wide')
    
    #   if(FragmentLength == 82){
    DataFilename <- paste(RDataPath, FragmentName.Sub, '_Length_', FragmentLength, 
                          '.RData', sep='')
    save(IntensityData.Aligned, file=DataFilename)
    #   }
    
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
    
    Model1.Norm <- lm(Intensity_Normalized ~ PixelFactor, data=IntensityData.Aligned)
    pValues1.Norm <- coef(summary(Model1.Norm))[,4]
    #anova(Model1.Norm)
    #summary(Model1.Norm)
    
    Model1.Norm.Log <- lm(log(Intensity_Normalized) ~ PixelFactor, data=IntensityData.Aligned)
    pValues1.Norm.Log <- coef(summary(Model1.Norm.Log))[,4]
    anova(Model1.Norm.Log)
    #summary(Model1.Norm.Log)
    
    pValue.Raw[Index] <- round(anova(Model1.Norm)$'Pr(>F)'[1], 4)
    pValue.Log[Index] <- round(anova(Model1.Norm.Log)$'Pr(>F)'[1], 4)
    # Model_RSS[Index] <- round(summary(Model1.Norm)$sigma^2, 4)
    # Model_RSquared[Index] <- round(summary(Model1.Norm)$r.squared, 4)
    
    Kruskal_Wallis <- kruskal.test(log(Intensity_Normalized) ~ PixelFactor, 
                                   data=IntensityData.Aligned)
    pValue.KW[Index] <- round(Kruskal_Wallis$p.value, 4)
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
                          FragmentName.Sub, '_Pixels_Trunc', nrow(IntensityData.Aligned.Wide), 
                          '_', Today, '.pdf', sep='')
    
    pdf(file=Filename.pdf, onefile=TRUE, pointsize=6)
    trellis.par.set(fontsize=list(text=8,points=8))
    
    PlotHist <- lattice::barchart(x=table(NumPixels$Pixels), col='gray50',
                                  xlab='Number of nMaps', ylab='Pixel length of nMaps',
                                  main=paste('Histogram of nMap Pixels aligned to', FragmentName.Sub))
    print(PlotHist)
    
    SeqComp <- fn_returnSeqComp(Chr=Chr, FragIndex=FragIndex, Interval=BasePairInterval, 
                                DataPath=DataPath.MFlorum, 
                                numPixels=(FragmentLength+2*TruncateLength), 
                                paste(DataPath.MFlorum, 'mesoplasma_florum_1.fasta', sep=''), 
                                FragBP_Start=FragBP_Start, FragBP_End=FragBP_End)
    SplitSeq_GCAT <- SeqComp[['SplitSeq_GCAT']]
    rownames(SplitSeq_GCAT) <- index(SplitSeq_GCAT)
    #SplitSeq_GCAT$PixelNum <- as.vector(rownames(SplitSeq_GCAT))
    
    my.Colors=c('olivedrab4', 'olivedrab1', 'gray25', 'gray56')
    KeepFrom <- 1 + TruncateLength
    KeepTo <- nrow(SplitSeq_GCAT) - TruncateLength
    
    PlotGC <- lattice::barchart(SplitSeq_GCAT, groups=TRUE, horizontal=FALSE,
                                xlab='Nucleotide Position',
                                scales=list(x=list(rot=90)),
                                # scales=list(x=list(draw=FALSE)),
                                main=paste('Sequence Composition of Fragment', FragIndex),  
                                par.settings = simpleTheme(col = my.Colors, border=my.Colors), 
                                # For changing colors in barplot and legend and no black 
                                # borders in the rectangles
                                auto.key = list(adj = 1, columns = 4),
                                panel=function(...){
                                  panel.barchart(...)
                                  panel.abline(v=c(KeepFrom, KeepTo), col.line="red", lwd=2) 
                                }
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
    print(c(FragmentName.Sub, FragmentLength))
  } #for(FragmentLength in FragmentLengths)
} else{ # if condition of fragmentLength
  print(FragmentName.Sub)
}
#} #for(FragIndex in 0:38)

Fragment <- rep(x=FragmentName.Sub, times=Index)

Output <- as.data.frame(cbind(Stretch, pValue.Raw, pValue.Log, pValue.KW, 
                              nMaps), stringsAsFactors=FALSE)
str(Output)

# Output <- within(data=Output, {
#   Stretch <- as.numeric(Stretch)
#   Model_pValue <- as.numeric(Model_pValue)
#   Model_RSS <- as.numeric(Model_RSS)
#   Model_RSquared <- as.numeric(Model_RSquared)
#   Sig_pValues <- round(as.numeric(Sig_pValues), 0)
#   Sig_pValues.BH <- round(as.numeric(Sig_pValues.BH), 0)
#   nMaps <- round(as.numeric(nMaps), 0)
# })

print(Output)

Output.xtable <- xtable(Output)
digits(Output.xtable) <- c(0, 0, 4, 4, 4, 0) ## Number of columns + 1
caption(Output.xtable) <- FragmentName.Sub
print(Output.xtable, include.rownames=FALSE)

xtable(t(count(seq=SeqComp[['SeqData']], wordsize=1)))
xtable(t(count(seq=SeqComp[['SeqData']], wordsize=2)))


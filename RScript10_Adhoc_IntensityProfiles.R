rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script is for plotting intensity profiles of multiple molecu- ##
## -les aligned to the same fragment location on the genome
########################################################################

########################################################################
## Initialize Header file and function library                        ##
########################################################################
FilePath <- '~/Project_GC_Content/RScripts_GC_Content/'
DataPath <- '~/Project_GC_Content/Adhoc_2013-09-25/'
Filename.Header <- paste('~/RScripts/HeaderFile_lmcg.R', sep='')
source(Filename.Header)
source(paste(FilePath, 'fn_Library_GC_Content.R', sep=''))
########################################################################

AlignedChr <- 11
AlignedFragIndex <- 4854

GCFilename <- paste(DataPath, 'GC_Content.txt', sep='')
GCFile <- scan(file=GCFilename, what=double())
str(GCFile)

IndexFilename <- paste(DataPath, 'GroupMoleculeInfo.txt', sep='')
IndexFile <- read.table(file=IndexFilename, header=TRUE, sep='\t', 
                        skip=0, stringsAsFactors=FALSE)
IndexFile <- within(IndexFile,{
  Type <- factor(Type)
  MoleculeName <- paste('Gr', GroupID, '_Mol', MoleculeID, sep='')
  MoleculeName <- factor(MoleculeName)
})
str(IndexFile)

fn_loadMoleculeFile <- function(Row=9, DataPath, IndexFile){
  MoleculeFilename <- paste(DataPath, 'Gr', IndexFile[Row,'GroupID'], '_Mol', 
                            IndexFile[Row,'MoleculeID'], '.txt', sep='')
  Intensity <- scan(file=MoleculeFilename, what=integer())
  MoleculeName <- IndexFile[Row,'MoleculeName']
  MoleculeType <- IndexFile[Row,'Type']
  if(MoleculeType == 'Reverse'){
    Intensity <- rev(Intensity)
  }
  IntensityData <- as.data.frame(Intensity)
  IntensityData$MoleculeName <- MoleculeName
  IntensityData$Intensity_Normalized <- IntensityData$Intensity/median(IntensityData$Intensity)
  IntensityData$PixelNum <- index(IntensityData)
  IntensityData$MoleculeType <- MoleculeType
  #return(list(MoleculeName=MoleculeName, MoleculeFile=MoleculeFile))
  return(IntensityData)
}

#########################################################################
## For One molecule                                                    ##
#########################################################################
# Row <- 1
# IntensityData <- rbind(fn_loadMoleculeFile(Row=1, DataPath=DataPath), 
#                        fn_loadMoleculeFile(Row=2, DataPath=DataPath))
# IntensityData$MoleculeName <- as.factor(IntensityData$MoleculeName)
# str(IntensityData)

Rows <- c(1:nrow(IndexFile))
IntensityData <- do.call(what=rbind, lapply(X=Rows, FUN=fn_loadMoleculeFile, 
                                            DataPath=DataPath, IndexFile=IndexFile))
str(IntensityData)
#IntensityData$MoleculeName <- as.factor(IntensityData$MoleculeName)

GCData <- as.data.frame(GCFile)
names(GCData) <- 'Intensity'
GCData$PixelNum <- index(GCData)
GCData$MoleculeName <- 'GC'
GCData$Intensity_Normalized <- GCData$Intensity
GCData$MoleculeType <- 'GC'

IntensityData <- rbind(IntensityData, GCData)

NumMolecules.Blue <- try(length(unique(subset(IntensityData, MoleculeType=='Blue')[,'MoleculeName'])) + 1) 
NumMolecules.Gold <- try(length(unique(subset(IntensityData, MoleculeType=='Gold')[,'MoleculeName'])) + 1) 
NumMolecules.Reverse <- try(length(unique(subset(IntensityData, MoleculeType=='Reverse')[,'MoleculeName'])) + 1)
## The +1 is to add the GC content plot on top of all the molecules

Today <- Sys.Date()
Filename.Gold <- paste('~/Project_GC_Content/Adhoc_2013-09-25/IntensityPlots_', Today, 
                      '_Gold.jpg', sep='')
jpeg(file = Filename.Gold)
# Filename.pdf <- paste('~/Project_GC_Content/Adhoc_2013-09-25/IntensityPlots_', Today, 
#                       '.pdf', sep='')
# pdf(file=Filename.pdf, onefile=TRUE, pointsize=6)
lattice::xyplot(Intensity_Normalized ~ PixelNum | MoleculeName, 
                data=subset(IntensityData, MoleculeType=='Gold' | MoleculeType=='GC'), 
                layout=c(1,NumMolecules.Gold), 
                scales=list(x="free", y="free"), 
                type=c("l", "g"), lwd=2, col='gold',
                panel = function(...) { 
                  panel.fill(col = 'gray7') 
                  panel.xyplot(...) 
                }, 
                main='Gold Straight')
dev.off()

Filename.Rev <- paste('~/Project_GC_Content/Adhoc_2013-09-25/IntensityPlots_', Today, 
                      '_Reverse.jpg', sep='')
jpeg(file = Filename.Rev)
lattice::xyplot(Intensity_Normalized ~ PixelNum | MoleculeName, 
                data=subset(IntensityData, MoleculeType=='Reverse' | MoleculeType=='GC'), 
                layout=c(1,NumMolecules.Reverse), 
                scales=list(x="free", y="free"), 
                type=c("l", "g"), lwd=2, col='gold',
                panel = function(...) { 
                  panel.fill(col = 'gray50') 
                  panel.xyplot(...) 
                },
                main='Gold Reversed')
dev.off()

Filename.Blue <- paste('~/Project_GC_Content/Adhoc_2013-09-25/IntensityPlots_', Today, 
                       '_Blue.jpg', sep='')
jpeg(file = Filename.Blue)
lattice::xyplot(Intensity_Normalized ~ PixelNum | MoleculeName, 
                data=subset(IntensityData, MoleculeType=='Blue' | MoleculeType=='GC'), 
                layout=c(1,NumMolecules.Blue), 
                scales=list(x="free", y="free"), 
                type=c("l", "g"), lwd=2, col='light blue',
                panel = function(...) { 
                  panel.fill(col = "gray7") 
                  panel.xyplot(...) 
                },
                main='Blue Straight')
dev.off()



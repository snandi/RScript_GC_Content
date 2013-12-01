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
DataPath <- '/exports/aspen/steveg/human_nMaps/GC_content/'
AlignmentFilename <- 'alignmentChunks.withLength.all7134Groups.goldOnly'
########################################################################

setwd(DataPath)
Filename <- paste(DataPath, AlignmentFilename, sep='')

Colnames <- c('refChr', 'refStartIndex', 'refEndIndex', 'molID', 'molStartIndex', 
                       'molEndIndex', 'refStartCoord', 'refEndCoord', 'molStartCoord', 
                       'molEndCoord', 'orient', 'lengthRatio')

Colnames.Steve <- c('refID', 'refStartIndex', 'refEndIndex', 'opID', 
'opStartIndex', 'opEndIndex', 'refStartCoord', 'refEndCoord', 'opStartCoord', 'opEndCoord', 
'orient', 'lengthRatio')

########################################################################
## Run this part only for a new alignment, else use the saved binary  ##
## RData file already produced by a previous run                      ##
########################################################################
# ## Read in Alignment Chunk file and save the binary
# AlChunk <- fn_readAlignmentChunks(Filename=Filename, Header=FALSE, Colnames=Colnames, 
#                                   Colnames.Steve=Colnames.Steve)
# 
# str(AlChunk)
# Filename.Output <- paste(OutputDataPath, 'alignmentChunks.withLength.all7134Groups.goldOnly.RData', sep='')
# save(AlChunk, file=Filename.Output)

########################################################################
## Load the binary RData file of the Alignment Chunks, previously cr- ##
## eated                                                              ##
########################################################################
Filename.Bin <- paste(OutputDataPath, 
                      'alignmentChunks.withLength.all7134Groups.goldOnly.RData', sep='')
load(Filename.Bin)

## How many unique molecules?
length(unique(AlChunk$molID))

## How many fragments in each chromosome?
# Chr_Frags <- unique(AlChunk[,c('refChr', 'refStartIndex')])
# head(Chr_Frags)
# table(Chr_Frags$refChr)

########################################################################
## How many molecules are aligned to any given chr & frag index?      ##
########################################################################
Chromosomes <- c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 
                 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY')

for(Chr in Chromosomes){
  TableName <- paste(Chr, '_Table', sep='')
  Filename.Out <- paste(OutputDataPath, Chr, '_Table.RData', sep='')

  ChrTable <- fn_numMolAlignedperLoc(AlChunk=subset(AlChunk, refChr==Chr))
  assign(x=paste(Chr, '_Table', sep=''), value=ChrTable)
  save(list=paste(Chr, '_Table', sep=''), file=Filename.Out)
  rm(ChrTable)
}
########################################################################

########################################################################
## This for loop does the following: For each chromosome, it lists    ##
## those fragments that have at least 20 molecules aligned to the     ##
## reference.                                                         ##
########################################################################
MinMolecules <- 20
Chr <- 'chr7'
for(Chr in Chromosomes){
  Filename.In <- paste(OutputDataPath, Chr, '_Table.RData', sep='')
  load(Filename.In)
  chrTable <- get(x=paste(Chr, 'Table', sep='_'))
  Filename.Out <- paste(OutputDataPath, Chr, '_fragIndexList_Min', MinMolecules, '.txt', sep='')
  write.table(x=subset(chrTable, numMolecules>=MinMolecules)$refStartIndex, file=Filename.Out, 
              row.names=FALSE, col.names=FALSE)
  rm(chrTable)
}

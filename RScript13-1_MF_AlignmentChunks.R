rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script is for reading and extracting important information    ##
## out of any alignmentchunk file, produces after an alignment. This  ##
## file is copied from RScript11, and modified for mflorum data       ##
########################################################################

########################################################################
## Run Path definition file                                           ##
########################################################################
source('~/Project_GC_Content/RScripts_GC_Content/Paths_Header.R')
DataPath <- '/aspen/nandi/MF_cap348/maps_inca34/'
OutputDataPath <- '~/Project_GC_Content/RData/mflorum/'
AlignmentFilename <- 'MF_cap348_inca34_cf209_minSize50_minFrag5_Aligned.alignmentChunks'
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
AlChunk <- fn_readAlignmentChunks(Filename=Filename, Header=FALSE, Colnames=Colnames, 
                                  Colnames.Steve=Colnames.Steve)

str(AlChunk)
#AlChunk$refChr <- 'chr1'

Filename.Output <- paste(OutputDataPath, 'MF_cap348_inca34_cf209_minSize50_minFrag5_alignmentChunks.RData', sep='')
save(AlChunk, file=Filename.Output)

########################################################################
## Load the binary RData file of the Alignment Chunks, previously cr- ##
## eated                                                              ##
########################################################################
Filename.Bin <- paste(OutputDataPath, 
                      'MF_cap348_inca34_cf209_minSize50_minFrag5_alignmentChunks.RData', sep='')
load(Filename.Bin)

## How many unique molecules?
length(unique(AlChunk$molID))

## How many fragments in each chromosome?
#  Chr_Frags <- unique(AlChunk[,c('refChr', 'refStartIndex')])
#  head(Chr_Frags)
#  table(Chr_Frags$refChr)
## Nothing interesting here: just 1 chr, 39 fragments

########################################################################
## How many molecules are aligned to any given chr & frag index?      ##
########################################################################
Chr <- 'chr1'
TableName <- paste(Chr, '_Table', sep='')
Filename.Out <- paste(OutputDataPath, Chr, '_Table.RData', sep='')

ChrTable <- fn_numMolAlignedperLoc(AlChunk=subset(AlChunk, refChr==Chr))
assign(x=paste(Chr, '_Table', sep=''), value=ChrTable)
save(list=paste(Chr, '_Table', sep=''), file=Filename.Out)
#rm(ChrTable)
########################################################################


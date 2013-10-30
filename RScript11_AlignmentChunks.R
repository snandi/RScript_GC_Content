rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script is for plotting intensity profiles of multiple molecu- ##
########################################################################

########################################################################
## Initialize Header file and function library                        ##
########################################################################
FilePath <- '~/Project_GC_Content/RScripts_GC_Content/'
DataPath <- '/exports/aspen/steveg/human_nMaps/GC_content/'
OutputDataPath <- '~/Project_GC_Content/Data/'
Filename.Header <- paste('~/RScripts/HeaderFile_lmcg.R', sep='')
source(Filename.Header)
source(paste(FilePath, 'fn_Library_GC_Content.R', sep=''))
########################################################################

setwd(DataPath)
Filename <- paste(DataPath, 'alignmentChunks.withLength.all7134Groups.goldOnly', sep='')

Colnames <- c('refChr', 'refStartIndex', 'refEndIndex', 'molID', 'molStartIndex', 
                       'molEndIndex', 'refStartCoor', 'refEndCoord', 'molStartCoor', 
                       'molEndCoord', 'orient', 'lengthRatio')

Colnames.Steve <- c('refID', 'refStartIndex', 'refEndIndex', 'opID', 
'opStartIndex', 'opEndIndex', 'refStartCoor', 'refEndCoord', 'opStartCoor', 'opEndCoord', 
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
Chr_Frags <- unique(AlChunk[,c('refChr', 'refStartIndex', 'refEndIndex')])
head(Chr_Frags)
table(Chr_Frags$refChr)


rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script reads in DNA sequence data in FASTA format and conducts##
## some very basic analysis                                           ##
########################################################################

########################################################################
## Initialize Header file and function library                        ##
########################################################################
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
## Enter Fragment Details                                             ##
########################################################################
Chr <- 'chr16'
FragIndex <- 26
BasePairInterval <- 200     ## Length of base pair interval to estimate gcat %      
FragBP_Start <- 218358      ## Starting base pair coordinate of fragment
FragBP_End <- 224690        ## Ending base pair coordinate of fragment
NumBP_Frag <- FragBP_End - FragBP_Start ## Length of frag in BP
NumSubFrag <- NumBP_Frag/BasePairInterval ## Number of sub fragments
CpG_Start <- 222369         ## Starting base pair coordinate of CpG island
CpG_End <- 223447           ## Ending base pair coordinate of CpG island
NumBP_CpG <- CpG_End - CpG_Start ## Length of CpG island in BP
CpG_Start_Relative <- max(1, (CpG_Start - FragBP_Start + 1)/BasePairInterval)
                                 ## max is used for those CpG regions that start before
                                 ## the start of the fragment
CpG_End_Relative <- min(NumSubFrag, (CpG_End - FragBP_Start + 1)/BasePairInterval)
                                 ## min is used for those CpG regions that continue
                                 ## beyond the length of this fragment
########################################################################
## DNA sequence data                                                  ##
########################################################################
SeqComp <- fn_returnSeqComp(Chr=Chr, FragIndex=FragIndex, Interval=BasePairInterval, 
                            DataPath.CpG=DataPath.CpG)

SplitSeq_GCAT <- SeqComp[['SplitSeq_GCAT']]
SeqData <- SeqComp[['SeqData']]
SeqData.CpGIsland <- SeqData[(CpG_Start_Relative*BasePairInterval):(CpG_End_Relative*BasePairInterval)]
count(seq=SeqData, wordsize=2)
count(seq=SeqData.CpGIsland, wordsize=2)
count(seq=SeqData.CpGIsland, wordsize=1)


names(SeqComp)

my.Colors=c('olivedrab4', 'olivedrab1', 'gray25', 'gray56')
lattice::barchart(SplitSeq_GCAT, groups=TRUE, horizontal=FALSE,
                  xlab='Nucleotide Position',
                  scales=list(x=list(rot=90)), # To rotate the x-tick labels
                  main='Sequence Composition', 
                  par.settings = simpleTheme(col = my.Colors, border=my.Colors), 
                  # For changing colors in barplot and legend and no black 
                  # borders in the rectangles
                  auto.key = list(adj = 1, columns = 4, pch=c(15, 15,15, 15)), 
                  panel=function(...){
                    panel.barchart(...)
                    panel.abline(v=c(CpG_Start_Relative, CpG_End_Relative), col.line="red", lwd=2) 
                  })

lattice::barchart(count(seq=SeqData, wordsize=2), groups=TRUE, horizontal=TRUE,
                  #                  col=c('red', 'green', 'blue', 'yellow'),
                  main='Sequence Composition', 
                  auto.key = list(adj = 1))

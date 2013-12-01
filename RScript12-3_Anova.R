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

########################################################################
## Enter Fragment Details                                             ##
########################################################################
Filename <- paste(DataPath.CpG, 'CpG_IslandsOfInterest.csv', sep='')
Fragments <- read.csv(file=Filename, header=TRUE, stringsAsFactors=F)
FragmentRow <- 4
#for(FragmentRow in 3:8){
#Chr <- 'chr19'
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

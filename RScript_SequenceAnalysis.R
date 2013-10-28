rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script is for conducting further analysis on intensity and gc ##
## content. This script conducts a group level model, a model on mol- ##
## -ecule level aggregated model and produces an added variable plot  ##
########################################################################

########################################################################
## Initialize Header file and function library                        ##
########################################################################
FilePath <- '~/Project_GC_Content/RScripts_GC_Content/'
DataPath <- '~/Project_GC_Content/'
DataPath.Steve <- '/aspen/steveg/human_nMaps/GC_content/'
Filename.Header <- paste('~/RScripts/HeaderFile_lmcg.R', sep='')
source(Filename.Header)
source(paste(FilePath, 'fn_Library_GC_Content.R', sep=''))
########################################################################

Filename <- paste(DataPath, 'mesoplasma_florum_1.fasta', sep='')
MF_Seq <- read.fasta(file = Filename)[[1]]

length(MF_Seq)

table(MF_Seq)
barplot(table(MF_Seq))

GC(MF_Seq)

rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script extract genomic sequence data for multiple intervals   ##
## and plots them, to notice any apparent difference in them          ##
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


FragNames <- c('chr1_frag1_Length_68.RData', 
               'chr1_frag2_Length_263.RData',
               'chr1_frag3_Length_47.RData', 
               'chr1_frag6_Length_39.RData', 
               'chr1_frag9_Length_102.RData',
               'chr1_frag11_Length_80.RData', 
               'chr1_frag12_Length_55.RData',
               'chr1_frag15_Length_33.RData', 
               'chr1_frag18_Length_106.RData',
               'chr1_frag19_Length_163.RData',
               'chr1_frag21_Length_128.RData',
               'chr1_frag22_Length_71.RData',
               'chr1_frag24_Length_134.RData',
               'chr1_frag25_Length_178.RData',
               'chr1_frag27_Length_58.RData',
               'chr1_frag28_Length_55.RData',
               'chr1_frag30_Length_155.RData',
               'chr1_frag31_Length_68.RData',
               'chr1_frag32_Length_133.RData',
               'chr1_frag33_Length_80.RData',
               'chr1_frag35_Length_48.RData',
               'chr1_frag36_Length_226.RData',
               'chr1_frag37_Length_59.RData')

Filename.bploc <- paste(RDataPath, 'mflorum.bploc', sep='')
bp.loc <- read.table(file=Filename.bploc, header=T, sep='\t')
bp.loc$BasePairLength <- bp.loc$refMapCoordEnd - bp.loc$refMapCoordStart
bp.loc$PixelLength <- round(bp.loc$BasePairLength/209, 0)

TruncateLength <- 10

FragName <- FragNames[13]
pdf(file='Documentation/forMeeting_2014-03-26/PlotSeqComps.pdf', onefile=TRUE)

for(FragName in FragNames[c(6, 8, 10, 13)]){
  Chr <- unlist(strsplit(x=FragName, split='_'))[1]
  FragIndex <- gsub(pattern='frag', replacement='', 
                    x=unlist(strsplit(x=FragName, split='_'))[2])
  BasePairInterval <- 209     ## Length of base pair interval to estimate gcat %
  FragmentLength <- as.numeric(gsub(pattern='.RData', replacement='', 
                                    x=unlist(strsplit(x=FragName, split='_'))[4]))
  FragBP_Start <- bp.loc[which(bp.loc$alignedFragIndex == FragIndex), 'refMapCoordStart'] 
  ## Ending base pair coordinate of fragment
  FragBP_End <- bp.loc[which(bp.loc$alignedFragIndex == FragIndex), 'refMapCoordEnd']
  
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
  
  PlotGC <- lattice::barchart(SplitSeq_GCAT[11:50,], groups=TRUE, horizontal=FALSE,
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
                                #panel.abline(v=c(KeepFrom, KeepTo), col.line="red", lwd=2) 
                              }
  )
  print(PlotGC)
}
dev.off()


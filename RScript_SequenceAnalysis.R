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
DataPath.Steve <- '/aspen/steveg/human_nMaps/GC_content/'
Filename.Header <- paste('~/RScripts/HeaderFile_lmcg.R', sep='')
source(Filename.Header)
source(paste(RScriptPath, 'fn_Library_GC_Content.R', sep=''))
########################################################################

Filename <- paste(DataPath.Nandi, 'mesoplasma_florum_1.fasta', sep='')
MF_Seq <- read.fasta(file = Filename)[[1]]


length(MF_Seq)
count(seq=MF_Seq, wordsize=1)

table(MF_Seq)
barplot(table(MF_Seq))

GC(MF_Seq)

TestSeq <- MF_Seq[1:100]
BasePosition <- seq(from=1, to=length(TestSeq), by=10)
SplitSeq <- fn_subseq(x=TestSeq, n=10, force.number.of.groups=TRUE)
SplitSeq_ACGT <- do.call(what=rbind, lapply(X=SplitSeq, FUN=count, wordsize=1))
rownames(SplitSeq_ACGT) <- BasePosition[1:nrow(SplitSeq_ACGT)]
colnames(SplitSeq_ACGT) <- c('A', 'C', 'G', 'T')
SplitSeq_ACGT <- SplitSeq_ACGT/rowSums(SplitSeq_ACGT)

lattice::barchart(SplitSeq_ACGT, groups=TRUE, horizontal=FALSE,
#                  col=c('red', 'green', 'blue', 'yellow'),
                  xlab='Nucleotide Position',
                  main='Sequence Composition', 
                  auto.key = list(adj = 1))

TestSeq <- MF_Seq[1:10000]
BasePosition <- seq(from=1, to=length(TestSeq), by=200)
SplitSeq <- fn_subseq(x=TestSeq, n=1000, force.number.of.groups=TRUE)
SplitSeq_ACGT <- do.call(what=rbind, lapply(X=SplitSeq, FUN=count, wordsize=1))
rownames(SplitSeq_ACGT) <- BasePosition[1:nrow(SplitSeq_ACGT)]
colnames(SplitSeq_ACGT) <- c('A', 'C', 'G', 'T')
SplitSeq_ACGT <- SplitSeq_ACGT/rowSums(SplitSeq_ACGT)
SplitSeq_GCAT <- SplitSeq_ACGT[,c('G', 'C', 'A', 'T')]

lattice::barchart(SplitSeq_ACGT, groups=TRUE, horizontal=FALSE,
                  #                  col=c('red', 'green', 'blue', 'yellow'),
                  xlab='Nucleotide Position',
                  main='Sequence Composition', 
                  auto.key = list(adj = 1))

lattice::barchart(SplitSeq_GCAT, groups=TRUE, horizontal=FALSE,
                  #                  col=c('red', 'green', 'blue', 'yellow'),
                  xlab='Nucleotide Position',
                  scales=list(x=list(rot=90)),
                  main='Sequence Composition', 
                  auto.key = list(adj = 1))

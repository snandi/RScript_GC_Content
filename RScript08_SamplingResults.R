2rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script reads in the results of the 4500 group sampling models ##
## and summarizes the finals model parameters.                        ## 
########################################################################

########################################################################
## Initialize Header file and function library                        ##
########################################################################
FilePath <- '~/RA_Genomics/Project_GC_Content/RScripts_GC_Content/'
DataPath <- '~/RA_Genomics/Project_GC_Content/'
DataPath.Steve <- DataPath
Filename.Header <- paste('~/RScripts/HeaderFile_Nandi.R', sep='')
source(Filename.Header)
source(paste(FilePath, 'fn_Library_GC_Content.R', sep=''))
########################################################################

########################################################################
## Data Input                                                         ##
########################################################################
Filename <- paste(DataPath, 'SamplingOutput.RData', sep='')
load(Filename)

summary(MC_Output)
MC_Output <- as.data.frame(MC_Output)

lattice::densityplot(~ ModelA_RSquared + ModelB_RSquared + ModelC_RSquared, 
                     data=MC_Output, type='density', auto.key=TRUE, lty=c(1:3),
                     xlab='Model RSquares', lines=T, lwd=2)

summary(MC_Output[,c('ModelA_RSquared', 'ModelB_RSquared','ModelC_RSquared')])

lattice::densityplot(~ ModelA_Intercept_Beta + ModelB_Intercept_Beta + ModelC_Intercept_Beta, 
                     data=MC_Output, type='density', auto.key=TRUE, lty=c(1:3),
                     xlab='Model Intercept_Beta ', lines=T, lwd=2)

summary(MC_Output[,c('ModelA_Intercept_Beta', 'ModelB_Intercept_Beta','ModelC_Intercept_Beta')])

lattice::densityplot(~ ModelA_GC_Content_Beta + ModelB_GC_Content_Beta + ModelC_GC_Content_Beta, 
                     data=MC_Output, type='density', auto.key=list(lines=T, lty=1:3), lty=c(1:3),
                     xlab='Model GC_Content_Beta ', lines=T, lwd=2)
summary(MC_Output[,c('ModelA_GC_Content_Beta', 'ModelB_GC_Content_Beta','ModelC_GC_Content_Beta')])


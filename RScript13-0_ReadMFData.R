m(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script reads in optimal conversion factor data for MFloram    ##
## and plots density/histogram                                        ##
########################################################################

########################################################################
## Initialize Header file and function library                        ##
########################################################################
RScriptPath <- '~/Project_GC_Content/RScripts_GC_Content/'
RDataPath <- '~/Project_GC_Content/RData/'
Filename.Header <- paste('~/RScripts/HeaderFile_lmcg.R', sep='')
source(Filename.Header)
source(paste(RScriptPath, 'fn_Library_GC_Content.R', sep=''))
########################################################################

Filename.In.Asp <- '/aspen/nandi/MF_asp74/maps_inca34/optimalFactor_asp74_inca34_cf200'
OptConvFactor.Asp <- read.table(file=Filename.In.Asp, header=FALSE, sep='\t', stringsAsFactors=FALSE)
names(OptConvFactor.Asp) <- c('MoleculeID', 'ConversionFactor')
str(OptConvFactor.Asp)

summary(OptConvFactor.Asp)
Median.Asp <- median(OptConvFactor.Asp$ConversionFactor)
Filename.Out <- '/aspen/nandi/MF_asp74/maps_inca34/optimalFactor_asp74_inca34_cf200.pdf'
pdf(file=Filename.Out)
lattice::densityplot(~ ConversionFactor, data=OptConvFactor.Asp, type='density', 
                     xlab=paste('Optimal Conversion Factor, Median =', Median.Asp), lines=T, lwd=2, 
                     main='Mesoplasm - Asp (74 groups)', 
                     panel=function(...){
                       panel.densityplot(...)
                       panel.abline(v=Median.Asp, col.line="red", lwd=2) 
                     })
dev.off()

Filename.In.Cap <- '/aspen/nandi/MF_cap348/maps_inca34/optimalFactor_cap348_inca34_cf200'
OptConvFactor.Cap <- read.table(file=Filename.In.Cap, header=FALSE, sep='\t', stringsAsFactors=FALSE)
names(OptConvFactor.Cap) <- c('MoleculeID', 'ConversionFactor')
str(OptConvFactor.Cap)

summary(OptConvFactor.Cap)
Median.Cap <- median(OptConvFactor.Cap$ConversionFactor)
Filename.Out <- '/aspen/nandi/MF_cap348/maps_inca34/optimalFactor_cap348_inca34_cf200.pdf'
pdf(file=Filename.Out)
lattice::densityplot(~ ConversionFactor, data=OptConvFactor.Cap, type='density', 
                     xlab=paste('Optimal Conversion Factor, Median =', Median.Cap), lines=T, lwd=2, 
                     main='Mesoplasm - Cap (348 groups)', 
                     panel=function(...){
                       panel.densityplot(...)
                       panel.abline(v=Median.Cap, col.line="red", lwd=2) 
                     })

dev.off()


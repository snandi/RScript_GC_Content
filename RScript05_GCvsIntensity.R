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
Fragments <- fn_subsetByCoverage(CoverageNum=40, filename='FragmentData_5000.RData')
names(Fragments)

Fragment.Keep <- Fragments[['Fragment.Keep']]
FragmentData <- Fragments[['FragmentData']]

Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.5000Groups', sep='')
#Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.1061Groups', sep='')
#Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.50Groups', sep='')
Data <- read.table(file=Filename, header=TRUE, sep='\t', 
                   skip=0, stringsAsFactors=FALSE, fill=TRUE)
Data$groupID <- substr(x=Data$moleculeID, start=1, stop=7)
str(Data)
dim(Data)

Data <- merge(x=Data, y=Fragment.Keep, all=T)

Data <- subset(Data, Keep==1)

rm(Fragments, Fragment.Keep, FragmentData)

NumRows <- 5000000
########################################################################
## Data Input                                                         ##
########################################################################
## Model4: Y = log(intensity), X = GC content + number of fragments + number of pixels
## + groupID (Number of pixels is a measure of the molecule length)
## load('Model4.RData')
# summary(Model4)$coefficients[1:4,]
# summary(Model4)$r.squared

## jpeg(filename='AVPlot_fractionGC_Model4.jpg')
## avPlot(model=Model4, variable='fractionGC')
## dev.off()

## Next model is to run the same data as in Model4, but weighted least squares
## Model5: Y = log(intensity), X = GC content + number of fragments + number of pixels
## + groupID (Weighted least square)
Data$Weights <- Data$numMolecules / sum(Data$numMolecules)
## Time1 <- Sys.time()
ModelA <- lm( log(intensityPerPixel) ~ numFrags + numPixels + as.factor(groupID),
             data=Data)
##jpeg(filename='AVPlot_fractionGC_Model4.jpg')
##avPlot(model=ModelA, variable='fractionGC', )
pdf(file='AVPlot_fractionGC_Model4.pdf', onefile=TRUE)
##layout(matrix(1:2, ncol=1, byrow=TRUE))
plot(residuals(ModelA) ~ Data$fractionGC, xlab='fractionGC', ylab='Residuals',
     main='Relationship between Intensity and GC-content \n after controlling for group and other variables')
lines(lowess(x = Data$fractionGC, y=residuals(ModelA)), col='red')
abline(h=c(-0.10, 0.10), lty=2, col='blue')

plot(residuals(ModelA) ~ Data$fractionGC, ylim=c(-0.10, 0.10), xlab='fractionGC',
     ylab='Residuals', main='Same plot on a magnified scale')
lines(lowess(x = Data$fractionGC, y=residuals(ModelA)), col='red')
dev.off()

##layout(matrix(1:1))

## Model5 <- biglm::biglm( log(intensityPerPixel) ~ fractionGC + numFrags + numPixels +
##              as.factor(groupID), weights=~Weights,
##              data=Data[1:NumRows,])
## Time3 <- Sys.time()

## print(Time3 - Time2)
## ## summary(ModelA)$coefficients[1:10,]
## ## summary(ModelA)$r.squared

## save(Model5, file='Model5_biglm.RData')
## ## Summary5 <- summary(Model5)
## ## save(Summary5, file='Model5_biglm_Summary.RData')
## rm(Model5)
## quit(save = 'no')

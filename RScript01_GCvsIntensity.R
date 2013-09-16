rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script is to conduct preliminary statistical analyses on int- ##
## -ensity and gc content                                             ##
########################################################################

########################################################################
## Initialize Header file and function library                        ##
########################################################################
FilePath <- '~/Project_GC_Content/RScripts_GC_Content/'
DataPath <- '~/Project_GC_Content/'
Filename.Header <- paste('~/RScripts/HeaderFile_lmcg.R', sep='')
source(Filename.Header)
source(paste(FilePath, 'fn_Library_GC_Content.R', sep=''))
########################################################################

########################################################################
## Data Input                                                         ##
########################################################################
Filename <- paste(DataPath, 'intensityAndGCForAlignedFrags.txt', sep='')
Data <- read.table(file=Filename, header=TRUE, sep='\t', 
                   skip=0, stringsAsFactors=FALSE, fill=TRUE)
str(Data)
dim(Data)
#View(head(Data))
#View(Data)
########################################################################

########################################################################
## Basic Plots                                                        ##
########################################################################
#jpeg(filename='~/Project_GC_Content/Plot1_HistofRelativeIntensity.jpg')
plot(log(Data$relativeIntensity) ~ fractionGC, data=Data, cex=0.2, col='blue')
layout(matrix(1:2, ncol=1, byrow=TRUE))
hist((Data$relativeIntensity), breaks=40, xlab='Relative Intensity', 
     main='Histogram of relative intensity', cex.main=0.85)
hist(log(Data$relativeIntensity), breaks=200, xlab='log(Relative Intensity)', 
     main='Histogram of log transformed relative intensity', cex.main=0.8)
layout(matrix(1:1))
#dev.off()

#jpeg(filename='~/Project_GC_Content/Plot2_scatterplot.jpg')
## Plot code from Steve's script
xyplot(log(relativeIntensity) ~ fractionGC, data=Data,
       ylim=c(-0.75,0.75),
       xlim=c(.3,.7),
       pch=3,
       cex=.2,
       panel = function(...) {
         panel.xyplot(...)
         panel.loess(...,col="red")
       },
       xlab = "Fraction GC of the reference fragment",
       ylab = "log(Relative intensity) of aligned nMap fragment",
       main = "GC content vs relative intensity for 50 MM52 groups\n(red curve is loess smoothing)"
)
#dev.off()

########################################################################
## Some preliminary regressions                                       ##
########################################################################
Model1 <- lm(relativeIntensity ~ fractionGC, data=Data)
summary(Model1)
layout(matrix(1:4, ncol=2, byrow=TRUE))
plot(Model1)
layout(matrix(1:1))

Model2 <- lm(log(relativeIntensity) ~ fractionGC, data=Data)
summary(Model2)
layout(matrix(1:4, ncol=2, byrow=TRUE))
plot(Model2)
layout(matrix(1:1))

Model2.1 <- lm(log(intensity) ~ fractionGC, data=Data[-c(17243, 8552, 6329), ])
summary(Model2.1)
layout(matrix(1:4, ncol=2, byrow=TRUE))
plot(Model2.1)
layout(matrix(1:1))
X <- model.matrix(Model2.1)
AIC(object=Model2.1, k=log(dim(X)[2]))

Model3 <- lm(log(relativeIntensity) ~ fractionGC + as.factor(moleculeConversionFactor), 
             data=Data)
summary(Model3)$coefficients[1:2,]
summary(Model3)$r.squared
layout(matrix(1:4, ncol=2, byrow=TRUE))
plot(Model3)
layout(matrix(1:1))

Model3.1 <- lm(log(intensity) ~ fractionGC + as.factor(moleculeConversionFactor), 
             data=Data[-c(17243, 8552, 6329), ])
summary(Model3.1)$coefficients[1:2,]
summary(Model3.1)$r.squared
# layout(matrix(1:4, ncol=2, byrow=TRUE))
# plot(Model3.1)
# layout(matrix(1:1))
X <- model.matrix(Model3.1)
AIC(object=Model3.1, k=log(dim(X)[2]))

anova(Model2.1, Model3.1)

Model4 <- lm(log(relativeIntensity) ~ fractionGC + as.factor(moleculeConversionFactor), 
             data=Data[-c(17243, 8552, 6329),])
summary(Model4)$coefficients[1:2,]
summary(Model4)$r.squared
layout(matrix(1:4, ncol=2, byrow=TRUE))
plot(Model4)
layout(matrix(1:1))

Model5a <- lm( fractionGC ~ log(intensity), data=Data)
summary(Model5a)

Model5 <- lm( fractionGC ~ log(relativeIntensity) + as.factor(moleculeConversionFactor), 
              data=Data, )
summary(Model5)$coefficients[1:2,]
summary(Model5)$r.squared
layout(matrix(1:4, ncol=2, byrow=TRUE))
plot(Model5)
layout(matrix(1:1))
X <- model.matrix(Model5)
AIC(object=Model5, k=log(dim(X)[2]))

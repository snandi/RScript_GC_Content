rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script is for conducting further analysis on intensity and gc ##
## content. This will use data produced by Steve on 07/10 and after   ##
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

########################################################################
## Data Input                                                         ##
########################################################################
Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.50Groups', sep='')
#Filename <- paste(DataPath.Steve, 'intensityAndGCForAlignedFrags.50Groups', sep='')
Data <- read.table(file=Filename, header=TRUE, sep='\t', 
                   skip=0, stringsAsFactors=FALSE, fill=TRUE)
str(Data)
dim(Data)
#View(head(Data))
#View(Data)
D#ata <- Data[1:100000,]
########################################################################

########################################################################
## Some Molecule level summary                                        ##
########################################################################
length(unique(Data$moleculeID))
Data.Molecule <- unique(Data[,c('moleculeID', 'moleculeConversionFactor', 'numFrags', 
                                'numPixels', 'totalIntensity')])
length(unique(Data.Molecule$moleculeID))
Corr.Molecule <- corr.test(Data.Molecule[,-1])

########################################################################
## Regression Models                                                  ##
########################################################################

########################################################################
## Model0: Y = GC content, X = log(intensity)
Model0 <- lm( fractionGC ~ log(intensityPerPixel), data=Data)
summary(Model0)
X <- model.matrix(Model0)
AIC(object=Model0, k=log(dim(X)[2]))

########################################################################
Time1 <- Sys.time()
## Model1: Y = GC content, X = log(intensity) + molecule level random effect
Model1 <- lm( fractionGC ~ log(intensityPerPixel) + as.factor(moleculeID), 
              data=Data)
print(paste('Time for fitting Model1:', Sys.time() - Time1))
summary(Model1)$coefficients[1:2,]
summary(Model1)$r.squared
# layout(matrix(1:4, ncol=2, byrow=TRUE))
# plot(Model1)
# layout(matrix(1:1))
X <- model.matrix(Model1)
AIC(object=Model1, k=log(dim(X)[2]))

plot(residuals(object=Model1) ~ Data$moleculeConversionFactor, main='Residual of Model1 vs 
     Molecule Conversion Factor', xlab='Molecule Conversion Factor', ylab='Residuals')
lines(lowess(y=residuals(object=Model1), x=Data$moleculeConversionFactor), col="Red")

########################################################################
## Model2: Y = GC content, X = log(intensity) + molecule level random effect + 
## moleculeConverstionFactor (as a cont variable)
hist(Data$moleculeConversionFactor, breaks=40)
Model2 <- lm( fractionGC ~ log(intensityPerPixel) + moleculeConversionFactor + 
                as.factor(moleculeID), data=Data)
summary(Model2)$coefficients[1:3,]
summary(Model2)$r.squared
X <- model.matrix(Model2)
AIC(object=Model2, k=log(dim(X)[2]))

anova(Model1, Model2)
## comments on Model2: R^2, BIC, anova confirm the addition of this new variable 
## moleculeConverstionFactor
plot(residuals(object=Model2) ~ log(Data$numFrags), main='Residual of Model2 vs 
     Number of Fragments', xlab='log(Number of Fragments)', ylab='Residuals')
lines(lowess(y=residuals(object=Model2), x=log(Data$numFrags)), col="Red")

########################################################################
## Model3: Y = GC content, X = log(intensity) + moleculeConverstionFactor 
## (as a cont variable)
Model3 <- lm( fractionGC ~ log(intensityPerPixel) + moleculeConversionFactor, 
                data=Data)
summary(Model3)
X <- model.matrix(Model3)
AIC(object=Model3, k=log(dim(X)[2]))

anova(Model0, Model3)

########################################################################
## Model4: Y = GC content, X = log(intensity) + moleculeConverstionFactor +
## numFrags (another molecule level variable)
Model4 <- lm( fractionGC ~ log(intensityPerPixel) + moleculeConversionFactor +
                numFrags, data=Data)
summary(Model4)
X <- model.matrix(Model4)
AIC(object=Model4, k=log(dim(X)[2]))
vif(object=Model4)
anova(Model3, Model4)

########################################################################
## Model5: Y = GC content, X = log(intensity) + moleculeConverstionFactor +
## numFrags + numPixels (measure of molecule length)
Model5 <- lm( fractionGC ~ log(intensityPerPixel) + moleculeConversionFactor +
                numFrags + numPixels, data=Data)
summary(Model5)
X <- model.matrix(Model5)
AIC(object=Model5, k=log(dim(X)[2]))
vif(object=Model5)
anova(Model4, Model5)

########################################################################
## Model6: Y = GC content, X = log(intensity) + moleculeConverstionFactor +
## numFrags + numPixels + molecule level random effect
Model6 <- lm( fractionGC ~ log(intensityPerPixel) + moleculeConversionFactor +
                numFrags + numPixels + as.factor(moleculeID), data=Data)
summary(Model6)$coefficients[1:5,]
summary(Model6)$r.squared
X <- model.matrix(Model6)
AIC(object=Model6, k=log(dim(X)[2]))
VIF.Model6 <- vif(object=Model6)
VIF.Model6[1:5]
anova(Model5, Model6)
anova(Model1, Model6)


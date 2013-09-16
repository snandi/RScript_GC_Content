rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script is for conducting further analysis on intensity and gc ##
## content. This will use data produced by RScript03.This script will ##
## incorporate fragment level effect first, as per discussion during  ##
## the meeting with M. Newton and D. Schwartz on 07/15                ##
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

########################################################################
## Regression Models on Fragment Data                                 ##
########################################################################
Model1.Frag <- lm(log(intensityPerPixel_mean) ~ fractionGC, data=FragmentData)
summary(Model1.Frag)
#plot(Model1.Frag)
anova(Model1.Frag)

#plot(log(intensityPerPixel_mean) ~ fractionGC, data=FragmentData.Subset)

Model2.Frag <- lm(log(intensityPerPixel_mean) ~ fractionGC + as.factor(alignedChr), 
             data=FragmentData)
summary(Model2.Frag)
# layout(matrix(1:4, ncol=2, byrow=TRUE))
# plot(Model2.Frag)
# layout(matrix(1:1))

anova(Model1.Frag, Model2.Frag)

########################################################################
## Regression Models on the whole Data                                ##
########################################################################
## Model1: Y = log(intensity), X = GC content
Model1 <- lm( log(intensityPerPixel) ~ fractionGC, data=Data)
save(Model1, file='Model1.RData')
Summary1 <- summary(Model1)
save(Summary1, file='Model1_Summary.RData')
X <- model.matrix(Model1)
AIC1 <- AIC(object=Model1, k=log(dim(X)[2]))

## Model2: Y = log(intensity), X = GC content + number of fragments
Time1 <- Sys.time()
Model2 <- lm( log(intensityPerPixel) ~ fractionGC + numFrags, 
              data=Data)
print(Sys.time() - Time1)
save(Model2, file='Model2.RData')
Summary2 <- summary(Model2)
save(Summary2, file='Model2_Summary.RData')
Summary2$coefficients[1:2,]
Summary2$r.squared

layout(matrix(1:4, ncol=2, byrow=TRUE))
#plot(Model2)
#layout(matrix(1:1))
#X <- model.matrix(Model2)
AIC(object=Model2, k=log(dim(X)[2]))
Anova2 <- anova(Model1, Model2)
save(Anova2, file='Model2_Anova.RData')

## Model3: Y = log(intensity), X = GC content + number of fragments + groupID
Time1 <- Sys.time()
Model3 <- lm( log(intensityPerPixel) ~ fractionGC + numFrags + as.factor(groupID), 
              data=Data)
print(Sys.time() - Time1)
save(Model3, file='Model3.RData')
Summary3 <- summary(Model3)
save(Summary3, file='Model3_Summary.RData')
Summary3$coefficients[1:4,]
Summary3$r.squared
#layout(matrix(1:4, ncol=2, byrow=TRUE))
#plot(Model3)
#layout(matrix(1:1))
X <- model.matrix(Model3)
AIC(object=Model3, k=log(dim(X)[2]))
Anova3 <- anova(Model2, Model3)
save(Anova3, file='Model3_Anova.RData')

## Model4: Y = log(intensity), X = GC content + number of fragments + number of pixels
## + groupID (Number of pixels is a measure of the molecule length)
Time1 <- Sys.time()
Model4 <- lm( log(intensityPerPixel) ~ fractionGC + numFrags + numPixels +
             as.factor(groupID),
             data=Data)
print(Sys.time() - Time1)
save(Model4, file='Model4.RData')
summary(Model4)$coefficients[1:4,]
summary(Model4)$r.squared
#layout(matrix(1:4, ncol=2, byrow=TRUE))
#plot(Model4)
#layout(matrix(1:1))
X <- model.matrix(Model4)
AIC(object=Model4, k=log(dim(X)[2]))
Anova4 <- anova(Model3, Model4)
save(Anova4, file='Model4_Anova.RData')


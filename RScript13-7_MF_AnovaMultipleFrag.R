rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script performs a within fragment and between fragment anova, ##
## between muliple fragments. This is an extension of RScript13-6     ##
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
I <- 1
Data <- NULL
MaxPixel <- 49

for(I in 1:length(FragNames)){
  Filename <- paste(RDataPath, FragNames[I], sep='')
  load(Filename)
  Data.Temp <- IntensityData.Aligned
  Data.Temp$Fragment <- unlist(strsplit(x=FragNames[I], split="_"))[2]
  Data <- rbind(Data, Data.Temp)
  rm(IntensityData.Aligned, Data.Temp)
}


Data$GroupID <- substr(Data$MoleculeID, 
                       start=1, stop=7)
GroupVariance <- as.data.frame(aggregate(x=Data$Intensity_Normalized, by=list(Data$GroupID), 
                                         FUN=sd))
names(GroupVariance) <- c('GroupID', 'GroupVariance')
Data <- merge(x=Data, y=GroupVariance, by='GroupID')

Data <- subset(Data, PixelNum <= MaxPixel)

Data <- within(data=Data, {
  Intensity_Weighted <- Intensity_Normalized/GroupVariance
  PixelFactor <- as.factor(PixelNum)
  GroupID <- as.factor(GroupID)
  Fragment <- as.factor(Fragment)
})
str(Data)

#lattice::histogram(~ Intensity_Normalized | Fragment, data = Data, breaks=20)
#lattice::histogram(~ log(Intensity_Normalized) | Fragment, data = Data, breaks=20)
#lattice::histogram(~ log(Intensity_Weighted) | Fragment, data = Data, breaks=20)

Model1.Norm <- lm(Intensity_Normalized ~ Fragment + PixelFactor, data=Data)
pValues1.Norm <- coef(summary(Model1.Norm))[,4]
anova(Model1.Norm)
summary(Model1.Norm)

Model1 <- lm(log(Intensity_Normalized) ~ Fragment + PixelFactor, data=Data)
pValues1 <- coef(summary(Model1))[,4]
anova(Model1)
#summary(Model1)
xtable(anova(Model1))

Model2 <- lm(log(Intensity_Normalized) ~ Fragment + PixelFactor + GroupID, 
                      data=Data)
pValues2 <- coef(summary(Model2))[,4]
anova(Model2)
xtable(anova(Model2))

Tukey2 <- TukeyHSD.aov(aov(log(Intensity_Normalized) ~ Fragment + PixelFactor, 
                           data=Data), conf.level=0.95)

Tukey2.pValues <- fn_return_pValueTukeyMatrix(TukeyObj=Tukey2)

Plot1 <- levelplot(Tukey2.pValues, cex.axis=2, cex.lab=2,
          col.regions=colorRampPalette(c("red", "green")), 
          scales=list(x=list(rot=90)),
          main="Tukey test of pairwise difference in fragments")
jpeg(filename='Documentation/forMeeting_2014-03-26/Plot1.jpg')
print(Plot1)
dev.off()
pdf(file='Documentation/forMeeting_2014-03-26/Plot1.pdf')
print(Plot1)
dev.off()

print(as.data.frame(Tukey2$Fragment))

Model3 <- lm(log(Intensity_Weighted) ~ Fragment + PixelFactor, 
                      data=Data)
anova(Model3)
xtable(anova(Model3))
summary(Model3)
Tukey3 <- TukeyHSD.aov(aov(log(Intensity_Weighted) ~ Fragment + PixelFactor, 
                           data=Data), conf.level=0.95)
str(Tukey3$Fragment)
print(Tukey3$Fragment)
Tukey3.pValues <- fn_return_pValueTukeyMatrix(TukeyObj=Tukey3)
Plot2 <- levelplot(Tukey3.pValues, cex.axis=2, cex.lab=2,
          col.regions=colorRampPalette(c("red", "green")), 
          scales=list(x=list(rot=90)),
          main="Tukey pValues of pairwise diff between intervals (normalized by group variability) ")
jpeg(filename='Documentation/forMeeting_2014-03-26/Plot2.jpg')
print(Plot2)
dev.off()

pdf(file='Documentation/forMeeting_2014-03-26/Plot2.pdf')
print(Plot2)
dev.off()


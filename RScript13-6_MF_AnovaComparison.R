rm(list=ls(all.names=TRUE))
rm(list=objects(all.names=TRUE))
#dev.off()

########################################################################
## This script performs a within fragment and between fragment anova, ##
## only between two fragments, frag33 and frag11. 
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

Filename1 <- paste(RDataPath, 'chr1_frag33_Length_82.RData', sep='')
load(Filename1)
Frag33 <- IntensityData.Aligned
Frag33$Fragment <- 'Frag33'
rm(IntensityData.Aligned)

Filename2 <- paste(RDataPath, 'chr1_frag11_Length_82.RData', sep='')
load(Filename2)
Frag11 <- IntensityData.Aligned
Frag11$Fragment <- 'Frag11'
rm(IntensityData.Aligned)

Data <- rbind(Frag11, Frag33)

Data$GroupID <- substr(Data$MoleculeID, 
                       start=1, stop=7)
Data <- within(data=Data, {
  PixelFactor <- as.factor(PixelNum)
  GroupID <- as.factor(GroupID)
  Fragment <- as.factor(Fragment)
})
str(Data)

Model1 <- lm(Intensity ~ Fragment + PixelFactor, data=Data)
pValues.1 <- coef(summary(Model1))[,4]
anova(Model1)
summary(Model1)

Model1.Norm <- lm(Intensity_Normalized ~ Fragment + PixelFactor, data=Data)
pValues1.Norm <- coef(summary(Model1.Norm))[,4]
anova(Model1.Norm)

#summary(Model1.Norm)

Model1.Norm.Log <- lm(log(Intensity_Normalized) ~ Fragment + PixelFactor, data=Data)
pValues1.Norm.Log <- coef(summary(Model1.Norm.Log))[,4]
anova(Model1.Norm.Log)
summary(Model1.Norm.Log)

xtable(anova(Model1.Norm.Log))

pValue.Raw[Index] <- round(anova(Model1.Norm)$'Pr(>F)'[1], 4)
pValue.Log[Index] <- round(anova(Model1.Norm.Log)$'Pr(>F)'[1], 4)
# Model_RSS[Index] <- round(summary(Model1.Norm)$sigma^2, 4)
# Model_RSquared[Index] <- round(summary(Model1.Norm)$r.squared, 4)

Kruskal_Wallis <- kruskal.test(log(Intensity_Normalized) ~ PixelFactor, 
                               data=Data)
pValue.KW[Index] <- round(Kruskal_Wallis$p.value, 4)

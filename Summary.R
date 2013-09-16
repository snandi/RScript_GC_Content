load('Model1.RData')
load('Model2.RData')
load('Model3.RData')
load('Model4.RData')

Summary1 <- summary(Model1)
save(Summary1, file='Model1_Summary.RData')

Summary2 <- summary(Model2)
save(Summary2, file='Model2_Summary.RData')
Anova2 <- anova(Model1, Model2)
save(Anova2, file='Model2_Anova.RData')

Summary3 <- summary(Model3)
save(Summary3, file='Model3_Summary.RData')
Anova3 <- anova(Model2, Model3)
save(Anova3, file='Model3_Anova.RData')

load('Model4.RData')
Summary4 <- summary(Model4)
save(Summary4, file='Model4_Summary.RData')
Anova4 <- anova(Model3, Model4)
save(Anova4, file='Model4_Anova.RData')


quit(save = 'no')

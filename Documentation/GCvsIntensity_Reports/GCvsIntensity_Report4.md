Fluorescence Intensity & Sequence Composition
=============================================
### Author: Subhrangshu Nandi
### Date: 2013-11-26

### Big Picture (Goal)

- Establish relationship between Fluorescence Intensity & Sequence Composition
- Construct acceptable band of intensities for as many locations on the genome as possible
- Use data from Multiple Miloma Sample with 7,134 groups

### Today's objectives

1. Discuss alignment algorithm
  - Needs improvement
2. Discuss anova result of genomic position & intensity
3. Discuss temporary detour
  - CpG Islands
  - Effect of Methylation on Fluorescence Intensity
  - Plots of CpG Islands
  - Plots of Alpha-globin CpG Islands
  - Statistical test results
4. Discuss getting back on track
  - Next steps




#### CpG Islands

- 6 separate CpG islands chosen by Adi (chr8, chr2, chr19, chr20, chr3, chr7)
- Example: chr20_frag6217


```
## 'data.frame':	430 obs. of  11 variables:
##  $ PixelNum            : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ Intensity           : num  32252 31144 57532 38549 49398 ...
##  $ MoleculeID          : Factor w/ 10 levels "2329767_734_2060207",..: 6 5 9 8 2 1 10 7 4 3 ...
##  $ Intensity_Normalized: num  0.961 0.949 1.234 1.108 0.902 ...
##  $ PixelFactor         : Factor w/ 45 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ GG                  : num  0.305 0.305 0.305 0.305 0.305 ...
##  $ CC                  : num  0.218 0.218 0.218 0.218 0.218 ...
##  $ AA                  : num  0.223 0.223 0.223 0.223 0.223 ...
##  $ TT                  : num  0.254 0.254 0.254 0.254 0.254 ...
##  $ CpG                 : Factor w/ 2 levels "0","1": 2 2 2 2 2 2 2 2 2 2 ...
##  $ GroupID             : Factor w/ 10 levels "2329767","2332829",..: 6 5 9 8 2 1 10 7 4 3 ...
```

```
##   PixelNum Intensity          MoleculeID Intensity_Normalized PixelFactor
## 1        1     32252 2392908_734_2060517               0.9610           1
## 2        1     31144 2391890_734_2060406               0.9486           1
## 3        1     57532 2408760_734_2060630               1.2341           1
## 4        1     38549 2400457_734_2060369               1.1081           1
## 5        1     49398 2332829_734_2060762               0.9020           1
## 6        1     41426 2329767_734_2060207               0.9448           1
##       GG     CC     AA     TT CpG GroupID
## 1 0.3046 0.2183 0.2234 0.2538   1 2392908
## 2 0.3046 0.2183 0.2234 0.2538   1 2391890
## 3 0.3046 0.2183 0.2234 0.2538   1 2408760
## 4 0.3046 0.2183 0.2234 0.2538   1 2400457
## 5 0.3046 0.2183 0.2234 0.2538   1 2332829
## 6 0.3046 0.2183 0.2234 0.2538   1 2329767
```

#### Model 1: Intensity ~ Genome location (as factor)

```
## Analysis of Variance Table
## 
## Response: Intensity
##              Df   Sum Sq  Mean Sq F value Pr(>F)
## PixelFactor  42 8.58e+08 20416984    0.23      1
## Residuals   387 3.38e+10 87328789
```

```
## 
## Call:
## lm(formula = Intensity ~ PixelFactor, data = IntensityData.Reg)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -11989  -6837  -3567   7206  23249 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    39626.0     2955.1   13.41   <2e-16 ***
## PixelFactor2     494.2     4179.2    0.12     0.91    
## PixelFactor3    1136.7     4179.2    0.27     0.79    
## PixelFactor4    1663.7     4179.2    0.40     0.69    
## PixelFactor5    2076.7     4179.2    0.50     0.62    
## PixelFactor6    2328.6     4179.2    0.56     0.58    
## PixelFactor7    2472.7     4179.2    0.59     0.55    
## PixelFactor8    2545.0     4179.2    0.61     0.54    
## PixelFactor9    2479.1     4179.2    0.59     0.55    
## PixelFactor10   2127.9     4179.2    0.51     0.61    
## PixelFactor11   2223.1     4179.2    0.53     0.60    
## PixelFactor12   2665.2     4179.2    0.64     0.52    
## PixelFactor13   2691.5     4179.2    0.64     0.52    
## PixelFactor14   2776.9     4179.2    0.66     0.51    
## PixelFactor15   2147.7     4179.2    0.51     0.61    
## PixelFactor16   1761.1     4179.2    0.42     0.67    
## PixelFactor17   1520.5     4179.2    0.36     0.72    
## PixelFactor18   1488.3     4179.2    0.36     0.72    
## PixelFactor19   1556.7     4179.2    0.37     0.71    
## PixelFactor20   1622.5     4179.2    0.39     0.70    
## PixelFactor21   1352.5     4179.2    0.32     0.75    
## PixelFactor22    725.2     4179.2    0.17     0.86    
## PixelFactor23    797.5     4179.2    0.19     0.85    
## PixelFactor24    118.5     4179.2    0.03     0.98    
## PixelFactor25     40.0     4179.2    0.01     0.99    
## PixelFactor26    -10.2     4179.2    0.00     1.00    
## PixelFactor27   -283.1     4179.2   -0.07     0.95    
## PixelFactor28    -51.0     4179.2   -0.01     0.99    
## PixelFactor29    265.7     4179.2    0.06     0.95    
## PixelFactor30    -16.7     4179.2    0.00     1.00    
## PixelFactor31   -181.1     4179.2   -0.04     0.97    
## PixelFactor32   -312.6     4179.2   -0.07     0.94    
## PixelFactor33   -371.7     4179.2   -0.09     0.93    
## PixelFactor34   -355.3     4179.2   -0.09     0.93    
## PixelFactor35   -358.1     4179.2   -0.09     0.93    
## PixelFactor36   -379.2     4179.2   -0.09     0.93    
## PixelFactor37   -492.8     4179.2   -0.12     0.91    
## PixelFactor38   -795.7     4179.2   -0.19     0.85    
## PixelFactor39  -1179.1     4179.2   -0.28     0.78    
## PixelFactor40  -1795.5     4179.2   -0.43     0.67    
## PixelFactor41  -1914.4     4179.2   -0.46     0.65    
## PixelFactor42  -2232.4     4179.2   -0.53     0.59    
## PixelFactor43  -2159.4     4179.2   -0.52     0.61    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 9340 on 387 degrees of freedom
## Multiple R-squared: 0.0247,	Adjusted R-squared: -0.0811 
## F-statistic: 0.234 on 42 and 387 DF,  p-value: 1
```


#### Model 3: Intensity ~ Genome location (as factor) + GroupID

```
## Analysis of Variance Table
## 
## Response: Intensity
##              Df   Sum Sq  Mean Sq F value Pr(>F)    
## GroupID       9 2.88e+10 3.20e+09  239.57 <2e-16 ***
## PixelFactor  42 8.58e+08 2.04e+07    1.53  0.022 *  
## Residuals   378 5.04e+09 1.33e+07                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## Analysis of Variance Table
## 
## Model 1: Intensity ~ PixelFactor
## Model 2: Intensity ~ GroupID + PixelFactor
##   Res.Df      RSS Df Sum of Sq   F Pr(>F)    
## 1    387 3.38e+10                            
## 2    378 5.04e+09  9  2.88e+10 240 <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## 
## Call:
## lm(formula = Intensity ~ GroupID + PixelFactor, data = IntensityData.Reg)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -10189  -1464     -2   1680  14113 
## 
## Coefficients:
##                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     42833.0     1270.0   33.73  < 2e-16 ***
## GroupID2332829  10529.4      787.6   13.37  < 2e-16 ***
## GroupID2336206   8189.4      787.6   10.40  < 2e-16 ***
## GroupID2390555  -6187.8      787.6   -7.86  4.1e-14 ***
## GroupID2391890 -10391.9      787.6  -13.19  < 2e-16 ***
## GroupID2392908  -9993.7      787.6  -12.69  < 2e-16 ***
## GroupID2398040 -10065.5      787.6  -12.78  < 2e-16 ***
## GroupID2400457  -7868.9      787.6   -9.99  < 2e-16 ***
## GroupID2408760   5929.2      787.6    7.53  3.8e-13 ***
## GroupID2428567 -12210.4      787.6  -15.50  < 2e-16 ***
## PixelFactor2      494.2     1633.2    0.30     0.76    
## PixelFactor3     1136.7     1633.2    0.70     0.49    
## PixelFactor4     1663.7     1633.2    1.02     0.31    
## PixelFactor5     2076.7     1633.2    1.27     0.20    
## PixelFactor6     2328.6     1633.2    1.43     0.15    
## PixelFactor7     2472.8     1633.2    1.51     0.13    
## PixelFactor8     2545.0     1633.2    1.56     0.12    
## PixelFactor9     2479.1     1633.2    1.52     0.13    
## PixelFactor10    2127.9     1633.2    1.30     0.19    
## PixelFactor11    2223.1     1633.2    1.36     0.17    
## PixelFactor12    2665.2     1633.2    1.63     0.10    
## PixelFactor13    2691.5     1633.2    1.65     0.10    
## PixelFactor14    2776.9     1633.2    1.70     0.09 .  
## PixelFactor15    2147.7     1633.2    1.32     0.19    
## PixelFactor16    1761.1     1633.2    1.08     0.28    
## PixelFactor17    1520.5     1633.2    0.93     0.35    
## PixelFactor18    1488.3     1633.2    0.91     0.36    
## PixelFactor19    1556.7     1633.2    0.95     0.34    
## PixelFactor20    1622.5     1633.2    0.99     0.32    
## PixelFactor21    1352.5     1633.2    0.83     0.41    
## PixelFactor22     725.2     1633.2    0.44     0.66    
## PixelFactor23     797.6     1633.2    0.49     0.63    
## PixelFactor24     118.5     1633.2    0.07     0.94    
## PixelFactor25      40.0     1633.2    0.02     0.98    
## PixelFactor26     -10.2     1633.2   -0.01     1.00    
## PixelFactor27    -283.1     1633.2   -0.17     0.86    
## PixelFactor28     -51.0     1633.2   -0.03     0.98    
## PixelFactor29     265.7     1633.2    0.16     0.87    
## PixelFactor30     -16.7     1633.2   -0.01     0.99    
## PixelFactor31    -181.1     1633.2   -0.11     0.91    
## PixelFactor32    -312.5     1633.2   -0.19     0.85    
## PixelFactor33    -371.7     1633.2   -0.23     0.82    
## PixelFactor34    -355.3     1633.2   -0.22     0.83    
## PixelFactor35    -358.1     1633.2   -0.22     0.83    
## PixelFactor36    -379.2     1633.2   -0.23     0.82    
## PixelFactor37    -492.8     1633.2   -0.30     0.76    
## PixelFactor38    -795.7     1633.2   -0.49     0.63    
## PixelFactor39   -1179.1     1633.2   -0.72     0.47    
## PixelFactor40   -1795.5     1633.2   -1.10     0.27    
## PixelFactor41   -1914.4     1633.2   -1.17     0.24    
## PixelFactor42   -2232.4     1633.2   -1.37     0.17    
## PixelFactor43   -2159.4     1633.2   -1.32     0.19    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 3650 on 378 degrees of freedom
## Multiple R-squared: 0.855,	Adjusted R-squared: 0.835 
## F-statistic: 43.5 on 51 and 378 DF,  p-value: <2e-16
```


#### Model 4: Intensity ~ GroupID + G + C + T

```
## 
## Call:
## lm(formula = Intensity ~ GroupID + GG + CC + TT, data = IntensityData.Reg)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -10780  -1400    146   1304  14596 
## 
## Coefficients:
##                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       42664       3371   12.65  < 2e-16 ***
## GroupID2332829    10529        797   13.21  < 2e-16 ***
## GroupID2336206     8189        797   10.28  < 2e-16 ***
## GroupID2390555    -6188        797   -7.77  6.4e-14 ***
## GroupID2391890   -10392        797  -13.04  < 2e-16 ***
## GroupID2392908    -9994        797  -12.54  < 2e-16 ***
## GroupID2398040   -10066        797  -12.63  < 2e-16 ***
## GroupID2400457    -7869        797   -9.87  < 2e-16 ***
## GroupID2408760     5929        797    7.44  5.8e-13 ***
## GroupID2428567   -12210        797  -15.32  < 2e-16 ***
## GG                 4302       4672    0.92     0.36    
## CC                 2282       3524    0.65     0.52    
## TT                -4678       6109   -0.77     0.44    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 3690 on 417 degrees of freedom
## Multiple R-squared: 0.836,	Adjusted R-squared: 0.831 
## F-statistic:  177 on 12 and 417 DF,  p-value: <2e-16
```

```
## Analysis of Variance Table
## 
## Model 1: Intensity ~ GroupID
## Model 2: Intensity ~ GroupID + GG + CC + TT
##   Res.Df      RSS Df Sum of Sq    F Pr(>F)   
## 1    420 5.90e+09                            
## 2    417 5.69e+09  3  2.06e+08 5.03  0.002 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


#### Model 5: Intensity ~ GroupID + CpG (Yes/No)

```
## 
## Call:
## lm(formula = Intensity ~ GroupID + CpG, data = IntensityData.Reg)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -10271  -1519     83   1526  14067 
## 
## Coefficients:
##                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       41937        601   69.73  < 2e-16 ***
## GroupID2332829    10529        774   13.60  < 2e-16 ***
## GroupID2336206     8189        774   10.58  < 2e-16 ***
## GroupID2390555    -6188        774   -7.99  1.3e-14 ***
## GroupID2391890   -10392        774  -13.42  < 2e-16 ***
## GroupID2392908    -9994        774  -12.91  < 2e-16 ***
## GroupID2398040   -10066        774  -13.00  < 2e-16 ***
## GroupID2400457    -7869        774  -10.16  < 2e-16 ***
## GroupID2408760     5929        774    7.66  1.3e-13 ***
## GroupID2428567   -12210        774  -15.77  < 2e-16 ***
## CpG1               2301        369    6.23  1.2e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 3590 on 419 degrees of freedom
## Multiple R-squared: 0.844,	Adjusted R-squared: 0.84 
## F-statistic:  227 on 10 and 419 DF,  p-value: <2e-16
```

```
## Analysis of Variance Table
## 
## Model 1: Intensity ~ GroupID
## Model 2: Intensity ~ GroupID + CpG
##   Res.Df     RSS Df Sum of Sq    F  Pr(>F)    
## 1    420 5.9e+09                              
## 2    419 5.4e+09  1     5e+08 38.8 1.2e-09 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


#### Comments

- Need more accurate molecule aligning algorithm
- Genomic position not significant for this fragment
- Need to analyze multiple molecules from the same groups, in order to better control for the experimental artifacts
- Methylated or not is a random variable (not sure how to proceed)

#### See other plots

#### Example: chr16_frag27 (Alpha-Globin CpG Island)


```
## 'data.frame':	516 obs. of  11 variables:
##  $ PixelNum            : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ Intensity           : num  48195 40728 34151 45641 31794 ...
##  $ MoleculeID          : Factor w/ 12 levels "2331094_734_2060343",..: 2 6 3 1 12 9 11 5 4 8 ...
##  $ Intensity_Normalized: num  0.961 1.022 0.97 0.881 0.972 ...
##  $ PixelFactor         : Factor w/ 43 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ GG                  : num  0.221 0.221 0.221 0.221 0.221 ...
##  $ CC                  : num  0.251 0.251 0.251 0.251 0.251 ...
##  $ AA                  : num  0.327 0.327 0.327 0.327 0.327 ...
##  $ TT                  : num  0.201 0.201 0.201 0.201 0.201 ...
##  $ CpG                 : Factor w/ 2 levels "0","1": 1 1 1 1 1 1 1 1 1 1 ...
##  $ GroupID             : Factor w/ 12 levels "2331094","2333529",..: 2 6 3 1 12 9 11 5 4 8 ...
```

```
##   PixelNum Intensity          MoleculeID Intensity_Normalized PixelFactor
## 1        1     48195 2333529_734_2060193               0.9610           1
## 2        1     40728 2391723_734_2060650               1.0218           1
## 3        1     34151 2374503_734_2061370               0.9702           1
## 4        1     45641 2331094_734_2060343               0.8814           1
## 5        1     31794 2401144_734_2060461               0.9718           1
## 6        1     33446 2392825_734_2060550               0.9466           1
##       GG     CC     AA    TT CpG GroupID
## 1 0.2211 0.2513 0.3266 0.201   0 2333529
## 2 0.2211 0.2513 0.3266 0.201   0 2391723
## 3 0.2211 0.2513 0.3266 0.201   0 2374503
## 4 0.2211 0.2513 0.3266 0.201   0 2331094
## 5 0.2211 0.2513 0.3266 0.201   0 2401144
## 6 0.2211 0.2513 0.3266 0.201   0 2392825
```

#### Model 1: Intensity ~ Genome location (as factor)

```
## Analysis of Variance Table
## 
## Response: Intensity
##              Df   Sum Sq  Mean Sq F value Pr(>F)
## PixelFactor  42 6.18e+08 14717068    0.31      1
## Residuals   473 2.26e+10 47874483
```

```
## 
## Call:
## lm(formula = Intensity ~ PixelFactor, data = IntensityData.Reg)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
##  -9385  -4875  -1511   1506  16249 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    35654.6     1997.4   17.85   <2e-16 ***
## PixelFactor2     232.5     2824.7    0.08     0.93    
## PixelFactor3     369.3     2824.7    0.13     0.90    
## PixelFactor4     762.3     2824.7    0.27     0.79    
## PixelFactor5    1187.3     2824.7    0.42     0.67    
## PixelFactor6    1577.0     2824.7    0.56     0.58    
## PixelFactor7    2144.3     2824.7    0.76     0.45    
## PixelFactor8    2463.8     2824.7    0.87     0.38    
## PixelFactor9    2591.3     2824.7    0.92     0.36    
## PixelFactor10   2608.5     2824.7    0.92     0.36    
## PixelFactor11   2719.7     2824.7    0.96     0.34    
## PixelFactor12   2805.9     2824.7    0.99     0.32    
## PixelFactor13   2617.2     2824.7    0.93     0.35    
## PixelFactor14   2278.3     2824.7    0.81     0.42    
## PixelFactor15   2384.5     2824.7    0.84     0.40    
## PixelFactor16   2592.1     2824.7    0.92     0.36    
## PixelFactor17   2700.7     2824.7    0.96     0.34    
## PixelFactor18   2638.7     2824.7    0.93     0.35    
## PixelFactor19   2948.0     2824.7    1.04     0.30    
## PixelFactor20   2855.7     2824.7    1.01     0.31    
## PixelFactor21   2879.8     2824.7    1.02     0.31    
## PixelFactor22   2841.3     2824.7    1.01     0.32    
## PixelFactor23   2740.1     2824.7    0.97     0.33    
## PixelFactor24   2538.8     2824.7    0.90     0.37    
## PixelFactor25   2629.2     2824.7    0.93     0.35    
## PixelFactor26   2747.0     2824.7    0.97     0.33    
## PixelFactor27   2926.1     2824.7    1.04     0.30    
## PixelFactor28   2914.6     2824.7    1.03     0.30    
## PixelFactor29   2930.5     2824.7    1.04     0.30    
## PixelFactor30   2989.3     2824.7    1.06     0.29    
## PixelFactor31   2839.3     2824.7    1.01     0.32    
## PixelFactor32   2516.9     2824.7    0.89     0.37    
## PixelFactor33   2198.4     2824.7    0.78     0.44    
## PixelFactor34   2019.8     2824.7    0.72     0.47    
## PixelFactor35   1661.7     2824.7    0.59     0.56    
## PixelFactor36   1135.0     2824.7    0.40     0.69    
## PixelFactor37    398.0     2824.7    0.14     0.89    
## PixelFactor38    336.4     2824.7    0.12     0.91    
## PixelFactor39    437.5     2824.7    0.15     0.88    
## PixelFactor40   -115.1     2824.7   -0.04     0.97    
## PixelFactor41   -690.4     2824.7   -0.24     0.81    
## PixelFactor42    -17.1     2824.7   -0.01     1.00    
## PixelFactor43    627.9     2824.7    0.22     0.82    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 6920 on 473 degrees of freedom
## Multiple R-squared: 0.0266,	Adjusted R-squared: -0.0599 
## F-statistic: 0.307 on 42 and 473 DF,  p-value: 1
```


#### Model 3: Intensity ~ Genome location (as factor) + GroupID

```
## Analysis of Variance Table
## 
## Response: Intensity
##              Df   Sum Sq  Mean Sq F value Pr(>F)    
## GroupID      11 2.12e+10 1.93e+09   629.6 <2e-16 ***
## PixelFactor  42 6.18e+08 1.47e+07     4.8 <2e-16 ***
## Residuals   462 1.42e+09 3.07e+06                   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## Analysis of Variance Table
## 
## Model 1: Intensity ~ PixelFactor
## Model 2: Intensity ~ GroupID + PixelFactor
##   Res.Df      RSS Df Sum of Sq   F Pr(>F)    
## 1    473 2.26e+10                            
## 2    462 1.42e+09 11  2.12e+10 630 <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```
## 
## Call:
## lm(formula = Intensity ~ GroupID + PixelFactor, data = IntensityData.Reg)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
##  -5816   -892     64   1150   3970 
## 
## Coefficients:
##                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)     49293.5      566.4   87.03  < 2e-16 ***
## GroupID2333529  -1251.9      377.6   -3.32  0.00099 ***
## GroupID2374503 -15703.5      377.6  -41.59  < 2e-16 ***
## GroupID2375652 -14323.3      377.6  -37.93  < 2e-16 ***
## GroupID2387634 -20832.5      377.6  -55.17  < 2e-16 ***
## GroupID2391723 -10667.8      377.6  -28.25  < 2e-16 ***
## GroupID2392072 -19948.4      377.6  -52.83  < 2e-16 ***
## GroupID2392327 -13931.8      377.6  -36.90  < 2e-16 ***
## GroupID2392825 -15492.6      377.6  -41.03  < 2e-16 ***
## GroupID2395072 -14871.4      377.6  -39.39  < 2e-16 ***
## GroupID2396530 -17996.3      377.6  -47.66  < 2e-16 ***
## GroupID2401144 -18647.6      377.6  -49.39  < 2e-16 ***
## PixelFactor2      232.5      714.7    0.33  0.74511    
## PixelFactor3      369.2      714.7    0.52  0.60567    
## PixelFactor4      762.3      714.7    1.07  0.28672    
## PixelFactor5     1187.2      714.7    1.66  0.09737 .  
## PixelFactor6     1577.0      714.7    2.21  0.02785 *  
## PixelFactor7     2144.3      714.7    3.00  0.00284 ** 
## PixelFactor8     2463.8      714.7    3.45  0.00062 ***
## PixelFactor9     2591.3      714.7    3.63  0.00032 ***
## PixelFactor10    2608.5      714.7    3.65  0.00029 ***
## PixelFactor11    2719.7      714.7    3.81  0.00016 ***
## PixelFactor12    2805.9      714.7    3.93  1.0e-04 ***
## PixelFactor13    2617.2      714.7    3.66  0.00028 ***
## PixelFactor14    2278.2      714.7    3.19  0.00153 ** 
## PixelFactor15    2384.5      714.7    3.34  0.00092 ***
## PixelFactor16    2592.1      714.7    3.63  0.00032 ***
## PixelFactor17    2700.7      714.7    3.78  0.00018 ***
## PixelFactor18    2638.7      714.7    3.69  0.00025 ***
## PixelFactor19    2948.0      714.7    4.12  4.4e-05 ***
## PixelFactor20    2855.7      714.7    4.00  7.5e-05 ***
## PixelFactor21    2879.8      714.7    4.03  6.5e-05 ***
## PixelFactor22    2841.2      714.7    3.98  8.2e-05 ***
## PixelFactor23    2740.1      714.7    3.83  0.00014 ***
## PixelFactor24    2538.7      714.7    3.55  0.00042 ***
## PixelFactor25    2629.2      714.7    3.68  0.00026 ***
## PixelFactor26    2747.0      714.7    3.84  0.00014 ***
## PixelFactor27    2926.1      714.7    4.09  5.0e-05 ***
## PixelFactor28    2914.6      714.7    4.08  5.4e-05 ***
## PixelFactor29    2930.5      714.7    4.10  4.9e-05 ***
## PixelFactor30    2989.2      714.7    4.18  3.5e-05 ***
## PixelFactor31    2839.2      714.7    3.97  8.3e-05 ***
## PixelFactor32    2516.9      714.7    3.52  0.00047 ***
## PixelFactor33    2198.4      714.7    3.08  0.00222 ** 
## PixelFactor34    2019.7      714.7    2.83  0.00492 ** 
## PixelFactor35    1661.7      714.7    2.32  0.02051 *  
## PixelFactor36    1135.0      714.7    1.59  0.11297    
## PixelFactor37     398.0      714.7    0.56  0.57790    
## PixelFactor38     336.4      714.7    0.47  0.63811    
## PixelFactor39     437.5      714.7    0.61  0.54073    
## PixelFactor40    -115.1      714.7   -0.16  0.87215    
## PixelFactor41    -690.4      714.7   -0.97  0.33460    
## PixelFactor42     -17.1      714.7   -0.02  0.98094    
## PixelFactor43     627.9      714.7    0.88  0.38012    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 1750 on 462 degrees of freedom
## Multiple R-squared: 0.939,	Adjusted R-squared: 0.932 
## F-statistic:  134 on 53 and 462 DF,  p-value: <2e-16
```


#### Model 4: Intensity ~ GroupID + G + C + T

```
## 
## Call:
## lm(formula = Intensity ~ GroupID + GG + CC + TT, data = IntensityData.Reg)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
##  -8157  -1052     69   1225   4828 
## 
## Coefficients:
##                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       45316       1333   34.00  < 2e-16 ***
## GroupID2333529    -1252        426   -2.94  0.00341 ** 
## GroupID2374503   -15704        426  -36.91  < 2e-16 ***
## GroupID2375652   -14323        426  -33.67  < 2e-16 ***
## GroupID2387634   -20832        426  -48.97  < 2e-16 ***
## GroupID2391723   -10668        426  -25.07  < 2e-16 ***
## GroupID2392072   -19948        426  -46.89  < 2e-16 ***
## GroupID2392327   -13932        426  -32.75  < 2e-16 ***
## GroupID2392825   -15493        426  -36.41  < 2e-16 ***
## GroupID2395072   -14871        426  -34.95  < 2e-16 ***
## GroupID2396530   -17996        426  -42.30  < 2e-16 ***
## GroupID2401144   -18648        426  -43.83  < 2e-16 ***
## GG                 7522       1822    4.13  4.3e-05 ***
## CC                 5473       1441    3.80  0.00016 ***
## TT                10935       2813    3.89  0.00012 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 1970 on 501 degrees of freedom
## Multiple R-squared: 0.916,	Adjusted R-squared: 0.914 
## F-statistic:  391 on 14 and 501 DF,  p-value: <2e-16
```

```
## Analysis of Variance Table
## 
## Model 1: Intensity ~ GroupID
## Model 2: Intensity ~ GroupID + GG + CC + TT
##   Res.Df      RSS Df Sum of Sq    F  Pr(>F)    
## 1    504 2.03e+09                              
## 2    501 1.95e+09  3  84485929 7.24 9.2e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


#### Model 5: Intensity ~ GroupID + CpG (Yes/No)

```
## 
## Call:
## lm(formula = Intensity ~ GroupID + CpG, data = IntensityData.Reg)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
##  -8016  -1096     34   1285   4716 
## 
## Coefficients:
##                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)       51072        306  166.99  < 2e-16 ***
## GroupID2333529    -1252        429   -2.92  0.00368 ** 
## GroupID2374503   -15704        429  -36.60  < 2e-16 ***
## GroupID2375652   -14323        429  -33.38  < 2e-16 ***
## GroupID2387634   -20832        429  -48.56  < 2e-16 ***
## GroupID2391723   -10668        429  -24.86  < 2e-16 ***
## GroupID2392072   -19948        429  -46.49  < 2e-16 ***
## GroupID2392327   -13932        429  -32.47  < 2e-16 ***
## GroupID2392825   -15493        429  -36.11  < 2e-16 ***
## GroupID2395072   -14871        429  -34.66  < 2e-16 ***
## GroupID2396530   -17996        429  -41.95  < 2e-16 ***
## GroupID2401144   -18648        429  -43.46  < 2e-16 ***
## CpG1                786        237    3.31  0.00098 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Residual standard error: 1990 on 503 degrees of freedom
## Multiple R-squared: 0.914,	Adjusted R-squared: 0.912 
## F-statistic:  448 on 12 and 503 DF,  p-value: <2e-16
```

```
## Analysis of Variance Table
## 
## Model 1: Intensity ~ GroupID
## Model 2: Intensity ~ GroupID + CpG
##   Res.Df      RSS Df Sum of Sq  F  Pr(>F)    
## 1    504 2.03e+09                            
## 2    503 1.99e+09  1  43479387 11 0.00098 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


### Comments
#### - Can't distinguish in behavior between CpG Islands chosen by Adi and those in the Alpha Globin sites
#### - Need a more thorough analysis, but need consensus on direction of project

### Proposed Next Steps
#### - Need to align (by pixel) the molecules accurately. Identify and discard outliers
#### - Need to conduct more thorough analysis, involving multiple groups
#### - Shelve CpG island analysis for now
#### - Need sample selection plan - how to group similar fragments in the same pool to observe the effect of sequence composition on intensity


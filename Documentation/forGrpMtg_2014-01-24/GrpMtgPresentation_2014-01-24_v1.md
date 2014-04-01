Fluorescence Intensity & Sequence Composition
=============================================
### Author: Subhrangshu Nandi
### Date: 2014-01-24

# Fundamental hypothesis
## - "Quantum yield and fluorescence lifetime for YO (oxazole yellow dye) complexed with GC is twice as large as that complexed with AT." (Larsson, et. al, 1994)
## - Some relationship between observed Fluorescence Intensity & Genomic Sequence Composition

<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />

# Motivation
## Motivation 1: Construct signature intensity (relative) patterns of as many intervals as possible, with an acceptable range of intensity values
![Image](figure/DNA-sequence_large_verge_medium_landscape.jpg)




![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 


## Motivation 2: Given any relative intensity profile, identify (with certain confidence level) possible genomic location(s) that could produce

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3.png) 

<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />

## Progress
## Step 1: Relationship between high GC areas and intensity
### Data: MM52, with all 7,134 groups
### - nMaps aligned to the reference
### - average 50 molecules (nMaps) per group
<br />
<br />
### - controlled for different surfaces and channels
<br />
### - used random effects regression technique to quantify relationship between fragment level GC content & total intensity
<br />
### - Established statistically significant result <br />
## (Higher GC content implies higher intensity)

<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />

## Step 2: How is the relative intensity affected by sequence composition (beyond just GC content)
### Same Data source (finer resolution)
<br />
### - 1 pixel ~ 200 bp subintervals
<br />
### - Consider sequence composition (beyond just % of GC)
<br />
### - Observe relative intensity profiles of as many intervals as possible
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4.png) 

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5.png) 


![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 

<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />

## Work in progress
### - Align intervals from multiple nMaps dynamically
<br />
### - Generate a consensus signal for as many intervals of the genome
<br />
### - Establish a methodology to compare new intensity signal with consensus signal

<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />
<br />

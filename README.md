# WDBC

This repository provides the R code used in the real data study of paper "Interval Estimation for the Youden Index of a Continuous Diagnostic Test with Verification Biased Data" by Shirui Wang, Shuangfei Shi & Gengsheng Qin (2024).

## About WDBC dataset

Wisconsin Diagnostic Breast Cancer (WDBC) dataset was created by Dr. William H. Wolberg, W. Nick Street and Olvi L. Mangasarian from University of Wisconsin in 1995. In this dataset, features were computed from a digitized image of a fine needle aspirate (FNA) of breast masses and the verified diagnosis (malignant or benign) of breast masses were recorded. They describe characteristics of the cell nuclei present in the image. There are 569 instances of breast masses (357 benign, 212 malignant), each has 32 attributes, including their IDs, true diagnoses and 30 real-valued input features. There is no missing value in the dataset. This database is also available on UCI Machine Learning Repository: https://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+%28Diagnostic%29.

Attribute Information:

Column 1: ID number

Column 2: Diagnosis (M = malignant, B = benign)

Column 3-32: Ten real-valued features computed for each cell nucleus, including:

a) radius (mean of distances from center to points on the perimeter)

b) texture (standard deviation of gray-scale values)

c) perimeter

d) area

e) smoothness (local variation in radius lengths)

f) compactness (perimeter^2 / area - 1.0)

g) concavity (severity of concave portions of the contour)

h) concave points (number of concave portions of the contour)

i) symmetry

j) fractal dimension ("coastline approximation" - 1)

The mean, standard error and "worst" or largest (mean of the three largest values) of these features were computed for each image, resulting in 30 features. For instance, field 3 is Mean Radius, field 13 is Radius SE, field 23 is Worst Radius.

All feature values are recoded with four significant digits.

**wdbc.new** is a subset of the WDBC dataset that we create to be used in our real data study, with 344 benign subjects and 86 malignant subjects, so the disease prevalence is 0.2. 

## About wdbc.ci file

**wdbc.ci** is the .R file that contains the R code for our study.

1. To run the code and replicate our results in the paper, first install and load R packages rocbc, ThresholdROC and dplyr.
2. Load dataset 'wdbc.new'.
3. **jc** contains the true J and cutoff point for all biomarkers, computed using complete data.
4. **wdbc.ci** is the main function to compute the CIs. It requires 5 arguments:
   * *seed*: The seed used to generate random numbers, use 123 (default) to replicate our results
   * *k*: The parameter used to control the missing proportion, 1 by default
   * *bio*: Biomarkers considered in the study, in the paper we use 4, 25, 29
   * *a*: Indicates the bias-corrected estimator (FI, MSI, IPW, SPE) used to compute the CIs, 'SPE' by default
   * *misv*: Indicates if the verification model is misspecified, 'FALSE' by default  
5. The function returns a list that contains the True J, verification proportion for both diseased and healthy groups, CIs computed from packages rocbc and ThresholdROC and CIs computed using our proposed methods.

################################################
# SAMPLE CODE FOR HOUSEMAN, ACCOMANDO, ET AL.
#
# November 6, 2011
#
################################################

##############################################
##### ##### Step 1: Fit Validation Model (S0)
##############################################

# Load test data
load("Example-WBC-Data.RData")

# Define validation model
theModel = y~PanTcell+CD8T+CD4T+NK+Bcell+Mono+Gran
sizeModel =              8 #Number of coefficients in "theModel"

M = 500 # Number of CpGs on array (test data is a subset of 27K)

# Linear transformation of coefficient vector
#  representing contrast to test F statistic
L.forFstat = diag(sizeModel)[-1,]  #All non-intercept coefficients

# Initialize various containers
sigmaResid = sigmaIcept = nObserved = nClusters = Fstat = rep(NA, M)
coefEsts = matrix(NA, M, sizeModel)
coefVcovs =list()

library(nlme)
for(j in 1:M){ # For each CpG

  #Remove missing methylation values
  ii = !is.na(validationData_Assay[j,])
  nObserved[j] = sum(ii)
  validationData_Pheno$y = validationData_Assay[j,]

  if(j%%10==0) cat(j,"\n") # Report progress

  try({ # Try to fit a mixed model to adjust for plate
    fit = try(lme(theModel, random=~1|PLATE, data=validationData_Pheno[ii,]))

    if(inherits(fit,"try-error")){ # If LME can't be fit, just use OLS
      fit = lm(theModel, data=validationData_Pheno[ii,])
      fitCoef = fit$coef
      sigmaResid[j] = summary(fit)$sigma
      sigmaIcept[j] = 0
      nClusters[j] = 0
    }
    else{ 
      fitCoef = fit$coef$fixed
      sigmaResid[j] = fit$sigma
      sigmaIcept[j] = sqrt(getVarCov(fit)[1])
      nClusters[j] = length(fit$coef$random[[1]])
    }
    coefEsts[j,] = fitCoef
    coefVcovs[[j]] = vcov(fit)

    useCoef = L.forFstat %*% fitCoef
    useV = L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
    Fstat[j] = (t(useCoef) %*% solve(useV, useCoef))/sizeModel
  })
}

# Name the rows so that they can be easily matched to the target data set
rownames(coefEsts) = rownames(validationData_Assay)
colnames(coefEsts) = names(fitCoef)

# Get P values corresponding to F statistics
Pval = 1-pf(Fstat, sizeModel, nObserved - nClusters - sizeModel + 1)

# NOTE:  test data consists of CpGs having 500 largest F statistics
#  in descending order of magnitude
# For the full 27K array, it is necessary to sort by Pvalue and/or Fstatistic
#  and choose the top CpGs
table(sign(diff(Fstat))) # Fstats decrease as j increases from 1 to M


# Save for future use
save(file="Validation-Coefficients.RData", compress=TRUE,
  list=c("coefEsts","coefVcovs","sigmaIcept","sigmaResid" ,"L.forFstat","Fstat",
    "nClusters", "nObserved")
)

##############################################
##### ##### Step 2: Fit Target Model (S1)
##############################################

# Load test data [*Footnote 1*]
load("Example-WBC-Data.RData")
load("Validation-Coefficients.RData")

library(nlme)
source("./rscript/wbcInference.R")

# Contrast matrix for PanT rearrangement [*Footnote 2*]
Lwbc = diag(8)[-(1:2),] 
Lwbc[1:2,2] = 1 ; Lwbc[,1] = 1
rownames(Lwbc) = colnames(coefEsts)[-(1:2)]
colnames(Lwbc) = colnames(coefEsts)

Lwbc # View contrast matrix

# Denominator degrees-of-freedom for parametric bootstrap
degFree = nObserved - nClusters - 7

CpGSelection = rownames(coefEsts)[1:100] # Use the top 100
#Note:  if the CpGs were scattered throughout the array,
#    you would want to select them by name as is indicated here.
#    For this sample version, it would be easier just to use
#    "[1:100]"
targetEst = inferWBCbyLme(
  targetDataHNSCC_Assay[1:100,],     # Target methylation (cpGs x subjects)
  targetDataHNSCC_Covariates,        # Target phenotype frame (subjects x covariates)
  y~case+gender+ageCtr,  # Target model (fixed effects)
  ~1|BeadChip,           # Target adjustment (random effects) [*Footnote 3*]
  coefEsts[CpGSelection,],      # Raw coefficient estimates for WBC 
  Lwbc                   # Contrast matrix [*Footnote 2*]
)

targetEst # View model estimates

### Get bootstraps [*Footnote 3*]
# Warning:  this can take a long time
targetBoot = bootInferWBCbyLme( 
  targetDataHNSCC_Assay[1:100,],     # Target methylation (cpGs x subjects)
  targetDataHNSCC_Covariates,        # Target phenotype frame (subjects x covariates)
  y~case+gender+ageCtr,     # Target model (fixed effects) [*Footnote 3*]
  ~1|BeadChip,              # Target adjustment (random effects) [*Footnote 4*] 
  coefEsts[CpGSelection,],           # Raw coefficient estimates for WBC 
  Lwbc,                     # Contrast matrix [*Footnote 2*]
  R=250,                    # Number of bootstrap samples to run
  vcovWBC=coefVcovs[1:100],   # WBC fixed effects v-cov estimates [*Footnote 5*]
  degFree=degFree[1:100]      # WBC degrees-of-freedom [*Footnote 6*]
)

targetBoot # View bootstrap summary

# View summary with bootstraps
summary(targetEst, targetBoot)

##############################################
##### ##### Step 3: View projections
##############################################

# Load test data [*Footnote 1*]
load("Example-WBC-Data.RData")
load("Validation-Coefficients.RData")

source("wbcInference.R")

# Contrast matrix for PanT rearrangement [*see details in help*]
Lwbc = diag(8)[-(1:2),] 
Lwbc[1:2,2] = 1 ; Lwbc[,1] = 1
rownames(Lwbc) = colnames(coefEsts)[-(1:2)]
colnames(Lwbc) = colnames(coefEsts)

Lwbc # View contrast matrix

CpGSelection = rownames(coefEsts)[1:100] # Use the top 100
#Note:  if the CpGs were scattered throughout the array,
#    you would want to select them by name as is indicated here.
#    For this sample version, it would be easier just to use
#    "[1:100]"

####### Projections for HNSCC data

unconstrainedCoefs = projectWBC(
  targetDataHNSCC_Assay[1:100,],
  coefEsts[CpGSelection,],    
  Lwbc, nonnegative = FALSE)

constrainedCoefs = projectWBC(
  targetDataHNSCC_Assay[1:100,],
  coefEsts[CpGSelection,],    
  Lwbc)

head(unconstrainedCoefs)
head(constrainedCoefs)

####### Projections for mixture experiment data

ObsMix = projectWBC(
  mixtureExperiment_Assay[1:100,],
  coefEsts[CpGSelection,],    
  Lwbc)

ExMix = matrix(0, 12, 5)
colnames(ExMix) = c("Pan-T", colnames(ObsMix)[-(1:2)])
for(i in 1:12){
  ExMix[i,"Bcell"] = mixtureExperiment_Design$B[i]
  ExMix[i,"Pan-T"] = mixtureExperiment_Design$T[i]
  ExMix[i,"Gran"] = mixtureExperiment_Design$Gran[i]
  ExMix[i,"Mono"] = mixtureExperiment_Design$Mono[i]
}

colnames(ObsMix) = c("T (CD8+)", "T (CD4+)", "NK", "B Cell", "Monocyte", "Granulocyte")
colnames(ExMix) = c("T Cell", "NK", "B Cell", "Monocyte", "Granulocyte")
rownames(ObsMix) = rownames(ExMix) = mixtureExperiment_Design$strip

colSortO = c(5,4,3,1,2,6)[6:1]
colSortE = c(4,3,2,1,5)[5:1]
rowSort = c(2*(1:6)-1, 2*(1:6))

# Print Observed and Expected
round(100*ObsMix[rowSort, colSortO],1)
round(100*ExMix[rowSort, colSortE],1)

# Print errors
ObsMixCollapse = ObsMix[rowSort, colSortO]
ObsMixCollapse = cbind(ObsMixCollapse[,1], 
   ObsMixCollapse[,2]+ObsMixCollapse[,3], ObsMixCollapse[,4:6])
colnames(ObsMixCollapse) = colnames(ExMix)
round(100*(ObsMixCollapse-ExMix[rowSort, colSortE]),1)

############################################

##### ##### FOOTNOTES
#
#[*Footnote 1*] 
#
# Test data consist of 500 CpGs chosen from Illumina Infinium 27K array
# These were the 500 most informative CpGs for distinguishing WBC type
# (obtained by ordering the ANOVA F-statistics from largest to smallest)
# These CpGs are already ordered from "most informative" to "least informative".
# In reconstructing a new validation set, it would be necessary to choose
# the CpGs in the same manner (ordering by F statistic) or using a more careful procedure.
#
#
#[*Footnote 2*]
#
# The hierarchical design reflected in the WBC phenotype file
# requires that an appropriate contrast matrix be defined to extract desired 
# WBC types.  
#
#
#[*Footnote 3*]
#
# Proper syntax for specifying a fixed effects model is as usual, except that
# the outcome must be specified as "y".  Consequently, "y" is a reserved name
# in this function (i.e. any variable in the phenotype file named "y" will
# not be accessible).  
#
#
#[*Footnote 4*]
#
# In this WBC function, adjustment for chip effects is *required*
# Another function implements a faster, OLS-based algorithm that ignores
# chip effects ("inferWBCbyLm").  The proper syntax is "~1|BeadChip",
# where "BeadChip" is the name of the variable in the phenotype file that
# specifies chip number.
#
#
#[*Footnote 5*]
#
# The argument "vcovWBC" may be omitted if no double-bootstrap estimates
# are desired.  However, the double-bootstrap estimates are not particularly
# processor-intensive, so there is little additional cost to computing them.
# They do require WBC fixed effects variance-covariance estimates.
# Make sure to select the sub-list that corresponds to the CpGs selected for
# the computation!
#
#
#[*Footnote 6*]
#
# The argument "degFree" may be omitted if Gaussian noise is to be used
# in the double-bootstrap, instead of t-distributed noise.  T-distributed noise
# is likely to be (slightly) more accurate.
#

# rspmssimulation

This repository contains five files. Before describing each file, I will first explain the notation used.

**Notation**

The following notation will be used in the comments (####) of each file.

[_A_][_B_][_C_][_D_][_E_]

where

- _A_ indicates whether the surrogate endpoint dataset is Normally distributed. It has two options: **Normal** (for Normally distributed), and **XNormal** (for non-Normally distributed).
- _B_ describes the ranges of the datasets for the control and treatment arms of the surrogate endpoint. It has two options: **Overlap** (the two ranges mostly overlap), and **Subset** (the range of the control arm data is mostly contained within the range of the treatment arm data).
- _C_ indicates the average sample size between the control and treatment arms. It has two options: **200**, and **1000**.
- _D_ indicates the 
- _E_ indicates the allocation ratio. It has two options: **1**, and **3**.

To illustrate this notation, the setting **XNormalSubset1000HM3** would indicate the following:

- the surrogate endpoint is non-Normally distributed;
- the range of the control arm of the surrogate endpoint data is mostly contained within the range of the treatment arm of the surrogate endpoint data;
- the average sample size between the two arms is 1000;
- the data-generating mechanism is HM,, and
- the allocation ratio is r=3,

The third and fifth items mean that the sample size of the treatment arm of the surrogate endpoint is 1500, and the sample size of the control arm of the surrogate endpoint is 500.

**Files**

The first file (rspmsillustration.R) is 

- **rspmsillustration.R**
  
As this is an illustration, extensive comments (#) to the code are provided. This code runs simulations in the case **NormalSubset200LL1**; that is, for the data-generating mechanism where both RM and RS (as defined in Proposition 2.3) are low (hence denoted **LL**), in the case of a Normal surrogate endpoint such that the range of surrogate values in the control arm lies within the range of surrogate values in the treatment arm (subsetted supports), under equal allocation (r=1) and a sample size of 200 per arm (total sample size 400).

The results correspond to the first two rows of Table 2.8 (point estimates of performance measures) and the first two rows of Table B.4 (Monte Carlo standard errors of performance measures; Appendix B) in my thesis.

- **NormalSubsetted.R**

This code contains all combinations of **NormalSubset[_C_][_D_][_E_]**. The code for each combination is analogous to rspmsillustration.R, so the comments are mostly removed.

- **NormalOverlapping.R**

This code contains all combinations of **NormalOverlap[_C_][_D_][_E_]**. The code for each combination is analogous to rspmsillustration.R, so the comments are mostly removed.

- **NonNormalSubsetted.R**

This code contains all combinations of **XNormalSubset[_C_][_D_][_E_]**. The code for each combination is analogous to rspmsillustration.R, so the comments are mostly removed.

- **NonNormalOverlapping.R**

This code contains all combinations of **XNormalOverlap[_C_][_D_][_E_]**. The code for each combination is analogous to rspmsillustration.R, so the comments are mostly removed.

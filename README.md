# MultiwayClassification

This is an R package to perform linear classification for data with multi-way structure.  The distance-weighted discrimination (DWD) or support vector machine (SVM) classification objectives are optimized under the assumption that the multi-way coefficients have low rank [1].  Additional functions perform multiway DWD with sparsity [2]. 
This package depends on the packages `DWD` (for DWD), `kernlab` (for SVM), and `sdwd` (for sparse DWD). `DWD` is not currently available on CRAN, and so will need to be installed via its url:
```
install.packages("https://cran.r-project.org/src/contrib/Archive/DWD/DWD_0.11.tar.gz",repos = NULL, type = "source")
```
The `MultiwayClassification` package can then be installed, directly from GitHub, using the devtools library:

```
install.packages('devtools')
library(devtools)
install_github("lockEF/MultiwayClassification")
``` 

Code for this package was written primarily by Tianmeng Lyu (for multiway DWD and SVM) and Bin Guo (for multiway sparse DWD).     

[1] Lyu, T., Lock, E.F., & Eberly, L. E. (2017). Discriminating sample groups with multi-way data. Biostatistics, 18 (3): 434–450. https://arxiv.org/abs/1606.08046 .

[2] Guo, B., Eberly, L.E., Henry, P.G., Lenglet, C. & Lock, E. F. (2020). Sparse multiway distance weighted discrimination. Preprint.

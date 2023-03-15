# ComBat and CovBat Families of Harmonization Methods

**Maintainer**: Andrew Chen, andrewac@pennmedicine.upenn.edu

**License**: Artistic License 2.0

The ComBat Family extends the original ComBat methodology to enable flexible
covariate modeling, leveraging efficient R implementations of regression models.
A method that belongs in the ComBat Family satisfies the following conditions:

1. Modeling of covariate effects in location and scale
2. Batch effects in location and scale of measurements
3. Empirical Bayes step for borrowing information across features

The CovBat Family likewise extends the original CovBat methodology to enable flexible
covariate modeling while removing batch effects in the mean and covariance of measurements.

**NOTE:** This package is still a work-in-progress and will be updated soon to include the following features:
- Nonparametric Empirical Bayes step
- Empirical Bayes step leveraging non-Gaussian data distributions
- Functions for evaluating batch effects before and after harmonization

## 1. Installation
The R package can be installed via devtools by running the following code

```
# install.packages("devtools")
devtools::install_github("andy1764/ComBatFamily", build_vignettes = TRUE)
```

Then, you can load this package via

```
library(ComBatFamily)
```

## 2. Usage
Vignettes are provided for both the ComBat family `comfam` and the CovBat family `covfam` which can be accessed via

```
vignette("comfam")
vignette("covfam")
```

Note that non-Gaussian data distributions are supported by functions such as `glm` and `gamlss`; however, the batch effect correction may produce harmonized data outside the original range of values. For now, specification of non-Gaussian distributions will generate a warning. This support is still a work-in-progress.

## 3. Citations
The original ComBat methodology is implemented in R, Matlab, and Python at https://github.com/Jfortin1/ComBatHarmonization. When using ComBat, please cite the following papers:

> Fortin, J.-P., Cullen, N., Sheline, Y. I., Taylor, W. D., Aselcioglu, I., Cook, P. A., Adams, P., Cooper, C., Fava, M., McGrath, P. J., McInnis, M., Phillips, M. L., Trivedi, M. H., Weissman, M. M., & Shinohara, R. T. (2018). Harmonization of cortical thickness measurements across scanners and sites. *NeuroImage*, 167, 104–120. https://doi.org/10.1016/j.neuroimage.2017.11.024
> 
> Fortin, J.-P., Parker, D., Tunç, B., Watanabe, T., Elliott, M. A., Ruparel, K., Roalf, D. R., Satterthwaite, T. D., Gur, R. C., Gur, R. E., Schultz, R. T., Verma, R., & Shinohara, R. T. (2017). Harmonization of multi-site diffusion tensor imaging data. *NeuroImage*, 161, 149–170. https://doi.org/10.1016/j.neuroimage.2017.08.047
> 
> Johnson, W. E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*, 8(1), 118–127. https://doi.org/10.1093/biostatistics/kxj037

The original CovBat method is available at https://github.com/andy1764/CovBat_Harmonization. If implemented, please cite the original article:

> Chen, A. A., Beer, J. C., Tustison, N. J., Cook, P. A., Shinohara, R. T., Shou, H., & Initiative, T. A. D. N. (2022). Mitigating site effects in covariance for machine learning in neuroimaging data. *Human Brain Mapping*, 43(4), 1179–1195. https://doi.org/10.1002/hbm.25688

For longitudinal ComBat, the original R package is available at https://github.com/jcbeer/longCombat with corresponding paper:

> Beer, J. C., Tustison, N. J., Cook, P. A., Davatzikos, C., Sheline, Y. I., Shinohara, R. T., & Linn, K. A. (2020). Longitudinal ComBat: A method for harmonizing longitudinal multi-scanner imaging data. *NeuroImage*, 220, 117129. https://doi.org/10.1016/j.neuroimage.2020.117129

For ComBat-GAM, the Python implementation is available via https://github.com/rpomponio/neuroHarmonize with corresponding paper:

> Pomponio, R., Shou, H., Davatzikos, C., et al., (2020). "Harmonization of large MRI datasets for the analysis of brain imaging patterns throughout the lifespan." *NeuroImage* 208. https://doi.org/10.1016/j.neuroimage.2019.116450.



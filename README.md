# ComBat Family of Harmonization Methods

WARNING: This is the developmental branch of ComBatFamily, new functions and modified functions are still in testing.

**Author**: Andrew A. Chen, chenandr@musc.edu 

**Maintainers**:  
Andrew A. Chen, chenandr@musc.edu  
Haochang Shou, hshou@pennmedicine.upenn.edu  
Margaret Gardner, margaret.gardner@pennmedicine.upenn.edu  
Zheng Ren, zheng.ren@pennmedicine.upenn.edu  
Randa Melham, randa.melhem@pennmedicine.upenn.edu  

**License**: Artistic License 2.0

The ComBat Family extends the original ComBat methodology to enable flexible covariate modeling, leveraging efficient R implementations of regression models. A method that belongs in the ComBat Family satisfies the following conditions:

- Modeling of covariate effects in location and scale
- Batch effects in location and scale of measurements
- Empirical Bayes step for borrowing information across features

ComBat Family methods include:

1. ComBat (Johnson et al., 2007)
2. ComBat-GAM (Pomponio et al., 2020)
3. Longitudinal ComBat (Beer et al., 2020)
4. Robust ComBat (Work-in-progress)
5. ComBatLS (Gardner et al., preprint)

This package also includes the CovBat Family, which likewise extends the original CovBat methodology to enable flexible covariate modeling while removing batch effects in the mean and covariance of measurements.

**NOTE:** This package is still a work-in-progress and will be updated soon to include the following features:

- Nonparametric Empirical Bayes step
- Empirical Bayes step leveraging non-Gaussian data distributions
- Functions for evaluating batch effects before and after harmonization

## 1. Installation
The developmental version can be installed via `remotes` by running the following code

```
# install.packages("remotes")
remotes::install_github("andy1764/ComBatFamily@dev", force = TRUE)
```

Then, you can load this package via

```
library(ComBatFamily)
```

## 2. Usage
Vignettes are provided for both the ComBat family `comfam` and the CovBat family `covfam`. To install with vignettes, first install the suggested dependencies via

```
remotes::install_github("jcbeer/longCombat")
remotes::install_github("andy1764/CovBat_Harmonization/R")
remotes::install_github("jfortin1/neuroCombat_Rpackage")
```

Then install the developmental version with vignettes via

```
remotes::install_github("andy1764/ComBatFamily@dev", force = TRUE, build_vignettes = TRUE)
```

Vignettes can then be accessed through

```
vignette("comfam")
vignette("covfam")
```

Example ComBat Family calls for `iris` data, treating `Species` as batch:

```
# Original ComBat
comfam(iris[,1:2], iris$Species, covar = iris[3:4], lm, y ~ Petal.Length + Petal.Width)

# ComBat-GAM
comfam(iris[,1:2], iris$Species, covar = iris[3:4], gam, y ~ s(Petal.Length) + Petal.Width)

# Alternative shorthand functions
combat(iris[,1:2], iris$Species, covar = iris[3:4], y ~ Petal.Length + Petal.Width)
combat_gam(iris[,1:2], iris$Species, covar = iris[3:4], y ~ s(Petal.Length) + Petal.Width)
```

Note that non-Gaussian data distributions are supported by functions such as `glm` and `gamlss`; however, the batch effect correction may produce harmonized data outside the original range of values. For now, specification of non-Gaussian distributions will generate a warning. This support is still a work-in-progress.

## 3. Additional features
On top of unifying existing harmonization packages, we include additional features in this package.

For out-of-sample harmonization, we provide `predict.comfam` to apply estimated harmonization to a specified sample. This function will estimate new batch adjustment parameters if needed, otherwise it will apply existing estimates. `predict.comfam` has been tested for linear models (`lm`) and generalized additive models (`gam`), but may give errors for other chosen models. Below is an example call:

```
com_out <- comfam(iris[1:75,1:2], iris$Species[1:75])

# out-of-sample with new batch
out_pred <- predict(com_out, iris[76:150,1:2], iris$Species[76:150])

# in-sample
in_pred <- predict(com_out, iris[1:25,1:2], iris$Species[1:25])
max(in_pred$dat.combat - com_out$dat.combat[1:25,])
```

EXPERIMENTAL: We also include the `predict.covfam` function for applying `covfam` fits to new data. Our implementation utilizes the same model fits and PCA fit from the original `covfam` harmonization.

```
cov_out <- covfam(iris[1:75,1:4], iris$Species[1:75])

# out-of-sample with new batch
out_pred <- predict(cov_out, iris[76:150,1:4], iris$Species[76:150])

# in-sample
in_pred <- predict(cov_out, iris[1:25,1:4], iris$Species[1:25])
max(in_pred$dat.covbat - cov_out$dat.covbat[1:25,1:4])
```

We also provide a wrapper to access model fit diagnostic plots, `plot.comfam`. Other additional features are in active development.

## 4. Citations
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

> Pomponio, R., Erus, G., Habes, M., Doshi, J., Srinivasan, D., Mamourian, E., Bashyam, V., Nasrallah, I. M., Satterthwaite, T. D., Fan, Y., Launer, L. J., Masters, C. L., Maruff, P., Zhuo, C., Völzke, H., Johnson, S. C., Fripp, J., Koutsouleris, N., Wolf, D. H., … Shou, H., Davatzikos, C. (2020). Harmonization of large MRI datasets for the analysis of brain imaging patterns throughout the lifespan. *NeuroImage*, 208, 116450. https://doi.org/10.1016/j.neuroimage.2019.116450

For ComBatLS, please refer the corresponding preprint:
> Gardner, M., Shinohara, R. T., Bethlehem, R. A. I., Romero-Garcia, R., Warrier, V., Dorfschmidt, L., Shanmugan, S., Seidlitz, J., Alexander-Bloch, A., & Chen, A. A. (2024). ComBatLS: A location- and scale-preserving method for multi-site image harmonization. bioRxiv, 2024.06.21.599875. https://doi.org/10.1101/2024.06.21.599875

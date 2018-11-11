# EPIC metabolomics
This repository puts together R scripts that contain various different workflows for untargeted metabolomics data in the EPIC study.

`Feature_match.R` matches features from different untargeted metabolomics feature tables given mass and retention time tolerances.

`Intake_correlation.R` is a function that takes untargeted metabolomics data from EPIC subjects and calculates partial correlations between spectral feature intensities and food intakes or lifestyle factors of interest. The output is a data frame consisting of one row per feature. This table can then be used to aid compound identification. Two further functions use this output to generate a correlation matrix or a Manhattan plot.

`Intake_correlation_453_obs.R` is a previous version of the above function that will remain archived.

`Coffee_biomarkers_CS.R` and `Coffee_biomarkers_HCC.R` calculate associations between known biomarkers and factors of interest or HCC status in the EPIC cross-sectional and HCC studies respectively.
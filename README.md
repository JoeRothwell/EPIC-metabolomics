# EPIC metabolomics
This repository puts together R scripts that contain various different workflows for untargeted metabolomics data in the EPIC study.

`Feature_match.R` matches features from different untargeted metabolomics feature tables given mass and retention time tolerances.

`Intake_correlation.R` is a function that takes EPIC study metadata, metabolomics data and food intake data and calculates partial correlations between spectral features and food intakes of interest. An output data frame is generated consisting of one row per feature. This table can then be used to aid compound identification or generate a Manhattan plot.

`Intake_correlation_453_obs.R` is a previous version of the above function that will remain archived.
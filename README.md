# EPIC metabolomics
This repository puts together R scripts that contain various different workflows for untargeted metabolomics data in the EPIC study. Analyses are for alcohol, coffee and acylcarnitine projects.

![](manhattan2.png)

_Manhattan plot showing results of a metabolome wide association study (MWAS) on coffee intake. Each point represents a spectral feature generated by the mass spectrometer after injection of serum samples and separation of components._

#### Description of files

`alcohol_study_untarg.R` runs two workflows for discovery of alcohol-related spectral features from untargeted metabolomics analyses. In the first, CS and HCC datasets are feature matched (with m/z and RT tolerances; `fuzzyjoin` package) and the remaining features subject to partial correlation with alcohol intake. In the second, only features from the CS study associated with alcohol intake are matched with all features from the HCC study, and the matched HCC features then tested for correlations with alcohol intake. Note: functions from `Intake_correlation_new.R` are required and must be run first.

`Feature_matching.R` matches features from different untargeted metabolomics feature tables given mass and retention time tolerances. (For alcohol biomarker study)

`Intake_correlation.R` contains a function that performs a metabolome-wide association study for food intake in the EPIC cross-sectional study. Partial correlations are computed between spectral feature intensities and food intakes or lifestyle factors of interest, adjusting for confounders. The output is a data frame consisting of one row per feature, with other information such as median intensity and number of detections also extracted from the feature table. This output is a basis for compound identification, and two further functions use it to generate a correlation matrix or a Manhattan plot.

`Intake_correlation_HCC.R` does the same as above but for the HCC case-control study.

`Intake_correlation_453_obs.R` is a previous version of the above function that will shortly be deleted.

`Intake_correlation_new.R` is a previous version of the above function that will shortly be deleted.

`Coffee_biomarkers_CS.R` and `Coffee_biomarkers_HCC.R` calculate associations between known coffee intake biomarkers and factors of interest or HCC status in the EPIC cross-sectional and HCC studies respectively.

`CS_descriptives_coffee.R` analyses coffee intake data in EPIC and extracts coffee intake data from the whole dataset.

`Acylcarnitines.R` is an analysis of acylcarnitine plasma intensities extracted from untargeted data. Compounds correlations with food intake are tested and a heatmap plotted. Blood acylcarnitine intensies are also tested for associations with WCRF score.
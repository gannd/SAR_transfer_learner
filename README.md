"Spatial_Transfer_Learning_RF.r"

This code conducts the spatial transfer learning for either pre- or post-Irma (tuning with variable "period"), using one of the three feature sets ("sar", "optical", "mix") (tuning with parameter "option").
Variable "prop" indicates the percentage for stratified random sampling for model training.

Important parameters for tuning:

period <- 1 ## option 1 is pre, option 2 is post
option <- 3 ## option 1 sar, 2 optical, 3 both
prop <- 0.8 # This is the proportion for selecting training samples

"Temporal_Transfer_Learning_RF.r"

This code conducts temporal transfer learning by reading pre-Irma data to predict post-Irma canopy heights. "prop" and "option" are similar to the spatial transfer learning. There is no "period" variable in this model. 
prop <- 0.01 # This is the proportion for selecting training samples
option <- 1 ## option 1 is sar, 2 is optical, 3 is both

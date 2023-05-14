# Multilayer_MS
Multilayer approach to diagnose and classify Multiple Sclerosis phenotypes using graph theory measures

## "HELPER" FILES
R files where some functions are defined:
1. helper_functions.R --> script where are defined some functions used in all R Markdowns files expcept functions for graph measurements
2. GraphMeasurements.R --> script where graph measurement functions are defined.

## MAIN FILES
We describe main files and the order intended to execute them.

1. Data_processing.Rmd. --> Preprocessing steps, filtering weights, correcting values for sex and age. 
2. Harmonization.ipynb --> Harmonization with combat. 
3. Separate_layers.Rmd --> Graph measurements conisdering layers as independent layers, without SVD normalization
4. Single_Layer_Measurements.Rmd --> Same R Markdown as previous one, but in this one we use SVD Normalization to appreciate differences between the two cases
5. Multiple_Layer_Measurements.Rmd --> Construct a multiplex graph and perform multiple layer measurements, includind selecting weights for inter-layer connections

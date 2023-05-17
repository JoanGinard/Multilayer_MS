# Multilayer_MS
Multilayer approach to diagnose and classify Multiple Sclerosis phenotypes using graph theory measures

## "HELPER" FILES
R files where some functions are defined:
1. **helper_functions.R** --> script where are defined some functions used in all R Markdowns files expcept functions for graph measurements
2. **GraphMeasurements.R** --> script where graph measurement functions are defined.

## MAIN FILES
We describe main files and the order intended to execute them.

1. **Data_processing.Rmd** --> Preprocessing steps, filtering weights, correcting values for sex and age. 
2. **Harmonization.ipynb** --> Harmonization with combat. 
3. **clinical_data.Rmd** --> file where data from subjects is summarized: age, sex, edss, dd.
### SEPARATE LAYERS
4. **Separate_layers.Rmd** --> Graph measurements conisdering layers as independent layers, without SVD normalization
5. **Separate_layers_results.ipynb** --> Notebook to apply machine learning models to graph measurements from "Separate_layers.Rmd"
### MULTILAYER NETWORK WITH 3 LAYERS
7. **Single_Layer_Measurements.Rmd** --> Same R Markdown as *Separate_layers.Rmd , but in this one we use SVD Normalization to appreciate differences between the two cases
8. **Multiple_Layer_Measurements.Rmd** --> Construct a multiplex graph and perform multiple layer measurements, includind selecting weights for inter-layer connections
9. **Mutilayer_network_results.ipynb** --> Almost same notebook as *Separate_layers_results.ipynb* but applying models to single and mutlilayer mesaurements of a multilayer network (From "Single_Layer_Measurements.Rmd" and "Multiple_Layer_Measurements.Rmd")
### MULTILAYER NETWORK WITH 2 LAYERS
As it GM layers does not seem to yield a great deal of information we reproduce *Multiple_Layer_Measurements.Rmd* and *Mutilayer_network_results.ipynb* but only considering FA and fMRI layer.

11. **Multiple_2Layer_Measurements.Rmd** --> Same R markdown file as *Multiple_Layer_Measurements.Rmd* but using only 2 layers
12. **Mutilayer_network_2layer_results.ipynb** --> Same notebook as *Mutilayer_network_results.ipynb* but using only 2 layers


## HTML folder
HTML Files corresponding to outputs of R markdowns and jupyter notebooks

# Multilayer_MS
Multilayer approach to diagnose and classify Multiple Sclerosis phenotypes using graph theory measures

![Flow_image](images/Flow_files_3.png)

## 1.	“Helper” files
During the process with R I found myself reusing some pieces of code, so I decided include those functions in two separate files:
1.	*helper_functions.R*
2.	*GraphMeasurements.R*

The 'helper_functions.R' file houses general-purpose functions, such as those used for data loading, Singular Value Decomposition (SVD) correction, conducting statistical tests, and correcting for sex and age. These functions are employed in both 'Data_Processing.Rmd' and 'Graph Measurements', as depicted in the figure
The 'GraphMeasurements.R' file, on the other hand, includes functions specific to graph measurements. 

## 2.	Data processing and Harmonization
In the 'Data_Processing.Rmd' file, we apply several steps to process our data:
1.	Eliminate matrices of participants for whom clinical data is not available.
2.	Apply thresholds to FA and GM matrices. Those thresholds are explained in the file and also will be detailed in final master‘s thesis document.
3.	Apply absolute value to fMRI matrices and remove the diagonals.
4.	Correct all matrices for sex and age.
5.	Save the modified matrices.

In the “Harmonization” file, we take into account that 2 different scanners were used to acquire participants’ data. We utilize the NeuroCombat Python library to apply necessary corrections. The analysis reveals that only FA matrices require such adjustments.

## 3.	Graph measurements
Graph measurements phase is executed through four files
1.	*Separate_Layers.Rmd*
2.	*Single_Layer_Measurements.Rmd*
3.	*Multiple_Layer_Measurements.Rmd*
4.	*Multiple_2Layer_Measurements.Rmd*

In the first two files, we perform measurements for individual layers. The distinction is that the 'Separate Layers' file does not perform the Singular Value Decomposition (SVD) since we are considering three independent networks rather than a multi-layer network. Both files are almost identical and were separated to serve as standalone files, rather than integrating them into one file or workflow.
In the third and fourth files, we apply multi-layer measurements. The third file considers 3 layers while the fourth considers only 2. Again, the two files are very similar, designed this way for the same reasons as the previous two.
The primary outputs of all these code files are CSV files, which contain participants’ clinical data along with the graph measures obtained through the code. It is worth noting that only measurements that “pass” statistical tests are included. We have considered the case for only two layers because from single layer we have seen that GM layer does not give meaningful information.  In addition to these outputs, the code also generates images and summary tables that could be useful in the final document.

## 4.	Machine learning models
This phase is comprised of these files:
1.	*Separate_layers_results.ipynb* -  utilizes data from “Separate_Layers.Rmd”
2.	*Multilayer_network_Results.ipynb* -  utilizes data from “Single_Layer_Measurements.Rmd” and “Multiple_Layer_Measurements.Rmd”
3.	*Multilayer_network_2layer_Results.ipynb* – employs data from “Single_Layer_Measurements.Rmd” and “Multiple_2Layer_Measurements.Rmd”
4.	Feature_importance: *"Separate_Layers._Permutation_Importance_HS_PwMS.ipynb", "2-multilayer._Permutation_Importance_HS_PwMS.ipynb", "3-multilayer._Permutation_Importance_HS_PwMS.ipynb"*

Files 1, 2 and 3 apply machine learning models to their respective datasets, while 4 is dedicated to investigating feature importance in HS vs PwMS in RF and XGBoost (it comprises 3 files) .  Like in previous phase, files 1,2 and 3 share code so they can be used as standalone files.

## 5.	Other files
Here we include files that were developed specifically purpose of generating tables, graphs or data summaries.
1.	*clinical_data.Rmd* - offers a summary of participant details such as age, sex, disease duration (dd), and expanded disability status scale (EDSS). It also generates graphs of data distribution.
2.	 *SVDNormalization_comparison.Rmd*  - produces a graph that compares histogram of weights before and after SVD correction
3.	 *Connectogram_circus_plot.ipynb* – creates a connectogram visualization using a circus plot. 

It is important to note that some code used in the main files is reused in files 2 and 3.


## HTML folder
HTML Files corresponding to outputs of R markdowns and jupyter notebooks

## Images Folder
Images produced in code.

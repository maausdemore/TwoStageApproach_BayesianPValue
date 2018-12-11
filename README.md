# TwoStageApproach_BayesianPValue
Code and data associated with "Two-Stage Approach for the Inference of the Source of High-Dimensional and Complex Chemical Data in Forensic Science" (paper currently available at https://arxiv.org/abs/1804.01049)

### Files in Repository
#### R files 
There are two primary .R files associated with this project, a SCRIPT file and a FUNCTION file. 
##### 1) 2_JASA_Outlier_detector_model_SCRIPT.R
This file contains the code for replicating Section 6. "Worked Example" for FTIR spectra of paint in the aforementioned paper. The code allows for replicating the various simulations (obtaining the distribution of the scores, obtaining the distribution of the test statistic, determining the power of the test, and assigning the random match probability for each source in the set of 166 paint sources). Each section of the file also contains the necessary code for reproducing the plots found in the paper. The file is constructed such that the only lines to be changed are those indicating the file paths at the beginning of each section. A guide is provided within the file to let the user know which lines should be updated with the user's file paths. 

##### 2) 2_JASA_Outlier_detector_model_FUNCTIONS.R
This file contains the necessary functions for running the code provided in 2_JASA_Outlier_detector_model_SCRIPT.R above. No changes should be made to this file. Simply save this file in a location of your choice, and update the lines corresponding to the variable named 'my.functions' in the 2_JASA_Outlier_detector_model_SCRIPT.R file with the path associated with this location. The 2_JASA_Outlier_detector_model_SCRIPT.R file will then source the functions.

#### Data 
The data for this paper was provided by Dr. Cyril Muehlethaler at the Université du Québec à Trois Rivières. The data is a folder with 166 individual .csv files, each associated with a unique household paint. The 166 paints in this study are of similar colour. 

Each file contains the FTIR spectra of seven distinct paint chips from a common source. Each .csv file in this folder consists in eight columns of data. The first column corresponds to the wave numbers over which the measurements were taken. The wave numbers range from 4000 cm^-1 to 550 cm^-1, and are captured by an ~7000-dimensional vector. The remaining columns record the resulting FTIR spectra of the paint chips, measured in transmittance.

This Github repository exists in service of open science, to reproduce the analysis from the following paper: 

Balantic, C. M., & Donovan, T. M. (2019).                                   
Statistical learning mitigation of false positives from template-detected data in automated acoustic wildlife monitoring. Bioacoustics. https://www.tandfonline.com/doi/full/10.1080/09524622.2019.1605309                                                                
                                    
Readers of this paper who are interested in reproducing the analysis should begin by downloading the following items from the Github repository:

* ammls [a folder]
   --> contains 'classifiers221.RDS' [an .RDS file that contains the trained/tested classifiers from this manuscript)
* ammonitor-functions.R [an R file that contains all AMMonitor functions necessary to produce the analysis/results in the paper]
* false-positives-script.R [an R script to produce data, models, analysis, results, and figures for this paper]

You should also download the SQLite database located on Zenodo at https://zenodo.org/record/2720105 (DOI for the SQLite database is 10.5281/zenodo.2720105
):
* ammonitor_bioacoustics.sqlite [a SQLite database containing data for the paper] 


Once these files are downloaded, open false-positives-script.R and follow the directions in the comments. Please reach out if you have any questions!
###############################################################################
#
# Script to accompany analysis and results from: 
#
# "Statistical learning mitigation of false positives from template-detected 
#  data in automated acoustic wildlife monitoring"
#
# by Cathleen M. Balantic & Therese M. Donovan
#
# Submitted to Bioacoustics - revised 03/13/2019
#
###############################################################################

Readers of this paper who are interested in reproducing the analysis should begin by downloading the following items from the Github repository:

* ammls [a folder]
   --> contains 'classifiers221.RDS' [an .RDS file that contains the trained/tested classifiers from this manuscript)
* ammonitor_bioacoustics.sqlite [a sqlite database containing data for the paper] **I AM NOT ABLE TO UPLOAD THIS TO GITHUB BECAUSE IT EXCEEDS THE FILESIZE LIMIT - FIGURE OUT AN ALTERNATIVE**
* ammonitor-functions.R [an R file that contains all AMMonitor functions necessary to produce the analysis/results in the paper]
* bioacoustics-submission-script.R [an R script to produce data, models, analysis, results, and figures for this paper]


Once these files are downloaded, open bioacoustics-submission-script.R.
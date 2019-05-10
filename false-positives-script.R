###############################################################################
#                                                                             #
# Script to accompany analysis and results from:                              #
#                                                                             #
# Balantic, C. M., & Donovan, T. M. (2019).                                   #
#                                                                             #
# Statistical learning mitigation of false positives                          #
# from template-detected data in automated acoustic wildlife monitoring       #
#                                                                             #
# Bioacoustics                                                                #
#                                                                             #
# https://www.tandfonline.com/doi/full/10.1080/09524622.2019.1605309          #
#                                                                             #
###############################################################################

# This code provides a means to reproduce analysis and results from the paper
# "Statistical learning mitigation of false positives from 
# template-detected data in automated acoustic wildlife monitoring" by 
# Cathleen M. Balantic & Therese M. Donovan.

# The AMMonitor R package referenced in this paper is not quite finalized at 
# this time, must undergo review, and is still private on Github. Thus, we have
# provided all necessary AMMonitor functions for this paper in an 
# accompanying file, 'ammonitor-functions.R', which can be sourced in to run 
# the following analysis. 


# SECTION 0: SET-UP ===========================================================
# BASIC SET-UP TO RUN THE ANALYSES IN THIS PAPER

# Install & call in necessary packages: 
library(AMModels)
library(caret)
library(data.table)
library(DBI)
library(monitoR) 
library(RSQLite)
library(ggplot2)

# Set your working directory to the folder that contains the files downloaded from Github
# setwd('your/folder/path/here')

# Source in AMMonitor functions from ammonitor-functions.R
source('ammonitor-functions.R')

# CONNECT TO THE SQLITE DATABASE: 
# Recall from the Github repository ReadMe that you need to download 
# the 'ammonitor_bioacoustics.sqlite' database located on Zenodo at
# https://zenodo.org/record/2720105 
# The DOI for the SQLite database is 10.5281/zenodo.2720105.
# Once downloaded into your current working directory, you will be 
# able to run the next three lines of code and connect to the database: 
db.name <- 'ammonitor_bioacoustics.sqlite'
db.path <- paste0(getwd(), '/', db.name)
conx <- RSQLite::dbConnect(drv = dbDriver('SQLite'), dbname = db.path)


# SECTION 1: DATA SUMMARIES ===================================================

# THIS SECTION GIVES A RUN-DOWN OF SOME BASIC DATA SUMMARIES FOR THE MANUSCRIPT
#  - total number of recordings
#  - total recording time in March 2016
#  - number of detections by each template
#  - spectrograms of detections by each template


# Total number of recordings: 40,094
recs <- data.table(RSQLite::dbReadTable(conx, 'recordings'))
nrow(recs) 

# Total recording time March 2016: 3256
march <- as.character(seq.Date(from = as.Date('2016-03-01', origin = '1970-01-01'), 
                               to = as.Date('2016-03-31', origin = '1970-01-01'),
                               by = 'day'))
recs[startDate %in% march, .N] # 3256/60 = 54.3 h

# Total number of detections: 10132
scores <- data.table(RSQLite::dbReadTable(conx, 'scores'))
scores[,.N] 

#  Detections for each focal template:
table(scores$templateID)
#  e7   g5   v1 
# 4427 1464 4241 


# CHECK ON THE VERIFICATIONS TO MAKE SURE ALL LOOKS WELL
template.ids <- sort(unique(scores$templateID))
thresholds <- c(0.43, 0.33, 0.23)

# The below loop uses plotVerifications() to plot individual verifications 
# for ECDO, GAQU, and VERD, with target signals plotted with green borders, and 
# false alarms plotted with red borders. The loop also uses plotVerificationsAvg()
# to plot the average spectrogram across all detections, all target signals, and
# all false alarms for each species. Note that these plots are operating for
# for ALL labeled data in the database now (not just March 2016 as in Fig. 5):
for (i in 1:3) {
  plotVerifications(db.path = db.path, 
                    templateID = template.ids[i], 
                    score.threshold = thresholds[i], 
                    label.type = 'speciesID', 
                    new.window = TRUE, box.lwd = 1)
  # note: x11() function may not work on Macs, try plot.new()
  x11(height = 6, width = 9) 
  plotVerificationsAvg(db.path,  
                       templateID = template.ids[i], 
                       score.threshold = thresholds[i], 
                       label.type = 'speciesID')
}


# SECTION 2: CREATE THE CLASSIFIER MODELS =====================================

# USE THIS SECTION TO TRAIN AND TEST THE CLASSIFIER MODELS FOR THE PAPER.

# Train and test strictly on the March 2016 data
recs <- data.table(
  dbGetQuery(conn = conx,
             statement = "SELECT recordingID, startDate FROM recordings"))
march <- as.character(seq.Date(from = as.Date('2016-03-01', origin = '1970-01-01'), 
                               to = as.Date('2016-03-31', origin = '1970-01-01'),
                               by = 'day'))
march.recs <- recs[startDate %in% march, recordingID]

# Find all scoreIDs from March 2016
scrs <- data.table(dbReadTable(conn = conx, 'scores'))
scrs.march <- scrs[recordingID %in% march.recs]

# Identify the March 2016 scoreIDs for each species
ecdo.ids <- scrs.march[templateID == 'e7', scoreID]
gaqu.ids <- scrs.march[templateID == 'g5', scoreID]
verd.ids <- scrs.march[templateID == 'v1', scoreID]

# Total number of detections in March
length(c(ecdo.ids, gaqu.ids, verd.ids)) #631

# Set seed to encourage reproducibility
#  Note that inherent stochasticity in the caret package's underlying
#  machine learning packages may not perfectly reproduce our results,
#  given the relatively small labeled training data sample sizes for
#  all three species
seed <- 221 

# Train and test classifier models for ECDO, GAQU, VERD:
# *** NOTE: THE TRAIN/TEST PROCESS MAY TAKE 10-20 MINUTES DEPENDING ON YOUR MACHINE
# *** SKIP TO SECTION 4 OF THE CODE TO SIMPLY READ IN OUR TRAINED/TESTED MODELS
# *** AND CONTINUE ON FROM THERE.

start <- Sys.time() # MAY TAKE 10 - 20 MINUTES!!
ecdo <- classifierModels(db.path = db.path,
                         scoreID = ecdo.ids,
                         label.type = 'speciesID',
                         split.proportion = 0.7,
                         seed = seed,
                         method = 'repeatedcv',
                         number = 10,
                         repeats = 1)

gaqu <- classifierModels(db.path = db.path,
                         scoreID = gaqu.ids,
                         label.type = 'speciesID',
                         split.proportion = 0.7,
                         seed = seed,
                         method = 'repeatedcv',
                         number = 10,
                         repeats = 1)

verd <- classifierModels(db.path = db.path, 
                         scoreID = verd.ids, 
                         label.type = 'speciesID', 
                         split.proportion = 0.7,      
                         seed = seed,
                         method = 'repeatedcv',
                         number = 10,
                         repeats = 1)
total.time <- Sys.time() - start
total.time 
# Took ~16 minutes to run. Much slower than our previous non-SQLite system,
# Perhaps due to serialization / unserialization process of scores in SQLite BLOB objects.

# Name all of these by hand and turn them into amModel objects
species.classifiers <- list(ecdo = ecdo, gaqu = gaqu, verd = verd)
all.mods <- list()
for (sc in seq_along(species.classifiers)) {
  classifiers <- species.classifiers[[sc]]
  mods <- list()
  for (i in seq_along(classifiers)) {
    mods[[i]] <- AMModels::amModel(model = classifiers[[i]], comment = '')
  }
  names(mods) <- names(classifiers)
  all.mods <- c(all.mods, mods)
}

# Create amModel library (amml) 
classifiers <- AMModels::amModelLib(description = 'Classifier models. Seed = 221.',
                          info = list(owner = 'Cathleen Balantic',
                                      email = 'cbalanti@uvm.edu'))

# Insert classifiers into amml:
classifiers <- AMModels::insertAMModelLib(models = all.mods, 
                                amml = classifiers)

# Save to ammls folder:
saveRDS(classifiers, paste0('ammls/classifiers221.RDS'))


# SECTION 3: READ IN OUR PAPER'S TRAINED/TESTED CLASSIFIER MODELS =============
classifiers <- readRDS('ammls/classifiers221.RDS')

# Initialize the model names for each species
ecdo.models <- c("e7_0.43_speciesID_glmnet", "e7_0.43_speciesID_svmLinear",
                 "e7_0.43_speciesID_svmRadial", "e7_0.43_speciesID_rf",
                 "e7_0.43_speciesID_kknn")
gaqu.models <- c("g5_0.33_speciesID_glmnet", "g5_0.33_speciesID_svmLinear",
                 "g5_0.33_speciesID_svmRadial", "g5_0.33_speciesID_rf",
                 "g5_0.33_speciesID_kknn")
verd.models <- c("v1_0.23_speciesID_glmnet", "v1_0.23_speciesID_svmLinear",
                 "v1_0.23_speciesID_svmRadial", "v1_0.23_speciesID_rf",
                 "v1_0.23_speciesID_kknn")

# SECTION 4: RECREATE TABLE 3 -- PERFORMANCE FROM TRAINING/CROSS-VALIDATION ====
# Create tables for performance results during the training (resampling/cross-val) phase:
ecdo.training <- gaqu.training <- verd.training <- list()
cnames <- c('glmnet', 'svmLinear', 'svmRadial', 'rf', 'kknn')

# Loop through all five models for all three species
for (i in 1:5) {
  ecdo.mod <- AMModels::getAMModel(amml = classifiers, x = ecdo.models[i])
  ecdo.training[[i]] <- confusionMatrix(data = ecdo.mod$training.fit$pred$pred, 
                                        reference = ecdo.mod$training.fit$pred$obs, 
                                        mode = "prec_recall")$byClass
  gaqu.mod <- AMModels::getAMModel(amml = classifiers, x = gaqu.models[i])
  gaqu.training[[i]] <- confusionMatrix(data = gaqu.mod$training.fit$pred$pred, 
                                        reference = gaqu.mod$training.fit$pred$obs, 
                                        mode = "prec_recall")$byClass
  verd.mod <- AMModels::getAMModel(amml = classifiers, x = verd.models[i])
  verd.training[[i]] <- confusionMatrix(data = verd.mod$training.fit$pred$pred, 
                                        reference = verd.mod$training.fit$pred$obs, 
                                        mode = "prec_recall")$byClass
}

# Rearrange results into tables
ecdo.tr <- matrix(data = unlist(ecdo.training), nrow = 5, ncol = 11, byrow = TRUE)
rownames(ecdo.tr) <- cnames
colnames(ecdo.tr) <- names(ecdo.training[[1]])

gaqu.tr <- matrix(data = unlist(gaqu.training), nrow = 5, ncol = 11, byrow = TRUE)
rownames(gaqu.tr) <- cnames
colnames(gaqu.tr) <- names(gaqu.training[[1]])

verd.tr <- matrix(data = unlist(verd.training), nrow = 5, ncol = 11, byrow = TRUE)
rownames(verd.tr) <- cnames
colnames(verd.tr) <- names(verd.training[[1]])

# Table 3 results: 
for.table <- c('Sensitivity', 'Specificity', 'Pos Pred Value', 'F1')
ecdo.tr[,for.table]
gaqu.tr[,for.table]
verd.tr[,for.table]


# SECTION 5: VARIABLE IMPORTANCE RESULTS ======================================

# This section uses the caret function varImp() to investigate variable 
# importance for all three species detection models, focusing on the 
# regularized logistic regression (glmnet) classifier and the random forests (rf)
# classifiers which provide the most readily interpretable outputs

# ECDO: 
imps <- list()
for (i in 1:5) {
  mod <- AMModels::getAMModel(amml = classifiers, x = ecdo.models[i])
  imps[[i]] <- caret::varImp(object = mod$training.fit, useModel = TRUE, scale = FALSE)
}
names(imps) <- ecdo.models

# Variable importance for glmnet
# sp.sem (1), 
# sp.sd, sp.mean, sp.cent, tc.83
# sp.mode, some more tcs
plot(imps[[1]], top = 50, col = 'black',
     main = 'ECDO Variable Importance (glmnet)')

# Variable importance for random forest
# sp.mode hugely the most important
# sp.cent and sp.mean contributing
plot(imps[[4]], top = 50, col = 'black',
     main = 'ECDO Variable Importance (rf)')

# GAQU
imps <- list()
for (i in 1:5) {
  mod <- AMModels::getAMModel(amml = classifiers, x = gaqu.models[i])
  imps[[i]] <- caret::varImp(object = mod$training.fit, useModel = TRUE, scale = FALSE)
}

names(imps) <- gaqu.models

# Variable importance for glmnet
# kurtosis & skewness, 
# zero-crossing rates
plot(imps[[1]], top = 50, col = 'black',
     main = 'GAQU Variable Importance (glmnet)')

# Variable importance for random forest
# amplitude values
plot(imps[[4]], top = 50, col = 'black',
     main = 'GAQU Variable Importance (rf)')


# VERD
imps <- list()
for (i in 1:5) {
  mod <- AMModels::getAMModel(amml = classifiers, x = verd.models[i])
  imps[[i]] <- caret::varImp(object = mod$training.fit, useModel = TRUE, scale = FALSE)
}
names(imps) <- verd.models

# Variable importance for glmnet
# spectral entropy, flatness, and kurtosis
# zc & tc values
plot(imps[[1]], top = 50, col = 'black',
     main = 'VERD Variable Importance (glmnet)')

# Variable importance for random forest
# skewness, correlation score
plot(imps[[4]], top = 50, col = 'black',
     main = 'VERD Variable Importance (rf)')



# SECTION 6: RECREATE TABLE 4 -- PERFORMANCE ON TEST DATA =====================
cols <- c('Sensitivity', 'Specificity', 'Pos Pred Value', 'F1')

# Gather weighted ensemble averages

# Subset the test scoreIDs (they are the same for each classifier)
ecdo.test.ids <- classifiers@models$e7_0.43_speciesID_glmnet@model$test.scoreID
gaqu.test.ids <- classifiers@models$g5_0.33_speciesID_glmnet@model$test.scoreID
verd.test.ids <- classifiers@models$v1_0.23_speciesID_glmnet@model$test.scoreID

# Create a species scoreIDs, templateIDs, and score.thresholds list for looping through
sp <- list(ecdo = list(test.id = ecdo.test.ids, templateID = 'e7', score.threshold = 0.43),
           gaqu = list(test.id = gaqu.test.ids, templateID = 'g5', score.threshold = 0.33),
           verd = list(test.id = verd.test.ids, templateID = 'v1', score.threshold = 0.23))
ensembles <- c('sensitivity', 'specificity', 'precision', 'f1')
cols <- c('Sensitivity', 'Specificity', 'Pos Pred Value', 'F1')

# Loop through species to make predictions and create ensemble results
ens.perf <- list()
for (i in 1:length(sp)) {
  
  # Gather predictions for this species
  preds.march <- classifierPredict(db.path = db.path,
                                   amml = classifiers,
                                   date.range = c('2016-03-01', '2016-03-31'),
                                   templateID = sp[[i]]$templateID, 
                                   label.type = 'speciesID',
                                   score.threshold = sp[[i]]$score.threshold,
                                   classifiers = c('glmnet', 'svmLinear', 
                                                   'svmRadial', 'rf', 'kknn'),
                                   db.insert = FALSE)
  
  # Gather classifications on the test data 
  test.DT <- preds.march[scoreID %in% sp[[i]]$test.id]
  
  # Loop through the four weighted average ensembles
  
  sp.ens.perf <- matrix(data = 0, nrow = 4, ncol = 4)
  for (j in 1:length(ensembles)) {
    
    # Create weighted average 'ensemble'
    ensemble.classifier <- classifierEnsemble(db.path = db.path, 
                                              amml = classifiers,
                                              classificationsDT = test.DT,
                                              ensemble = ensembles[j])
    # Assess performance of the ensemble
    sp.ens.perf[j,] <- caret::confusionMatrix(data = ensemble.classifier$predicted,
                                       reference = ensemble.classifier$class,
                                       positive = '1')$byClass[cols]
  }
  rownames(sp.ens.perf) <- paste0('wt.by.', ensembles)
  colnames(sp.ens.perf) <- cols
  
  # Store all ensemble performance for this species
  ens.perf[[i]] <- sp.ens.perf
}
names(ens.perf) <- c('ECDO', 'GAQU', 'VERD')

# View performance on the test data for the 5 classifiers
classifierPerformance(amml = classifiers,
                      model.names = ecdo.models)[,cols, with = FALSE]
classifierPerformance(amml = classifiers,
                      model.names = gaqu.models)[,cols, with = FALSE]
classifierPerformance(amml = classifiers,
                      model.names = verd.models)[,cols, with = FALSE]

# View performanace on the test data for the 4 'ensembles'
ens.perf


# SECTION 7: RECREATE FIGURE 6 -- ROC CURVES ==================================

# ROC curves of training and testing data
dev.off()
x11(height = 6, width = 7) # note: x11() may not work on Macs
windowsFonts(Times = windowsFont("Times New Roman"))
par(mfrow = c(2,3), family = 'Times')
font <- 6 
model.names <- list(ECDO = ecdo.models, GAQU = gaqu.models, VERD = verd.models)
for (i in c('train', 'test')) {
  if (i == 'train') ens <- NULL; if (i == 'test') ens <- 'f1'
  for (j in 1:3) {
    plotROC(db.path = db.path, amml = classifiers, 
            model.names = model.names[[j]],
            curve.type = 'roc', data.type = i, main = names(model.names)[j],
            ensembles = ens, lwd = 1, cex.lab = 1.25)
  }
  if (i == 'train')
    mtext('ROC Curves of Training Data', side = 3, 
          line = -2.25, outer = TRUE, font = font)
  else
    mtext('ROC Curves of Test Data', side = 3, 
          line = -25, outer = TRUE, font = font)
}


# SECTION 8: RECREATE FIGURE 7 -- PRECISION-RECALL CURVES =====================

par(family = 'Times')
layout(mat = matrix(c(1,2,3, 
                      4,5,6)), 
       widths = rep(2.66666667, 6), 
       heights = rep(4,6))
x11(height = 8, width = 8) # note: x11() may not work on Macs
windowsFonts(Times = windowsFont("Times New Roman"))
par(mfrow = c(2,3), family = 'Times')
font <- 6 
model.names <- list(ECDO = ecdo.models, GAQU = gaqu.models, VERD = verd.models)
for (i in c('train', 'test')) {
  if (i == 'train') ens <- NULL; if (i == 'test') ens <- 'f1'
  for (j in 1:3) {
    plotROC(db.path = db.path, amml = classifiers, 
            model.names = model.names[[j]],
            curve.type = 'pr', data.type = i, ensembles = ens,
            main = names(model.names)[j], lwd = 1, cex.lab = 1.25)
  }
  if (i == 'train')
    mtext('Precision-Recall Curves of Training Data', side = 3, 
          line = -2.25, outer = TRUE, font = font)
  else
    mtext('Precision-Recall Curves of Test Data', side = 3, 
          line = -25, outer = TRUE, font = font)
}


# SECTION 9: GATHER PREDICTIONS ON UNSEEN DATA -- APRIL 2016 - MAY 2017: =========

# When using the classifierPredict() function below: 
#
# Set db.insert = TRUE to insert classifications into the database
# (this has already been done and you don't need to do it again - you will 
# encounter a foreign key constraint scolding from the SQLite database if you
# attempt to insert classifications again)
#
# Set db.insert = FALSE to simply return a data.table of predictions

# This section takes 25-30 seconds to run

# ECDO 
pred.ecdo <- classifierPredict(db.path = db.path,
                               amml = classifiers,
                               date.range = c('2016-04-01', '2017-05-31'),
                               templateID = 'e7', 
                               label.type = 'speciesID',
                               score.threshold = 0.43,
                               classifiers = c('glmnet', 'svmLinear', 
                                               'svmRadial', 'rf', 'kknn'),
                               db.insert = FALSE)
# GAQU
pred.gaqu <- classifierPredict(db.path = db.path,
                               amml = classifiers,
                               date.range = c('2016-04-01', '2017-05-31'),
                               templateID = 'g5',
                               label.type = 'speciesID',
                               score.threshold = 0.33,
                               classifiers = c('glmnet', 'svmLinear', 
                                               'svmRadial', 'rf', 'kknn'),
                               db.insert = FALSE)
# VERD 
pred.verd <- classifierPredict(db.path = db.path, 
                               amml = classifiers, 
                               date.range = c('2016-04-01', '2017-05-31'),
                               templateID = 'v1', 
                               label.type = 'speciesID', 
                               score.threshold = 0.23, 
                               classifiers = c('glmnet', 'svmLinear', 
                                               'svmRadial', 'rf', 'kknn'), 
                               db.insert = FALSE)


# SECTION 10: CREATE ENSEMBLE ON UNSEEN DATA -- APRIL 2016 - MAY 2017 =========
ens.pred.ecdo <- classifierEnsemble(db.path = db.path,
                                    classificationsDT = pred.ecdo,
                                    amml = classifiers,
                                    ensemble = 'f1')

ens.pred.gaqu <- classifierEnsemble(db.path = db.path,
                                    classificationsDT = pred.gaqu,
                                    amml = classifiers,
                                    ensemble = 'f1')

ens.pred.verd <- classifierEnsemble(db.path = db.path, 
                                    classificationsDT = pred.verd, 
                                    amml = classifiers, 
                                    ensemble = 'f1')

# SECTION 11: ASSESS PERFORMANCE ON UNSEEN DATA -- APRIL 2016 - MAY 2017: =====

# This section generates the results used in Figure 8

cols <- c('Sensitivity', 'Specificity', 'Pos Pred Value', 'F1')

# PERFORMANCE OF THE CLASSIFICATION SYSTEM:
results.ecdo <- confusionMatrix(data = ens.pred.ecdo$predicted,
                                reference = ens.pred.ecdo$class,
                                positive = '1')$byClass[cols]
results.gaqu <- confusionMatrix(data = ens.pred.gaqu$predicted,
                                reference = ens.pred.gaqu$class,
                                positive = '1')$byClass[cols]
results.verd <- confusionMatrix(data = ens.pred.verd$predicted, 
                                reference = ens.pred.verd$class, 
                                positive = '1')$byClass[cols]
results.ecdo
results.gaqu
results.verd


# PERFORMANCE OF THE TEMPLATE-ONLY SYSTEM:
# Sensitivity == 1
# Specificity == 0

# precision: TP / TP + FP   
# f1: 2*(sensitivity*precision)/(sensitivity + precision)

# Precision calculcations for ecdo, gaqu, verd templates:
(tmp.prec.ecdo <- nrow(ens.pred.ecdo[class == 1 ])/nrow(ens.pred.ecdo))
(tmp.prec.gaqu <- nrow(ens.pred.gaqu[class == 1 ])/nrow(ens.pred.gaqu))
(tmp.prec.verd <- nrow(ens.pred.verd[class == 1 ])/nrow(ens.pred.verd))

# F1 score calculations for ecdo, gaqu, verd templates:
(tmp.f1.ecdo <- 2*(1*tmp.prec.ecdo)/(1 + tmp.prec.ecdo))
(tmp.f1.gaqu <- 2*(1*tmp.prec.gaqu)/(1 + tmp.prec.gaqu))
(tmp.f1.verd <- 2*(1*tmp.prec.verd)/(1 + tmp.prec.verd))


# SECTION 12: RECREATE FIGURE 8 -- Classification vs. Template-based Results ====

# Gather classification-based results
results.cl <- rbind(results.ecdo, results.gaqu, results.verd)

# Gather template-based results (sens, spec, prec, f1)
#   (template-based sensitivity is 1, specificity is 0)
f1.tmps <- c(tmp.f1.ecdo, tmp.f1.gaqu, tmp.f1.verd)
prec.tmps <- c(tmp.prec.ecdo, tmp.prec.gaqu, tmp.prec.verd)

nms <- c('ECDO', 'GAQU', 'VERD')
mets <- c(rep('Sensitivity',3),
          rep('Pos. Pred. Value',3),
          rep('F1', 3),
          rep('Specificity', 3))
dat <- data.frame(Template = rep(c('ECDO', 'GAQU', 'VERD'),8),
                  Type = c(rep('Classification Phase',12), 
                           rep('Template Screening Phase',12)),
                  variable = rep(mets,2),
                  value = c(results.cl[,'Sensitivity'], # classification system
                            results.cl[,'Pos Pred Value'], # classification system
                            results.cl[,'F1'], # classification system
                            results.cl[,'Specificity'], # classification system
                            rep(1,3), # template only
                            prec.tmps, # template only
                            f1.tmps, # template only
                            rep(0,3))) # template only
dat$variable <- factor(dat$variable, levels(dat$variable)[c(3,4,2,1)])

# Comparison of Metrics Plot 
mytheme <- theme(axis.title.y = element_blank(),
                 axis.ticks.y = element_blank(), panel.border = element_blank(),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 axis.line = element_line(color = 'black'), text = element_text(family = 'Times'),
                 plot.title = element_text(hjust = 0.5, size = 13),
                 legend.text = element_text(size = 12), strip.text.x = element_text(size = 13),
                 legend.position = 'top',
                 legend.title = element_blank(), 
                 axis.text = element_text(color = 'black'))

# https://github.com/tidyverse/ggplot2/wiki/labeller
windowsFonts(Times = windowsFont("Times New Roman"))
fig8 <- ggplot(dat, aes(Template, value)) +
  geom_bar(aes(fill = Type), position = 'dodge', stat = 'identity') +
  facet_grid(.~ variable) +
  scale_fill_manual(values = c("gray33", 'gray66')) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  ggtitle('Performance Metrics to be Maximized') +
  ylab('') + xlab('') + theme_bw() + mytheme 
x11(height = 2.5, width = 6)
fig8





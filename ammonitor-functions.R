###############################################################################
#                                                                             #
# Functions to accompany analysis and results from:                           #
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

# NOTE: This file does not include all AMMonitor package functions. 
# It includes only those necessary to run the manuscript experiment located in
# false-positives-script.R
# Package data not included for function example documentation at this time; 
# contact me if this is of interest to you.

# Function list: 
# classifierAssess
# classifierEnsemble
# classifierModels
# classifierPerformance
# classifierPredict
# classifierTest
# classifierTrain
# plotROC
# plotVerifications
# plotVerificationsAvg
# pr
# qryPkCheck
# qryDeployment
# scoresVerify
# templatesUnserialize


#' @name classifierAssess
#' @title Assess classifier performance on test data
#' @description Internal function. In practice, please use \code{\link{classifierModels}} instead.
#' @param test.predictions Test predictions list object produced by \code{\link{classifierTest}}.
#' @return A length-two nested list of the following data:
#' \enumerate{
#' \item{\strong{summary}: a data.frame of 18 columns and a number of rows equal to the number of classifiers used. The 18 columns contain classifier assessment metrics produced by \strong{caret}, which are: Accuracy, Kappa, AccuracyLower, AccuracyUpper, AccuracyNull, AccuracyPValue, McnemarPValue, Sensitivity, Specificity, Pos Pred Value, Neg Pred Value, Precision, Recall, F1, Prevalence, Detection Rate, Detection Prevalence, and Balanced Accuracy.}
#' \item{\strong{confusion.matrices}: a list of confusion matrices for each classifer.}
#' }
#' @seealso \code{\link{classifierModels}}
#' @export
#' @examples
#' \dontrun{
#' 
#' # ------------------------------------------------------------
#' # Set up a demo AMMonitor database
#' # ------------------------------------------------------------
#' 
#' # Create the database (this will be deleted):
#' db.name <- 'demo.sqlite'
#' 
#' dbCreateSample(db.name = db.name, 
#'                file.path = paste0(getwd(),'/'), 
#'                tables = c('recordings', 'templates', 'locations', 
#'                           'equipment', 'people', 'accounts', 
#'                           'library', 'species', 'scores'))
#'                 
#' # Verify that the database exists in your current working directory:
#' file.exists(db.name)
#' 
#' # Connect to the db: 
#' db.path <- paste0(getwd(), '/', db.name)
#' conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
#' 
#' # ---------------------------------------------------------------
#' # Add in verifications at the libraryID and speciesID level.
#' # Below, we manually input verification labels. (Note that in 
#' # practice, users should acquire verifications using scoresVerify(). 
#' # We update the verification labels manually here to easily 
#' # demonstrate the function.)
#' #----------------------------------------------------------------
#' 
#' # Update verifications in the manualVerifyLibraryID column
#' # of the scores table:
#' lib.verifications <- c(0,0,0,0,0,1,0,1,1,1,0,1,0,1,0,1,1,1,1,1,1,1)
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyLibraryID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = lib.verifications, scoreID = scoreIDs))
#'                      
#' # Update verifications in the manualVerifySpeciesID column
#' # of the scores table (remember, these will be different from the libraryID)
#' sp.verifications <- rep(1, length(scoreIDs))
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyspeciesID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = sp.verifications, scoreID = scoreIDs))
#'
#' #----------------------------------------------------------------
#' # Create Train/Test Data Partition
#' #----------------------------------------------------------------
#' 
#' # Create Data Partition (train/test)
#' set.seed(3)
#' training.index <- createDataPartition(y = verifications, p = 0.7, 
#'                                       list = FALSE, times = 1)
#' train.scoreID <- scoreIDs[training.index]
#' test.scoreID <- scoreIDs[-training.index]
#' 
#' #----------------------------------------------------------------
#' # Train the classifiers
#' #----------------------------------------------------------------
#' fits <- classifierTrain(db.path = db.path, 
#'                         train.scoreID = train.scoreID, 
#'                         label.type = 'libraryID',
#'                         classifiers = c('glmnet', 'svmLinear', 
#'                                         'svmRadial', 'rf', 'kknn'),
#'                         seed = 10)
#' 
#' #----------------------------------------------------------------
#' # Test the classifiers 
#' #----------------------------------------------------------------
#' test <- classifierTest(db.path = db.path, 
#'                        test.scoreID = test.scoreID, 
#'                        training.fits = fits)
#'                        
#' #----------------------------------------------------------------
#' # Assess the classifiers 
#' #----------------------------------------------------------------
#' assess <- classifierAssess(test.predictions = test)
#' 
#' # ------------------------------------------------------------
#' # Disconnect from demo database and delete from working directory
#' # ------------------------------------------------------------
#'
#' dbDisconnect(conx)
#' unlink(db.path)
#'}
#'

classifierAssess <- function(test.predictions) {
  
  # Grab the scoreIDs used from testing:
  test.scoreID <- test.predictions$test.scoreID
  
  # Acquire scores based on input scoreIDs
  y.test <- test.predictions$y.test
  
  # Extract the raw predictions
  raw.preds <- test.predictions$test.predictions
  
  # Confusion matrix summary data for the individual classifiers:
  summ <- data.frame(matrix(0, nrow = length(raw.preds), ncol = 18))
  cm <- list()
  for (i in 1:length(raw.preds)) {
    cm[[i]] <- confusionMatrix(data = raw.preds[[i]]$Class, reference = y.test, positive = 'TS')
    summ[i, 1:7] <- cm[[i]]$overall
    summ[i, 8:18] <- cm[[i]]$byClass
  }
  colnames(summ) <- c(names(cm[[1]]$overall), names(cm[[1]]$byClass))
  rownames(summ) <- names(raw.preds)
  summ <- round(summ, 4)
  
  # Reorganize individual classifier confusion matrices
  mats <- lapply(cm, '[[', 'table') 
  names(mats) <- rownames(summ)
  
  # Return metrics performance summary and confusion matrices
  assess.list <- list(summary = summ, confusion.matrices = mats)
  return(assess.list)
  
}





#' @name classifierEnsemble
#' @title Create weighted-average ensemble classifications
#' @description Generate a data.table of 'ensemble' classifications averaged across multiple classifiers. The 'ensemble' may be weighted according to classifier performance on accuracy, sensitivity, specificity, precision, F1 score, or a simple average across all classifiers. \code{classifierEnsemble} requires either the amml or model.list argument (not both). In addition to the amml or model.list arguments, users can optionally input a classificationsDT to directly specify which records should have an ensemble applied to them. If using the amml argument without a classificationsDT, \code{classifierEnsemble} returns ensemble results for \strong{all} records in the database according to the user-specified model.names. If using the model.list argument without a classificationsDT, \code{classifierEnsemble} returns ensemble results strictly for the 'test' data contained within the model.list. 
#' @param db.path The file path that houses the SQLite database.
#' @param amml Default = NULL. An \code{amModelLib} object. Not necessary if using the \code{model.list} argument.
#' @param model.list Default = NULL. A list of classification models returned by \code{\link{classifierModels}}. Not necessary if using the \code{amml} argument. 
#' @param model.names Default = NULL. Character vector of the model names from which to generate ensembles. Model names should all come from the same templateID_score.threshold_label.type origin (e.g.: c('verd1_0.2_speciesID_glmnet', 'verd1_0.2_speciesID_svmLinear')). Not required if using classificationsDT argument.
#' @param ensemble Default = 'f1'. Character string indicating a single ensemble type to use. Options: c('accuracy', 'sensitivity', 'specificity', 'precision', 'f1', 'simple').
#' @param classificationsDT Default = NULL. Optional: a data.table of classifications from desired classifier models, from which to generate ensembles. In practice, the classificationsDT data.table is produced by a call to classifierPredict() (see \code{\link{examples}}). The 'modelNames' column in this classifications data.table should contain model names that come from the same templateID_score.threshold_label.type origin (e.g.: c('verd1_0.2_speciesID_glmnet', 'verd1_0.2_speciesID_svmLinear')). The model.names argument is not required if using classificationsDT. The amml argument OR the model.list argument IS required if using classificationsDT. 
#' @return A classifications data.table providing ensemble results. Columns are \strong{scoreID}, \strong{classifier} (which is the ensemble.type), \strong{modelProbability} (the weighted average probability of target signal returned by all classifiers input to model.names), \strong{class} (the true class for this signal, as 1: target signal, or 0: false alarm), and \strong{predicted} (the predicted class for this signal based on the ensemble model probability). The \strong{class} and \strong{predicted} columns provide easy inputs to the \code{\link[caret]{confusionMatrix}} function for assessing performance of the ensemble classifier (see \strong{Examples}).
#' @seealso \code{\link{shapeOccupancy}}
#' @export
#' @examples
#' 
#' \dontrun{
#' # ------------------------------------------------------------
#' # Set up a demo AMMonitor database
#' # ------------------------------------------------------------
#'
#' # Create the database (this will be deleted):
#' db.name <- 'demo.sqlite'
#' 
#' dbCreateSample(db.name = db.name,
#'                file.path = paste0(getwd(),'/'),
#'                tables = c('recordings', 'templates', 'locations',
#'                           'equipment', 'people', 'accounts',
#'                           'library', 'species', 'scores',
#'                           'classifications'))
#'
#' # Verify that the database exists in your current working directory:
#' file.exists(db.name)
#'
#' # Connect to the db:
#' db.path <- paste0(getwd(), '/', db.name)
#' conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
#'
#' # ---------------------------------------------------------------
#' # Add in verifications at the libraryID and speciesID level.
#' # Below, we manually input verification labels. (Note that in 
#' # practice, users should acquire verifications using scoresVerify(). 
#' # We update the verification labels manually here to easily 
#' # demonstrate the function.)
#' #----------------------------------------------------------------
#' 
#' # Update verifications in the manualVerifyLibraryID column
#' # of the scores table:
#' lib.verifications <- c(0,0,0,0,0,1,0,1,1,1,0,1,0,1,0,1,1,1,1,1,1,1)
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyLibraryID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = lib.verifications, scoreID = scoreIDs))
#'                      
#' # Update verifications in the manualVerifySpeciesID column
#' # of the scores table (remember, these will be different from the libraryID)
#' sp.verifications <- rep(1, length(scoreIDs))
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyspeciesID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = sp.verifications, scoreID = scoreIDs))
#'
#'
#' #----------------------------------------------------------------
#' # classifierEnsemble for classifier models contained within an 
#' # amml, using classifications that already exist in the database
#' #----------------------------------------------------------------
#',
#' # Read in a sample amml object:
#' data(classifiers_amml)
#' 
#' ens <- classifierEnsemble(db.path = db.path,
#'                    amml = classifiers_amml,
#'                    model.list = NULL,
#'                    model.names = c('verd1_0.2_libraryID_glmnet',
#'                                    'verd1_0.2_libraryID_svmRadial',
#'                                    'verd1_0.2_libraryID_rf'),
#'                    ensemble = 'simple',
#'                    classificationsDT = NULL)
#' head(ens)
#' 
#' # Apply the confusionMatrix function to the ensemble results
#' caret::confusionMatrix(data = ens$predicted, reference = ens$class)
#'
#' #----------------------------------------------------------------
#' # classifierEnsemble using a model.list (where the ensembles
#' # are computed strictly on the test data used in classifierModels)
#' #----------------------------------------------------------------
#'
#' # Read in a model list object produced by classifierModels():
#' data(classifier_practice)
#'
#' ens <- classifierEnsemble(db.path = db.path,
#'                           amml = NULL,
#'                           model.list = classifier_practice,
#'                           model.names = c('verd1_0.2_libraryID_svmRadial',
#'                                           'verd1_0.2_libraryID_svmLinear',
#'                                           'verd1_0.2_libraryID_kknn'),
#'                           ensemble = 'precision',
#'                           classificationsDT = NULL)
#' head(ens)
#' 
#' # Apply the confusionMatrix function to the ensemble results
#' caret::confusionMatrix(data = ens$predicted, reference = ens$class)
#' 
#' #----------------------------------------------------------------
#' # classifierEnsemble using the classificationsDT argument
#' #----------------------------------------------------------------
#' 
#' # Read in sample classifier models amml:
#' data(classifiers_amml)
#' 
#' # Predict class of new detections:
#' new.preds <- classifierPredict(db.path = db.path, 
#'                                amml = classifiers_amml,
#'                                date.range = c('2016-03-01', 
#'                                               '2016-03-30'),
#'                                templateID = 'verd1', 
#'                                label.type = 'libraryID', 
#'                                score.threshold = 0.2, 
#'                                classifiers =  c('glmnet', 'svmRadial',  
#'                                                 'svmLinear', 'rf', 
#'                                                 'kknn'), 
#'                                db.insert = FALSE)
#'                                
#' ens <- classifierEnsemble(db.path = db.path,
#'                           amml = classifiers_amml,
#'                           model.list = NULL,
#'                           model.names = NULL,
#'                           ensemble = 'sensitivity',
#'                           classificationsDT = new.preds)
#' head(ens)  
#' 
#' # Apply the confusionMatrix function to the ensemble results
#' caret::confusionMatrix(data = ens$predicted, reference = ens$class)       
#' 
#' # ------------------------------------------------------------
#' # Disconnect from demo database and delete from working directory
#' # ------------------------------------------------------------
#'
#' dbDisconnect(conx)
#' unlink(db.path)
#'         
#' }
#'
#' 
classifierEnsemble <- function(db.path,
                               amml = NULL,
                               model.list = NULL,
                               model.names = NULL,
                               ensemble = 'f1',
                               classificationsDT = NULL
){
  
  # Connect to the database:
  conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
  
  # Turn the SQLite foreign constraints on
  rs <- dbSendQuery(conn = conx, statement = "PRAGMA foreign_keys = ON;" )
  dbClearResult(rs)
  
  # Check inputs
  if (!is.null(model.list) & !is.null(amml)) 
    stop('Please use either the amml OR model.list argument. Not both.')
  if (is.null(classificationsDT) & is.null(model.names))
    stop('Please use the model.names argument if not using the classificationsDT argument.')
  
  # If using an amml
  if (!is.null(amml)) {
    if (!is.null(classificationsDT)) {
      model.names <- unique(classificationsDT$modelName)
    } else {
      classificationsDT <-
        data.table(
          dbGetQuery(conn = conx, statement =
                       'SELECT * FROM classifications
                     WHERE modelName = $modelName',
                     param = list(modelName = model.names)))
    }
  }
  
  # If using model.list
  if (!is.null(model.list)) {
    
    # If not inputting a classificationsDT, create one
    if (!is.null(classificationsDT)) {
      model.names <- unique(classificationsDT$modelName)
    } else {
      # Create a classificationsDT
      model.list <- model.list[model.names]
      spl <- sapply(model.names, 
                    function(x) strsplit(x, split = '_', fixed = TRUE))
      classifiers <- unique(unlist(lapply(spl, '[[', 4)))
      ids <- unlist(lapply(model.list,'[[', 'test.scoreID'))
      classificationsDT <- data.table(scoreID = ids,
                                      classifier = rep(classifiers, 
                                                       each = length(unique(ids))),
                                      modelName = rep(model.names, 
                                                      each = length(unique(ids))),
                                      modelProbability = as.numeric(NA))
      
      # Add test predictions for target signal modelProbability
      for (l in 1:length(model.names)) {
        classificationsDT[modelName == model.names[l],
                          modelProbability := model.list[[model.names[l]]]$test.prediction$TS]
      }
    } # end if null classifications DT
  }
  
  # Trim classificationsDT to keep only relevant columns
  classif <- classificationsDT[,c('scoreID', 'classifier', 'modelProbability')]
  
  # Parse model.names for useful info
  spl <- sapply(model.names, 
                function(x) strsplit(x, split = '_', fixed = TRUE))
  score.threshold <- unique(unlist(lapply(spl, '[[', 2)))
  label.type <- unique(unlist(lapply(spl, '[[', 3)))
  classifiers <- unique(unlist(lapply(spl, '[[', 4)))
  
  # Ensure scoresIDS will all be in the same order for each classifier
  setkeyv(classif, c('classifier', 'scoreID')) 
  scoreIDs <- sort(unique(classif$scoreID))
  
  # Get the performance summary for all desired classifiers
  if (!is.null(amml))
    perf.summary <- classifierPerformance(amml = amml, model.names = model.names)
  if (!is.null(model.list))
    perf.summary <- classifierPerformance(model.list = model.list, 
                                          model.names = model.names)
  
  # Compute a simple probability average across classifiers if 'simple'
  if (ensemble == 'simple') {
    ens.prob <- classif[,mean(modelProbability), by = scoreID][,V1]
  } else {
    
    # Otherwise, perform weighted average
    
    # Convert caret-generated column names to lower case
    colnames(perf.summary) <- tolower(colnames(perf.summary))
    
    # Find each classifier's score on this performance metric:
    metric.scores <- perf.summary[,get(ensemble)]
    names(metric.scores) <- classifiers
    
    # If any classifier has a NaN or NA for this metric,
    #   (e.g. if the metric is precision, we can get a NaN on precision due to a 
    #   classifier scoring 0 on sensitivity) 
    # then we convert NaNs and NAs to zero to make sure that classifier doesn't
    # get to contribute to the weighted average vote on this metric. 
    # (Not sure if this is an acceptable way to deal with this, but it makes sense to me). 
    metric.scores[is.nan(metric.scores)] <- 0
    metric.scores[is.na(metric.scores)] <- 0
    
    # Compute a vector representing how proportionally
    #  close each score is to the highest score,
    proportion.of.max <- metric.scores/max(metric.scores)
    
    # Compute a vector of weights normalized to add to 1
    wts <- proportion.of.max/(sum(proportion.of.max))
    
    # Reshape classif to wide format in order to compute dot-product
    #   (dcast/reshape aren't working as expected, so a loop it is)
    pr.mat <- matrix(0, nrow = length(unique(classif$scoreID)), ncol = length(classifiers))
    colnames(pr.mat) <- classifiers
    for (i in 1:length(classifiers)) {
      pr.mat[,classifiers[i]] <- classif[classifier == classifiers[i], modelProbability]
    }
    
    # Compute dot-product of wts & prob vectors to
    #  acquire single ensemble output probability value
    ens.prob <- pr.mat %*% wts
  }
  
  # Begin a data.table of classifications to return
  dt <- data.table(scoreID = scoreIDs,
                   classifier = paste0('ensemble.', ensemble),
                   modelProbability = ens.prob)
  
  # For some reason data.table is not respecting the name "modelProbability" 
  # and is changing it to modelProbability.V1. Change it back
  colnames(dt)[colnames(dt) == 'modelProbability.V1'] <- 'modelProbability'
  
  # Add verified labels
  scrs <- data.table(
    dbGetQuery(conn = conx,
               statement = "SELECT scoreID, templateID, scoreThreshold, score, 
               manualVerifyLibraryID, manualVerifySpeciesID, features 
               FROM scores 
               WHERE scoreID = $scoreID",
               params = list(scoreID = dt$scoreID)))
  ifelse(label.type == 'speciesID', 
         type <- 'manualVerifySpeciesID', 
         type <- 'manualVerifyLibraryID')
  verified.labels <- scrs[,c('scoreID', type), with = FALSE]
  colnames(verified.labels)[colnames(verified.labels) == type] <- 'class'
  dt <- merge(x = dt, y = verified.labels, by = 'scoreID', all = TRUE)
  dt[,predicted := ifelse(modelProbability < 0.5, 0, 1)]
  
  # Turn class & predicted into factors for confusionMatrix function
  dt[,class := factor(class, levels = c('1', '0'), labels = c('1', '0'))]
  dt[,predicted := factor(predicted, levels = c('1', '0'), labels = c('1', '0'))]
  
  return(dt)
}




#' @name classifierModels
#' @title Train and test a suite of classifier models
#' @description Use verified scores to train and test a suite of classifier models that distinguish between target signals and false alarms in an automated acoustic monitoring program. Each classifier model is tuned to produce the probability that a detected event is a target signal. Models that show promise may be added to an AMModels database and used again later.
#' @param db.path The file path that houses the SQLite database.
#' @param templateID Character string of the templateID for which to create classifier models. 
#' @param label.type Character string indicating whether to generate classifiers using the speciesID or the libraryID. Options = c('speciesID', 'libraryID').
#' @param score.threshold Numeric score.threshold used for this templateID. This value should match a score.threshold you have already used for this templateID to populate the scores table. 
#' @param scoreID Optional vector of scoreIDs that grants finer control over which data should be used for creating the classifier models. This option is allowed if not using the templateID, label.type, and score.threshold arguments (which will automatically use ALL relevant labeled data for creating the classifier models). This vector of scoreIDs enables the user to directly select which labeled scoreIDs should be used for training and testing. Note that any scoreIDs input to this argument must all be linked to the same template and score threshold, and the scoreID must be 'labeled' with a 0 or 1 in the manualVerifyLibraryID or manualVerifySpeciesID column based on which label.type you have indicated in the label.type argument. Note that when using the scoreID argument, the label.type argument is still required. The templateID and score.threshold arguments are not required. 
#' @param split.proportion Numeric value input to \code{\link[caret]{createDataPartition}} indicating the proportion of input scoreIDs that should be used for training (one minus this value will be used during the test phase). For example, if 'split.proportion' = 0.7, 70\% of the scoreIDs will be used during the training phase, and the remaining 30\% will be used during the testing phase. 
#' @param classifiers Default = c('glmnet', 'svmLinear', 'svmRadial', 'kknn', 'rf'). Character vector indicating which of five classifiers to train and test. Users may select between one and all of the options. See \strong{Details}.
#' @param seed Optional integer seed for reproducibility. Seed is called each time an individual classifer is trained.
#' @param method Default = 'repeatedcv'. Character string specifying the resampling method to be used in \code{\link[caret]{trainControl}}. Options: 'boot', 'boot632', 'optimism_boot', 'boot_all', 'cv', 'repeatedcv', 'LOOCV', 'none' (only fits one model to the entire training set), 'oob' (only for random forest, bagged trees, bagged earth, bagged flexible discriminant analysis, or conditional tree forest models), timeslice, 'adaptive_cv', 'adaptive_boot' or 'adaptive_LGOCV'. See \href{http://topepo.github.io/caret/model-training-and-tuning.html#control}{http://topepo.github.io/caret/model-training-and-tuning.html#control}.
#' @param number Default = 10. Integer specifying the number or folders or number of resampling iterations to feed to \code{\link[caret]{trainControl}}. See \href{http://topepo.github.io/caret/model-training-and-tuning.html#control}{http://topepo.github.io/caret/model-training-and-tuning.html#control}.
#' @param repeats Default = 5. If using repeated k-fold cross-validation (\code{method = 'repeatedcv'}), the number of complete sets of folds to compute in \code{\link[caret]{trainControl}}. See \href{http://topepo.github.io/caret/model-training-and-tuning.html#control}{http://topepo.github.io/caret/model-training-and-tuning.html#control}.
#' @param search Default = 'grid'. Character string specifying how the tuning parameter grid is determined in  \code{\link[caret]{trainControl}}. See \href{http://topepo.github.io/caret/model-training-and-tuning.html#control}{http://topepo.github.io/caret/model-training-and-tuning.html#control}.
#' @param tuneGrids Default = NULL. If desired, users may input a list of tuning grids specific to each model, to be used within \code{\link[caret]{train}}. If using \code{tuneGrids}, list must be in same order as classifiers. Each list item should be a data.frame containing possible tuning values. The columns are named the same as the tuning parameters. Use \code{\link[caret:modelLookup]{getModelInfo}} to get a list of tuning parameters for each classifier model or see \href{http://topepo.github.io/caret/available-models.html}{http://topepo.github.io/caret/available-models.html}
#' @return \code{classifierModels} returns a nested list, where each classifier model has a unique name, and the model elements are stored in a list of nine. Model names for the classifiers are automatically generated as a string based on the templateID, score threshold, label type, and classifier name, with each element separated by an underscore (e.g. 'verd1_0.2_libraryID_glmnet').
#' 
#' Each classifier model's list contains nine elements defined as follows:
#' \enumerate{
#' \item{\strong{templateID}: the user-specified templateID.}
#' \item{\strong{label.type}: the user-specified label.type.}
#' \item{\strong{score.threshold}: the user-specified score.threshold.}
#' \item{\strong{train.scoreID}: vector of scoreIDs of verified scores used during the training phase.}
#' \item{\strong{test.scoreID}: vector of scoreIDs of verified scores used during the testing phase.}
#' \item{\strong{training.fit}: large list object generated by \strong{caret}, which stores all model training information (such as the 'trainingData', which are the features of each detected event). See \href{http://topepo.github.io/caret/index.html}{\strong{caret} documentation} for more information.}
#' \item{\strong{test.prediction}: a data.frame storing the predicted class (target signal, TS; or false alarm, FA) and target signal probabilities for each test event. }
#' \item{\strong{performance}: a data.frame storing this classifier's test phase performance on a variety of classifier assessment metrics, such as accuracy, sensitivity, specificity, precision, and F1 score.}
#' \item{\strong{confusion.matrix}: a table of the confusion matrix on the test data.} 
#'  }
#' @details 
#' 
#' The workflow in \code{classifierModels} depends heavily on the R package  \href{http://topepo.github.io/caret/index.html}{\strong{caret}}, which stands for \strong{c}lassiciation \strong{a}nd \strong{re}gression \strong{t}raining.
#' 
#'\code{classifierModels} splits verified scores from the selected \code{templateID}, \code{label.type}, and \code{score.threshold} into training and testing data according to the \code{split.proportion} argument, preserving any existing class imbalances between target signals and false alarms. Next, the training dataset is passed to a suite of statistical learning classifiers (also known as machine learning classifiers). The five available models are regularized logistic regression ('glmnet'), random forests ('rf'),  kernelized k-nearest neighbors ('kknn'), and two types of support vector machine, radial and linear ('svmRadial', 'svmLinear'). Each of the five classifiers is trained on the training data using a default of repeated 10-fold cross validation. All five classification algorithms fit models that map acoustic features (the features column of the scores table; see \strong{Details} of \code{\link{scoresDetect}}) to the labeled outputs (target signal or false alarm) in the verification data. Each classifier produces a probability that a given detection is a target signal.
#' 
#' Consult the \href{http://topepo.github.io/caret/index.html}{\strong{caret} documentation} for a wealth of additional details, particularly if venturing outside of the defaults for \code{method}, \code{number}, \code{repeats}, \code{search}, and \code{tuneGrids}. Note that not all \strong{caret} options are available in \code{classifierModels}, nor do we guarantee that all argument possibilities will work in \strong{AMMonitor} outside of the defaults. 
#' @seealso \href{placeholderurl}{AMMonitor Documentation Chapter 17: Classifications}, \code{\link{classifierPerformance}}, \code{\link{plotROC}}
#' @family classifier
#' @export
#' @examples
#' \dontrun{
#' 
#' # ------------------------------------------------------------
#' # Set up a demo AMMonitor database
#' # ------------------------------------------------------------
#' 
#' # Create the sample database containing necessary pre-populated data:
#' db.name <- 'demo.sqlite'
#' dbCreateSample(db.name = db.name, 
#'                file.path = paste0(getwd(),'/'), 
#'                tables = c('accounts', 'lists', 
#'                           'people', 'species',
#'                           'equipment', 'locations', 
#'                           'library', 'listItems',
#'                           'recordings', 'templates', 
#'                            'scores'))
#'                 
#' # Verify that the database exists in your current working directory:
#' file.exists(db.name)
#' 
#' # Connect to the db: 
#' db.path <- paste0(getwd(), '/', db.name)
#' conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
#'
#' # ---------------------------------------------------------------
#' # Add in verifications at the libraryID and speciesID level.
#' # Below, we manually input verification labels. (Note that in 
#' # practice, users should acquire verifications using scoresVerify(). 
#' # We update the verification labels manually here to easily 
#' # demonstrate the function.)
#' #----------------------------------------------------------------
#' 
#' # Update verifications in the manualVerifyLibraryID column
#' # of the scores table:
#' lib.verifications <- c(0,0,0,0,0,1,0,1,1,1,0,1,0,1,0,1,1,1,1,1,1,1)
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyLibraryID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = lib.verifications, scoreID = scoreIDs))
#'                      
#' # Update verifications in the manualVerifySpeciesID column
#' # of the scores table (remember, these will be different from the libraryID)
#' sp.verifications <- rep(1, length(scoreIDs))
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyspeciesID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = sp.verifications, scoreID = scoreIDs))
#'
#'           
#' #----------------------------------------------------------------
#' # Run classifierModels (perform data splitting, training, and 
#' # testing) and return items as a list. Use the templateID, 
#' # label.type, and score.threshold arguments to identify which
#' # labeled data should be used for training. 
#' #----------------------------------------------------------------
#' 
#' classifier.list <- classifierModels(db.path = db.path, 
#'                        templateID = 'verd1', 
#'                        label.type = 'libraryID', 
#'                        score.threshold = 0.2, 
#'                        split.proportion = 0.7, 
#'                        classifiers =  c('glmnet', 'svmLinear', 'rf'), 
#'                        seed = 3)
#'                        
#' # View information in the classifier.list object about the nine 
#' # returned elements for the first classifier:
#' str(classifier.list[[1]], max.level = 1)
#' 
#' #----------------------------------------------------------------
#' # Run classifierModels (perform data splitting, training, and 
#' # testing) and return items as a list. This time, use the 
#' # scoreID argument instead of the templateID and score.threshold
#' # arguments to directly identify which labeled data should be used for 
#' # training. The label.type argument is still required. 
#' #----------------------------------------------------------------
#' 
#' classifier.list <- classifierModels(db.path = db.path, 
#'                        scoreID = c(2,3,4,5,6,7,8,9,10,11,12,13,14,
#'                                    15,33,34,35,36,37,38), 
#'                        label.type = 'libraryID', 
#'                        split.proportion = 0.7, 
#'                        classifiers =  'rf', 
#'                        seed = 3)
#'                        
#' # View information in the classifier.list object about the nine 
#' # returned elements for the first classifier:
#' str(classifier.list[[1]], max.level = 1)
#'                  
#' # ------------------------------------------------------------
#' # Disconnect from demo database and delete from working directory
#' # ------------------------------------------------------------
#'
#' dbDisconnect(conx)
#' unlink(db.path)
#'         
#'}
#'

classifierModels <- function(db.path, 
                             templateID,
                             label.type,
                             score.threshold,
                             scoreID,
                             split.proportion = 0.7,
                             classifiers = c('glmnet', 'svmLinear', 'svmRadial', 'rf', 'kknn'),
                             seed,
                             method = 'repeatedcv',
                             number = 10,
                             repeats = 5,
                             search = 'grid',
                             tuneGrids = NULL
){
  
  # Perform checks
  if (missing(label.type) | !(label.type == 'speciesID' | label.type == 'libraryID')) 
    stop("Please use the label.type argument to specify how you would like to train these classifiers. \n  Train based on the call type ('libraryID') or the species ('speciesID').")
  
  if (!missing(templateID) & missing(score.threshold) & missing(scoreID)) 
    stop("If you are using the templateID argument, please also input a score.threshold.")
  
  if (missing(templateID) & !missing(score.threshold) & missing(scoreID)) 
    stop("If you are using the score.threshold argument, please also input a templateID.")
  
  if ((!missing(scoreID) & !missing(templateID)) | (!missing(scoreID) & !missing(score.threshold)))
    stop("If you are using the scoreID argument to directly select which labeled data you would like to use, you do not need to input a templateID or score.threshold.")
  
  # Assign correct label.type
  ifelse(label.type == 'speciesID', type <- 'manualVerifySpeciesID', type <- 'manualVerifyLibraryID')
  
  # Connect to the database:
  conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
  
  # Turn the SQLite foreign constraints on
  rs <- dbSendQuery(conn = conx, statement = "PRAGMA foreign_keys = ON;" )
  dbClearResult(rs)
  
  # Acquire scores for the desired label.type based on input templateID and score.threshold
  if (!missing(templateID) & !missing(score.threshold)) {
    scrs <- data.table(
      dbGetQuery(conn = conx,
                 statement = "SELECT scoreID, templateID, scoreThreshold, score, 
                 manualVerifyLibraryID, manualVerifySpeciesID, features 
                 FROM scores 
                 WHERE templateID = $templateID
                 AND scoreThreshold = $scoreThreshold",
                 params = list(templateID = templateID,
                               scoreThreshold = score.threshold)))
  }
  
  # Subset by scoreID directly, if desired by the user
  if (!missing(scoreID)) {
    scrs <- data.table(
      dbGetQuery(conn = conx,
                 statement = "SELECT scoreID, templateID, scoreThreshold, score, 
                 manualVerifyLibraryID, manualVerifySpeciesID, features 
                 FROM scores 
                 WHERE scoreID = $scoreID",
                 params = list(scoreID = scoreID)))
    
    # Check that all the same templateID & score.threshold, stop if not
    if (length(unique(scrs$templateID)) > 1)
      stop('The scoreID vector you have entered contains scoreIDs that correspond to more than one templateID. This is not allowed. See ?classifierModels and view the scoreID argument explanation.')
    if (length(unique(scrs$scoreThreshold)) > 1)
      stop('The scoreID vector you have entered contains scoreIDs that correspond to more than one score.threshold. This is not allowed. See ?classifierModels and view the scoreID argument explanation.')
    
    # If these are fine, then assign the templateID and score.threshold
    # based on the scoreID inputs:
    templateID <- unique(scrs$templateID)
    score.threshold <- unique(scrs$scoreThreshold)
  }
  
  # Subset the scores table so that it only contains rows with verifications
  # for the specified label.type
  scrs <- scrs[!(is.na(get(type))),]
  
  if (nrow(scrs) == 0 ) 
    return(message('No verifications for this template, label.type, and score.threshold.'))
  
  # Grab info associated with these scores
  vers <- scrs[,get(type)] # verification 0 1 data
  scoreIDs <- scrs[,scoreID] # scoreIDs
  
  # Split data into training and testing
  if (missing(seed)) seed <- NULL
  if (!(missing(seed))) set.seed(seed)
  training.index <- createDataPartition(y = vers, p = split.proportion, 
                                        list = FALSE, times = 1)
  train.scoreID <- scoreIDs[training.index]
  test.scoreID <- scoreIDs[-training.index]
  
  # Train classifiers
  fits <- classifierTrain(db.path = db.path,
                          train.scoreID = train.scoreID,
                          label.type = label.type,
                          classifiers = classifiers,
                          seed = seed,
                          tuneGrids = tuneGrids)
  
  # Test classifiers
  preds <- classifierTest(db.path = db.path,
                          test.scoreID = test.scoreID,
                          training.fits = fits)
  
  # Assess classifiers
  ens <- classifierAssess(test.predictions = preds)
  
  # Create raw classifier ammodel objects:
  classifier.list <- list()
  for (i in 1:length(fits$training.fits)) {
    training.time <- fits$training.time[[i]]
    classifier.list[[i]] <- list(templateID = templateID,
                                 label.type = label.type, 
                                 score.threshold = score.threshold, 
                                 train.scoreID = train.scoreID,
                                 test.scoreID = test.scoreID, 
                                 training.fit = fits$training.fits[[i]], # has the x.train data
                                 test.prediction = preds$test.predictions[[i]], 
                                 performance = ens$summary[names(fits$training.fits)[i], ],
                                 confusion.matrix = ens$confusion.matrices[names(fits$training.fits)[i]][[1]])
    
  }
  
  # Name the classifier models
  names(classifier.list) <- paste(templateID, score.threshold, label.type, 
                                  names(fits$training.fits), sep = '_')
  
  # Return classifier.list
  return(classifier.list)
}




#' @name classifierPerformance
#' @title Summarize classifier performance
#' @description Summarize how each classifier performed on the test data, using either an \code{amml} or \code{model.list} returned by \code{\link{classifierModels}}. 
#' @param amml Default = NULL. An \code{amModelLib} object. Not necessary if using the \code{model.list} argument.
#' @param model.list Default = NULL. A list of classification models returned by \code{\link{classifierModels}}. Not necessary if using the \code{amml} argument. 
#' @param model.names Default = NULL. Character vector of the model names for which to return performance metrics. If left NULL, function summarizes performance for all models in the \code{amml} or \code{model.list} object. 
#' @details In the returned data.table of performance metrics, we pay special attention to five key metrics. The most intuitive of these is  "accuracy", which represents the overall number of prediction cases correctly identified as target signals and false alarms by the classifier. 
#' 
#' However, accuracy is often a poor metric to rely upon when there is a class imbalance in the prediction data. In an example with 100 new prediction cases, if 95 events are false alarms, and 5 are target signals, a classifier that predicts everything to be a False Alarm will be 95\% accurate. If our real goal is to successfully find the five target signals, this classifier is useless to us despite its high accuracy. This phenomenon is known as the \href{https://en.wikipedia.org/wiki/Accuracy_paradox}{Accuracy Paradox}.
#' 
#' Instead of accuracy, we may strive to maximize a metric called \href{https://en.wikipedia.org/wiki/Sensitivity_and_specificity}{sensitivity} (also known as "recall" or "true positive rate"), which calculates the proportion of target signals correctly identified by the classifier. 
#' 
#' Conversely, we are often also interested in a classifier's \href{https://en.wikipedia.org/wiki/Sensitivity_and_specificity}{specificity} (also known as "true negative rate"), which speaks to its ability to correctly identify false alarms and label them as such. 
#' 
#' The interplay between specificity and sensitivity is captured by the metric \href{https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values}{precision} (also known as "positive predictive value"), which reflects the proportion of predicted target signals that are actually target signals. 
#' 
#' Finally, the \href{https://en.wikipedia.org/wiki/F1_score}{F1 score} represents a weighted average of precision and sensitivity, quantifying the tradeoff between a desire for high sensitivity and high precision.
#' 
#' For additional information on the columns in the returned data.table, see \href{http://topepo.github.io/caret/index.html}{\strong{caret} documentation}.

#' @return Data.table of 19 columns summarizing classifier performance. Columns include the model name (Model),  Accuracy, Kappa, AccuracyLower, AccuracyUpper, AccuracyNull, AccuracyPValue, McnemarPValue, Sensitivity, Specificity, Pos Pred Value, Neg Pred Value, Precision, Recall, F1, Prevalence, Detection Rate, Detection Prevalence, and Balanced Accuracy. 
#' @seealso \href{placeholderurl}{AMMonitor Documentation Chapter 17: Classifications}, \code{\link{classifierModels}}, \code{\link{plotROC}}
#' @family classifier
#' @export
#' @examples
#' \dontrun{
#'
#' #----------------------------------------------------------------
#' # Assess performance for classifier models in an amml:
#' #----------------------------------------------------------------
#'
#' # Read in a sample amml object:
#' data(classifiers_amml)
#'
#' # Assess performance for selected models in an amml:
#' classifierPerformance(amml = classifiers_amml,
#'                       model.list = NULL,
#'                       model.names = c("verd1_0.2_libraryID_glmnet",
#'                                       "verd1_0.2_libraryID_svmRadial"))
#'
#' # Assess performance for all classifier models in an amml:
#' classifierPerformance(amml = classifiers_amml,
#'                       model.list = NULL,
#'                       model.names = NULL)
#'
#' #----------------------------------------------------------------
#' # Assess performance for classifier models in a model.list
#' # generated by classifierModels():
#' #----------------------------------------------------------------
#'
#' # Read in a model list object produced by classifierModels():
#' data(classifier_practice)
#'
#' # Assess performance for selected models in a model.list:
#' classifierPerformance(amml = NULL,
#'                       model.list = classifier_practice,
#'                       model.names = c("verd1_0.2_libraryID_kknn",
#'                                       "verd1_0.2_libraryID_rf"))
#'
#' # Assess performance for all classifier models in a model.list:
#' classifierPerformance(amml = NULL, 
#'                       model.list = classifier_practice,
#'                       model.names = NULL)
#'         
#' }
#'

classifierPerformance <- function(amml = NULL,
                                  model.list = NULL,
                                  model.names = NULL)
{
  
  # Check amml & model.list inputs
  if (!(is.null(amml)) & !(is.null(model.list)))
    stop("Please use either the 'amml' or 'model.list' argument; not both.")
  if (is.null(amml) & is.null(model.list))
    stop("Please input either an 'amml' or 'model.list' argument; these arguments may not both be NULL.")
  
  # Check and inform user whether the model.names they've input are legal:
  if (!(is.null(model.names))) {
    if (!(is.null(amml))) mod.names <- names(amml@models)
    if (!(is.null(model.list))) mod.names <- names(model.list)
    
    dont.have <- model.names[!(model.names %in% mod.names)]
    if (length(dont.have) > 0 ) {
      stop('The following model(s) are not found: "', capture.output(dont.have), '".\n Check the spelling of inputs to model.names, and check that model(s) are present in the amml or model.list input you have provided.')
    }
  }
  
  # Create a mods list for storing models
  mods <- list()
  
  # If input object is an amml
  if (!(is.null(amml))) {
    
    # Get the names of all models in this amml
    if (is.null(model.names)) {model.names <- names(modelMeta(amml))}
    
    for (i in 1:length(model.names)) {
      mods[[i]] <- getAMModel(amml = amml, x = model.names[i])
    }
  }
  
  # If input object is merely a model list used for testing
  if (!(is.null(model.list))) {
    if (is.null(model.names)) model.names <- names(model.list)
    for (i in 1:length(model.names)) {
      # Into the mods object, extract only the model.names requesting by the user
      mods[[i]] <- model.list[[model.names[i]]]
    }
  }
  
  # Restructure each model's performance into a single table
  performance <- rbindlist(l = lapply(mods, '[[', 'performance'))
  performance <- cbind(model.names, performance)
  colnames(performance)[1] <- 'Model'
  
  # View performance summary 
  performance
}





#' @name classifierPredict
#' @title Use classifiers to predict on new detected events
#' @description Use calibrated classification models generated by \code{\link{classifierModels}} to predict the target signal probability and the class (TS: Target Signal, or FA: False Alarm) of new, incoming detected events produced by \code{\link{scoresDetect}}.
#' @param db.path The file path that houses the SQLite database.
#' @param amml An \code{amModelLib} object containing calibrated classification models generated by \code{\link{classifierModels}}. If following the recommended AMMonitor program protocol, this model library should be stored under the 'ammls' directory of your project. The library can be read into R as an object using \code{readRDS('ammls/classifiers.RDS')}.
#' @param date.range Length-two character vector of date ranges (inclusive) over which to run predictions. Dates should be given in YYYY-mm-dd format: e.g. c('2016-03-04', '2016-03-12'). 
#' @param templateID Character string of the templateID for which to make predictions. 
#' @param label.type Character string indicating whether to make predictions using the speciesID or the libraryID. Options = c('speciesID', 'libraryID').
#' @param score.threshold Numeric score.threshold used for this templateID. This value should match a score.threshold you have already used for this templateID to populate the scores table. 
#' @param classifiers Character vector indicating which classifiers to use for prediction. Note that these options will only be available if they were used during the calibration phase for this templateID_score.threshold_label.type combination when using the \code{\link{classifierModels}} function. Options are: c('glmnet', 'svmRadial', 'svmLinear', 'rf', 'kknn').
#' @param db.insert Default = FALSE. Logical flag for whether to insert records generated by this function into the database. If TRUE, events are added to the database in the classifications table. If FALSE, events are returned as a data.table for examination by the user, but not added to the database.
#' @return Data.table of inputs that would be added to the classifications table of an AMMonitor database. Columns include scoreID (integer), amml (character), classifier (character), modelName (character), modelProbability (numeric) and timestamp (character). The modelProbability value reflects the target signal probability assigned by a classifier. 
#' @seealso \href{placeholderurl}{AMMonitor Documentation Chapter 17: Classifications}, \code{\link{classifierModels}}, \code{\link{classifierPerformance}}, \code{\link{plotROC}}
#' @family classifier
#' @export
#' @examples
#' \dontrun{
#' 
#' # ------------------------------------------------------------
#' # Set up a demo AMMonitor database
#' # ------------------------------------------------------------
#' 
#' # Create the database (this will be deleted):
#' db.name <- 'demo.sqlite'
#' 
#' dbCreateSample(db.name = db.name, 
#'                file.path = paste0(getwd(),'/'), 
#'                tables = c('recordings', 'templates', 'locations', 
#'                           'equipment', 'people', 'accounts', 
#'                           'library', 'species', 'scores'))
#'                 
#' # Verify that the database exists in your current working directory:
#' file.exists(db.name)
#' 
#' # Connect to the db: 
#' db.path <- paste0(getwd(), '/', db.name)
#' conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
#' 
#' #----------------------------------------------------------------
#' # Predict on new data
#' #----------------------------------------------------------------
#' 
#' # Read in sample classifier models amml:
#' data(classifiers_amml)
#' 
#' # Predict class of new detections:
#' new.preds <- classifierPredict(db.path = db.path, 
#'                                amml = classifiers_amml,
#'                                date.range = c('2016-03-01', 
#'                                               '2016-03-30'),
#'                                templateID = 'verd1', 
#'                                label.type = 'libraryID', 
#'                                score.threshold = 0.2, 
#'                                classifiers =  c('glmnet', 'svmRadial',  
#'                                                 'svmLinear', 'rf', 
#'                                                 'kknn'), 
#'                                db.insert = FALSE)
#'                                
#' # ------------------------------------------------------------
#' # Disconnect from demo database and delete from working directory
#' # ------------------------------------------------------------
#'
#' dbDisconnect(conx)
#' unlink(db.path)
#' 
#'}
#'

classifierPredict <- function(db.path,
                              amml, 
                              date.range, 
                              templateID,
                              label.type,
                              score.threshold,
                              classifiers,  
                              db.insert = FALSE
) {
  
  # Perform checks
  if (missing(label.type) | !(label.type == 'speciesID' | label.type == 'libraryID')) {
    stop("Please use the label.type argument to specify which models should be used for prediction. \n Predict based on the call type ('libraryID') or the species ('speciesID').")
  }
  
  # Assign correct label.type
  ifelse(label.type == 'speciesID', type <- 'manualVerifySpeciesID', type <- 'manualVerifyLibraryID')
  
  # Connect to the database:
  conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
  
  # Turn the SQLite foreign constraints on
  rs <- dbSendQuery(conn = conx, statement = "PRAGMA foreign_keys = ON;" )
  dbClearResult(rs)
  
  # Check classifications table to see if there are already any for the desired 
  # label.type, templateID, score.threshold, classifiers. 
  # We do this by checking the model names.
  predict.mods <- paste(templateID, score.threshold, label.type, classifiers, sep = '_')
  
  # Check amml object for these model.names
  amml.names <- names(amml@models)
  dont.have <- predict.mods[!(predict.mods %in% amml.names)]
  if (length(dont.have) > 0 ) {
    stop('The following model(s) are not found: "', capture.output(dont.have), '".\n You likely did not add these models to your amml as an input to the "classifiers" argument of classifierModels(), or they have been deleted from the amml. Alternatively, check spelling and inputs to the "templateID", "score.threshold", "label.type", and "classifiers" arguments of classifierPredict().')
  }
  
  # Acquire scores for desired label.type based on templateID, score.threshold, date.range
  scrs <- data.table(
    dbGetQuery(conn = conx,
               statement = "SELECT scores.*, recordings.startDate
               FROM (scores INNER JOIN recordings
               ON scores.recordingID = recordings.recordingID)
               WHERE startDate BETWEEN date($firstDate) AND date($secondDate)
               AND templateID = $templateID
               AND scoreThreshold = $scoreThreshold",
               params = list(templateID = templateID,
                             scoreThreshold = score.threshold,
                             firstDate = date.range[1],
                             secondDate = date.range[2])))
  
  # Now we check the classifications table for these model names and scoreIDs
  predict.for <- expand.grid(scoreID = scrs$scoreID, 
                             modelName = predict.mods,
                             stringsAsFactors = FALSE)
  setDT(predict.for)
  already.have <- data.table(
    dbGetQuery(conn = conx,
               statement = "SELECT *
               FROM classifications 
               WHERE scoreID = $scoreID
               AND modelName = $modelName",
               params = as.list(predict.for)))
  
  # Remove model/scoreIDs combos for which we have already made predictions: 
  if (nrow(already.have) > 0) {
    predict.for <- predict.for[!already.have, on = .(scoreID, modelName)]
  } 
  
  # If we already have these predictions, stop
  if (nrow(predict.for) == 0 & db.insert == TRUE) {
    return(message('Predictions have already been made for these scoreID, templateID, score.threshold, and classifier combinations and will not be overwritten. Manually delete these records from the classifications table if needed.'))
  }
  
  # Unserialize features associated with these scores
  unserialized.feats <- lapply(scrs$features, 'unserialize')
  
  # If any failed to unserialize, try again: 
  check.ser <- sapply(unserialized.feats, class)
  bad.indices <- which(check.ser == 'raw')
  if (length(bad.indices)) {
    for (i in bad.indices) {
      unserialized.feats[[i]] <- unserialize(unserialized.feats[[i]])
    }
  }
  
  # Create matrix of prediction features
  x.predict <- rbindlist(l = lapply(unserialized.feats, 'as.data.table'))
  
  # Loop through all models and make predictions
  classifs <- list()
  for (m in 1:length(predict.mods)) {
    
    m.name <- predict.mods[m]
    mod <- getAMModel(amml = amml, x = m.name)
    fit <- mod$training.fit
    x.train <- mod$training.fit$trainingData
    
    # Make predictions
    probs <- predict(fit, x.predict, type = 'prob')
    
    # Arrange all prTS predictions
    all.probs <- data.table(scoreID = scrs$scoreID, modelProbability = probs$TS)
    
    # Extract the classifier type
    which.classifier <- strsplit(x = m.name, split = '_', fixed = TRUE)[[1]][[4]]
    
    # Generate the classifications table
    classifs[[m]] <- data.table(scoreID = all.probs$scoreID,
                                amml = 'classifiers',  # ??? should we not even bother with this col?
                                classifier = which.classifier, 
                                modelName = m.name, 
                                modelProbability = all.probs$modelProbability, 
                                timestamp = as.character(Sys.time())) # convert to char or else gets coerced
  }
  
  all.classifs <- rbindlist(l = classifs)
  
  # Insert classifications into database
  if (db.insert == TRUE) {
    
    # Get rid of any predictions we already have in the db that may have been tacked on
    all.classifs <- all.classifs[!already.have, on = .(scoreID, modelName)]
    
    # Append new data to table:
    dbWriteTable(conn = conx, name = 'classifications', value = all.classifs,
                 row.names = FALSE, overwrite = FALSE,
                 append = TRUE, header = FALSE)
    
    # If adding new results directly to database: 
    message('Added new data to classifications table.')
    
  }
  return(all.classifs)
}




#' @name classifierTest
#' @title Test classifiers on a subset of acoustic event data
#' @description Internal function. In practice, please use \code{\link{classifierModels}} instead.
#' @param db.path The file path that houses the SQLite database.
#' @param test.scoreID Character vector of scoreIDs to use for testing.
#' @param training.fits List of fitted classifier models produced by \code{\link{classifierTrain}}.
#' @return  A list of five, containing the following: 
#' \enumerate{
#' \item{\strong{test.predictions}: list of data.frames for each classifier, each containing the Class (TS: Target Signal or FA: False Alarm) and class probabilities for each test event.}
#' \item{\strong{label.type}: the user-specified label.type.}
#' \item{\strong{test.scoreID}: vector of scoreIDs of verified scores used during the testing phase.}
#' \item{\strong{x.test}: data.table of acoustic features for each testing event.}
#' \item{\strong{y.test}: factor vector of classes (TS: Target Signal; FA: False Alarm) for each event used in testing.}
#' }
#' @seealso \code{\link{classifierModels}}
#' @export
#' @examples
#' \dontrun{
#' 
#' # ------------------------------------------------------------
#' # Set up a demo AMMonitor database
#' # ------------------------------------------------------------
#' 
#' # Create the database (this will be deleted):
#' db.name <- 'demo.sqlite'
#' 
#' dbCreateSample(db.name = db.name, 
#'                file.path = paste0(getwd(),'/'), 
#'                tables = c('recordings', 'templates', 'locations', 
#'                           'equipment', 'people', 'accounts', 
#'                           'library', 'species', 'scores'))
#'                 
#' # Verify that the database exists in your current working directory:
#' file.exists(db.name)
#' 
#' # Connect to the db: 
#' db.path <- paste0(getwd(), '/', db.name)
#' conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
#' 
#' # ---------------------------------------------------------------
#' # Add in verifications at the libraryID and speciesID level.
#' # Below, we manually input verification labels. (Note that in 
#' # practice, users should acquire verifications using scoresVerify(). 
#' # We update the verification labels manually here to easily 
#' # demonstrate the function.)
#' #----------------------------------------------------------------
#' 
#' # Update verifications in the manualVerifyLibraryID column
#' # of the scores table:
#' lib.verifications <- c(0,0,0,0,0,1,0,1,1,1,0,1,0,1,0,1,1,1,1,1,1,1)
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyLibraryID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = lib.verifications, scoreID = scoreIDs))
#'                      
#' # Update verifications in the manualVerifySpeciesID column
#' # of the scores table (remember, these will be different from the libraryID)
#' sp.verifications <- rep(1, length(scoreIDs))
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyspeciesID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = sp.verifications, scoreID = scoreIDs))
#'
#' #----------------------------------------------------------------
#' # Create Train/Test Data Partition
#' #----------------------------------------------------------------
#' 
#' # Create Data Partition (train/test)
#' set.seed(3)
#' training.index <- createDataPartition(y = verifications, p = 0.7, 
#'                                       list = FALSE, times = 1)
#' train.scoreID <- scoreIDs[training.index]
#' test.scoreID <- scoreIDs[-training.index]
#' 
#' #----------------------------------------------------------------
#' # Train the classifiers
#' #----------------------------------------------------------------
#' fits <- classifierTrain(db.path = db.path, 
#'                         train.scoreID = train.scoreID, 
#'                         label.type = 'libraryID',
#'                         classifiers = c('glmnet', 'svmLinear', 
#'                                         'svmRadial', 'rf', 'kknn'),
#'                         seed = 10)
#' 
#' #----------------------------------------------------------------
#' # Test the classifiers 
#' #----------------------------------------------------------------
#' test <- classifierTest(db.path = db.path, 
#'                        test.scoreID = test.scoreID, 
#'                        training.fits = fits)
#'                        
#' # ------------------------------------------------------------
#' # Disconnect from demo database and delete from working directory
#' # ------------------------------------------------------------
#'
#' dbDisconnect(conx)
#' unlink(db.path)
#' 
#' }
#'

classifierTest <- function(db.path, 
                           test.scoreID, 
                           training.fits
) {
  
  # Connect to the database:
  conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
  
  # Turn the SQLite foreign constraints on
  rs <- dbSendQuery(conn = conx, statement = "PRAGMA foreign_keys = ON;" )
  dbClearResult(rs)
  
  # Acquire scores based on input scoreIDs
  scrs <- data.table(
    dbGetQuery(conn = conx,
               statement = "SELECT * FROM scores WHERE scoreID = $scoreID",
               params = list(scoreID = test.scoreID)))
  
  # Disconnect from database
  dbDisconnect(conx)
  
  # Get the test labels:
  ifelse(training.fits$label.type == 'speciesID', 
         type <- 'manualVerifySpeciesID', type <- 'manualVerifyLibraryID')
  y.test <- scrs[,get(type)]
  y.test[y.test == 0] <- 'FA'
  y.test[y.test == 1] <- 'TS'
  y.test <- factor(y.test, levels = c('TS', 'FA'))
  
  # Unserialize features associated with these scores
  unserialized.feats <- lapply(scrs$features, 'unserialize')
  
  # If any failed to unserialize, try again: 
  check.ser <- sapply(unserialized.feats, class)
  bad.indices <- which(check.ser == 'raw')
  if (length(bad.indices)) {
    for (i in bad.indices) {
      unserialized.feats[[i]] <- unserialize(unserialized.feats[[i]])
    }
  }
  
  # Create matrix of testing features
  x.test <- rbindlist(l = lapply(unserialized.feats, 'as.data.table'))
  
  # Extract the raw fit objects
  raw.fits <- training.fits$training.fits
  
  # Generate and store target signal probabilities and class predictions:
  preds <- list()
  for (i in 1:length(raw.fits)) {
    Class <- factor(predict(raw.fits[[i]], x.test), levels = c('TS', 'FA'))
    probs <- predict(raw.fits[[i]], x.test, type = 'prob')
    preds[[i]] <- cbind(Class, probs)  
  }
  names(preds) <- names(raw.fits)
  return(list(test.predictions = preds, label.type = training.fits$label.type, 
              test.scoreID = test.scoreID, x.test = x.test, y.test = y.test))
}





#' @name classifierTrain
#' @title Train classifiers on a subset of acoustic event data
#' @description Internal function. In practice, please use \code{\link{classifierModels}} instead. 
#' @param db.path The file path that houses the SQLite database.
#' @param train.scoreID Character vector of scoreIDs upon which to train, generated by \code{\link[caret]{createDataPartition}}.
#' @param label.type Character string indicating whether to train classifiers using the speciesID or the libraryID. Options = c('speciesID', 'libraryID').
#' @param classifiers Default = c('glmnet', 'svmLinear', 'svmRadial', 'kknn', 'rf'). Character vector indicating which of five classifiers to train and test. Users may select between one and all of the options. See \strong{Details}.
#' @param seed Optional integer seed for reproducibility. Seed is called each time an individual classifer is trained and each time a data partition is created to split data into training and testing data. 
#' @param method Default = 'repeatedcv'. Character string specifying the resampling method to be used in \code{\link[caret]{trainControl}}. Options: 'boot', 'boot632', 'optimism_boot', 'boot_all', 'cv', 'repeatedcv', 'LOOCV', 'none' (only fits one model to the entire training set), 'oob' (only for random forest, bagged trees, bagged earth, bagged flexible discriminant analysis, or conditional tree forest models), timeslice, 'adaptive_cv', 'adaptive_boot' or 'adaptive_LGOCV'. See \href{http://topepo.github.io/caret/model-training-and-tuning.html#control}{http://topepo.github.io/caret/model-training-and-tuning.html#control}.
#' @param number Default = 10. Integer specifying the number or folders or number of resampling iterations to feed to \code{\link[caret]{trainControl}}. See \href{http://topepo.github.io/caret/model-training-and-tuning.html#control}{http://topepo.github.io/caret/model-training-and-tuning.html#control}.
#' @param repeats Default = 5. If using repeated k-fold cross-validation (\code{method = 'repeatedcv'}), the number of complete sets of folds to compute in \code{\link[caret]{trainControl}}. See \href{http://topepo.github.io/caret/model-training-and-tuning.html#control}{http://topepo.github.io/caret/model-training-and-tuning.html#control}.
#' @param search Default = 'grid'. Character string specifying how the tuning parameter grid is determined in  \code{\link[caret]{trainControl}}. See \href{http://topepo.github.io/caret/model-training-and-tuning.html#control}{http://topepo.github.io/caret/model-training-and-tuning.html#control}.
#' @param tuneGrids Default = NULL. If desired, users may input a list of tuning grids specific to each model, to be used within \code{\link[caret]{train}}. If using \code{tuneGrids}, list must be in same order as classifiers. Each list item should be a data.frame containing possible tuning values. The columns are named the same as the tuning parameters. Use \code{\link[caret:modelLookup]{getModelInfo}} to get a list of tuning parameters for each classifier model or see \href{http://topepo.github.io/caret/available-models.html}{http://topepo.github.io/caret/available-models.html}
#' @return A list of six, containing the following: 
#' \enumerate{
#' \item{\strong{training.fits}: large list object generated by \strong{caret}, which stores model training information for each model. See \href{http://topepo.github.io/caret/index.html}{\strong{caret} documentation} for more information.}
#' \item{\strong{training.time}: numeric vector storing the number of seconds it took to train each model.}
#' \item{\strong{label.type}: the user-specified label.type.}
#' \item{\strong{train.scoreID}: vector of scoreIDs of verified scores used during the training phase.}
#' \item{\strong{x.train}: data.table of acoustic features for each training event.}
#' \item{\strong{y.train}: factor vector of classes (TS: Target Signal; FA: False Alarm) for each event used in training.}
#' }
#' @seealso \code{\link{classifierModels}}
#' @export
#' @examples
#' \dontrun{
#' 
#' # ------------------------------------------------------------
#' # Set up a demo AMMonitor database
#' # ------------------------------------------------------------
#' 
#' # Create the database (this will be deleted):
#' db.name <- 'demo.sqlite'
#' 
#' dbCreateSample(db.name = db.name, 
#'                file.path = paste0(getwd(),'/'), 
#'                tables = c('recordings', 'templates', 'locations', 
#'                           'equipment', 'people', 'accounts', 
#'                           'library', 'species', 'scores'))
#'                 
#' # Verify that the database exists in your current working directory:
#' file.exists(db.name)
#' 
#' # Connect to the db: 
#' db.path <- paste0(getwd(), '/', db.name)
#' conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
#' 
#' # ---------------------------------------------------------------
#' # Add in verifications at the libraryID and speciesID level.
#' # Below, we manually input verification labels. (Note that in 
#' # practice, users should acquire verifications using scoresVerify(). 
#' # We update the verification labels manually here to easily 
#' # demonstrate the function.)
#' #----------------------------------------------------------------
#' 
#' # Update verifications in the manualVerifyLibraryID column
#' # of the scores table:
#' lib.verifications <- c(0,0,0,0,0,1,0,1,1,1,0,1,0,1,0,1,1,1,1,1,1,1)
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyLibraryID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = lib.verifications, scoreID = scoreIDs))
#'                      
#' # Update verifications in the manualVerifySpeciesID column
#' # of the scores table (remember, these will be different from the libraryID)
#' sp.verifications <- rep(1, length(scoreIDs))
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyspeciesID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = sp.verifications, scoreID = scoreIDs))
#'
#' #----------------------------------------------------------------
#' # Create Train/Test Data Partition
#' #----------------------------------------------------------------
#' 
#' # Create Data Partition (train/test)
#' set.seed(3)
#' training.index <- createDataPartition(y = verifications, p = 0.7, 
#'                                       list = FALSE, times = 1)
#' train.scoreID <- scoreIDs[training.index]
#' test.scoreID <- scoreIDs[-training.index]
#' 
#' #----------------------------------------------------------------
#' # Train the classifiers
#' #----------------------------------------------------------------
#' fits <- classifierTrain(db.path = db.path, 
#'                         train.scoreID = train.scoreID, 
#'                         label.type = 'libraryID',
#'                         classifiers = c('glmnet', 'svmLinear', 
#'                                         'svmRadial', 'rf', 'kknn'),
#'                         seed = 10)
#'                         
#' # ------------------------------------------------------------
#' # Disconnect from demo database and delete from working directory
#' # ------------------------------------------------------------
#'
#' dbDisconnect(conx)
#' unlink(db.path)
#' }
#'

classifierTrain <- function(db.path, 
                            train.scoreID,
                            label.type,
                            classifiers = c('glmnet', 'svmLinear', 'svmRadial', 'rf', 'kknn'),
                            seed,
                            method = 'repeatedcv',
                            number = 10,
                            repeats = 5,
                            search = 'grid',
                            tuneGrids = NULL,
                            preProcess = NULL,
                            verboseIter = TRUE
){
  
  if (missing(label.type) | !(label.type == 'speciesID' | label.type == 'libraryID')) {
    stop("Please use the label.type argument to specify how you would like to train these classifiers. \n  Train based on the call type ('libraryID') or the species ('speciesID').")
  }
  
  # Assign correct label.type
  ifelse(label.type == 'speciesID', type <- 'manualVerifySpeciesID', type <- 'manualVerifyLibraryID')
  
  # Connect to the database:
  conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
  
  # Turn the SQLite foreign constraints on
  rs <- dbSendQuery(conn = conx, statement = "PRAGMA foreign_keys = ON;" )
  dbClearResult(rs)
  
  # Acquire scores based on input train.scoreIDs
  scrs <- data.table(
    dbGetQuery(conn = conx,
               statement = "SELECT * FROM scores WHERE scoreID = $scoreID",
               params = list(scoreID = train.scoreID)))
  
  # Unserialize features associated with these scores
  unserialized.feats <- lapply(scrs$features, 'unserialize')
  
  # If any failed to unserialize, try again: 
  check.ser <- sapply(unserialized.feats, class)
  bad.indices <- which(check.ser == 'raw')
  if (length(bad.indices)) {
    for (i in bad.indices) {
      unserialized.feats[[i]] <- unserialize(unserialized.feats[[i]])
    }
  }
  
  # Create matrix of training features
  # x.train Numeric matrix of features, with rows as observations and columns as features
  x.train <- rbindlist(l = lapply(unserialized.feats, 'as.data.table'), use.names = TRUE)
  
  # Create factor vector of class labels
  # y.train Factor vector of class labels. These must be converted from
  # TRUE/FALSE into something akin to "TS"/"FA" in order to avoid issues 
  # with how R reads logical objects when training the classifiers.
  y.train <- scrs[,get(type)]
  y.train[y.train == 0] <- 'FA'
  y.train[y.train == 1] <- 'TS'
  y.train <- factor(y.train, levels = c('TS', 'FA'))# ordered = TRUE)
  
  fits <- list()
  training.time <- matrix(0, nrow = length(classifiers), ncol = 1)
  
  # Set caret trainControl
  trc <- trainControl(method = method,
                      number = number,
                      repeats = repeats, # only relevant for 'repeatedcv' method
                      summaryFunction = twoClassSummary,
                      search = search,
                      classProbs = TRUE,
                      savePredictions = TRUE,
                      # seeing if this helps reproduce
                      seeds = NA, #A value of NA will stop the seed from being set within the worker processes 
                      verboseIter = TRUE)
  
  # Loop through classifiers to train the models
  for (i in 1:length(classifiers)) {
    
    start <- Sys.time()
    
    if (classifiers[i] %in% c('svmLinear', 'svmRadial')) {
      # Use scale arg FALSE for SVM -- this is not actually an argument to train(),
      #  but an argument passed to the SVM functions that needs to be defined as FALSE
      #  here in order to avoid errors
      if (!(missing(seed))) set.seed(seed)
      fits[[i]] <- train(x.train, y.train,
                         method = classifiers[i],
                         trControl = trc,
                         scale = FALSE,
                         metric = 'ROC',
                         preProcess = NULL,
                         tuneGrid = tuneGrids[[i]])
    } else {
      # No scale arg needed for other classifier types
      if (!(missing(seed))) set.seed(seed)
      fits[[i]] <- train(x.train, y.train,
                         method = classifiers[i],
                         trControl = trc,
                         metric = 'ROC',
                         preProcess = NULL,
                         tuneGrid = tuneGrids[[i]])
    }
    dif <- Sys.time() - start
    units(dif) <- 'secs'
    
    # Store the time used for training
    training.time[i,1] <- round(dif,2)
  }
  
  names(fits) <- paste0(classifiers)
  rownames(training.time) <- classifiers
  list(training.fits = fits, training.time = training.time, 
       label.type = label.type, train.scoreID = train.scoreID, 
       x.train = x.train, y.train = y.train)
}





#' @name plotROC
#' @title Plot ROC or precision-recall curves for testing data
#' @description Plot receiver operating characteristic (ROC) curves or precision-recall curves, using either an \code{amml} or \code{model.list} returned by \code{\link{classifierModels}}. 
#' @param db.path The file path that houses the SQLite database.
#' @param amml Default = NULL. An \code{amModelLib} object. Not necessary if using the \code{model.list} argument.
#' @param model.list Default = NULL. A list of classification models returned by \code{\link{classifierModels}}. Not necessary if using the \code{amml} argument. 
#' @param model.names Default = NULL. Character vector of the model names to plot. If left NULL, function plots all models in the \code{amml} or \code{model.list} object. \strong{Note that \code{plotROC} is meant to be used for models generated from a specific templateID_score.threshold_label.type combination} in order to compare that specific combination's glmnet, svmLinear, svmRadial, rf, and kknn classifiers. Users should take care to ensure model.names input to this function all come from the same templateID_score.threshold_label.type combination. 
#' @param curve.type Default = 'roc'. Character string specifying the type of curve to plot; 'roc' for ROC curve, 'pr' for precision-recall curve. 
#' @param data.type Default = 'test'. Character string specifying which data should be used in the plot: 'train' to plot curves of results from the training phase, 'test' to plot curves of results from the testing phase.
#' @param ensembles Only allowed if data.type = 'test'. A character vector of the ensemble classifier results to plot, if plotting based on the test data. Options are: c('accuracy', 'sensitivity', 'specificity', 'precision', 'f1', 'simple')).
#' @param main Character string indicating the plot title, if desired. If not provided, this function plots "ROC Plot" or "Precision-Recall Plot" as the title. 
#' @param ... Additional arguments to pass to \code{\link[graphics]{plot}}.
#' @seealso \href{placeholderurl}{AMMonitor Documentation Chapter 17: Classifications}, \code{\link{classifierModels}}, \code{\link{classifierPerformance}}
#' @return Plot of ROC or Precision-Recall curves. Each classifier's performance is plotted with a different color in the legend. The area under the ROC curve (AUC ROC) or area under the precision-recall curve (AUC PR) is stated in brackets next to each classifier's name in the legend. 
#' @references 
#' Code for precision-recall curve calculations is from:  
#' \href{https://github.com/kboyd/raucpr/blob/master/precision_recall.r}{https://github.com/kboyd/raucpr/blob/master/precision_recall.r}
#' @export
#' @examples
#' \dontrun{
#' # ------------------------------------------------------------
#' # Set up a demo AMMonitor database
#' # ------------------------------------------------------------
#'
#' # Create the database (this will be deleted):
#' db.name <- 'demo.sqlite'
#' 
#' dbCreateSample(db.name = db.name,
#'                file.path = paste0(getwd(),'/'),
#'                tables = c('recordings', 'templates', 'locations',
#'                           'equipment', 'people', 'accounts',
#'                           'library', 'species', 'scores',
#'                           'classifications'))
#'
#' # Verify that the database exists in your current working directory:
#' file.exists(db.name)
#'
#' # Connect to the db:
#' db.path <- paste0(getwd(), '/', db.name)
#' conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
#'
#' # ---------------------------------------------------------------
#' # Add in verifications at the libraryID and speciesID level.
#' # Below, we manually input verification labels. (Note that in 
#' # practice, users should acquire verifications using scoresVerify(). 
#' # We update the verification labels manually here to easily 
#' # demonstrate the function.)
#' #----------------------------------------------------------------
#' 
#' # Update verifications in the manualVerifyLibraryID column
#' # of the scores table:
#' lib.verifications <- c(0,0,0,0,0,1,0,1,1,1,0,1,0,1,0,1,1,1,1,1,1,1)
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyLibraryID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = lib.verifications, scoreID = scoreIDs))
#'                      
#' # Update verifications in the manualVerifySpeciesID column
#' # of the scores table (remember, these will be different from the libraryID)
#' sp.verifications <- rep(1, length(scoreIDs))
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyspeciesID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = sp.verifications, scoreID = scoreIDs))
#'
#'
#' #----------------------------------------------------------------
#' # plotROC() for classifier models in an amml:
#' #----------------------------------------------------------------
#',
#' # Read in a sample amml object:
#' data(classifiers_amml)
#'
#' # Plot ROC curve for selected classifiers in an amml, using the 
#' # testing data: 
#' plotROC(db.path = db.path, 
#'         amml = classifiers_amml,
#'         model.list = NULL,
#'         model.names = c("verd1_0.2_libraryID_glmnet",
#'                         "verd1_0.2_libraryID_svmRadial",
#'                         "verd1_0.2_libraryID_rf",
#'                         "verd1_0.2_libraryID_kknn",
#'                         "verd1_0.2_libraryID_svmLinear"),
#'         curve.type = 'roc',
#'         data.type = 'test',
#'         ensembles = c('sensitivity', 'f1'))
#'
#' # Plot precision-recall curve for selected classifiers in an amml,
#' # using the training data:
#' plotROC(db.path = db.path, 
#'         amml = classifiers_amml,
#'         model.list = NULL,
#'         model.names = c("verd1_0.2_libraryID_glmnet",
#'                         "verd1_0.2_libraryID_svmRadial",
#'                         "verd1_0.2_libraryID_svmLinear",
#'                         "verd1_0.2_libraryID_kknn",
#'                         "verd1_0.2_libraryID_rf"),
#'         curve.type = 'pr',
#'         data.type = 'train')
#'
#' #----------------------------------------------------------------
#' # plotROC() for classifier models in a model.list
#' # generated by classifierModels():
#' #----------------------------------------------------------------
#'
#' # Read in a model list object produced by classifierModels():
#' data(classifier_practice)
#'
#' # Plot ROC curve for selected classifiers in a model.list,
#' # using the training data
#' plotROC(db.path = db.path, 
#'         amml = NULL,
#'         model.list = classifier_practice,
#'         model.names = c("verd1_0.2_libraryID_rf",
#'                         "verd1_0.2_libraryID_kknn"),
#'         curve.type = 'roc',
#'         data.type = 'train')
#'
#' # Plot precision-recall curve for all classifiers in a model.list,
#' # using the testing data
#' plotROC(db.path = db.path, 
#'         amml = NULL,
#'         model.list = classifier_practice,
#'         model.names = c("verd1_0.2_libraryID_glmnet",
#'                         "verd1_0.2_libraryID_svmRadial",
#'                         "verd1_0.2_libraryID_svmLinear",
#'                         "verd1_0.2_libraryID_kknn",
#'                         "verd1_0.2_libraryID_rf"),
#'         curve.type = 'pr',
#'         data.type = 'test',
#'         ensembles = c('simple', 'accuracy', 'sensitivity'))
#'         
#' # ------------------------------------------------------------
#' # Disconnect from demo database and delete from working directory
#' # ------------------------------------------------------------
#'
#' dbDisconnect(conx)
#' unlink(db.path)
#'         
#' }
#'


plotROC <- function(db.path, 
                    amml = NULL,
                    model.list = NULL,
                    model.names = NULL,
                    curve.type = 'roc',
                    data.type = 'test',
                    ensembles = NULL, 
                    main = NULL,
                    ...)
{
  
  # Check inputs
  if (missing(db.path)) stop("Please input 'db.path'.")
  if (!(is.null(amml)) & !(is.null(model.list)))
    stop("Please use either the 'amml' or 'model.list' argument; not both.")
  if (is.null(amml) & is.null(model.list))
    stop("Please input either an 'amml' or 'model.list' argument; these arguments may not both be NULL.")
  if (data.type == 'train' & !is.null(ensembles))
    stop("Ensemble curves may only be computed if using data.type = 'test'.")
  
  # Check and inform user whether the model.names they've input are legal:
  if (!(is.null(model.names))) {
    if (!(is.null(amml))) mod.names <- names(amml@models)
    if (!(is.null(model.list))) mod.names <- names(model.list)
    
    dont.have <- model.names[!(model.names %in% mod.names)]
    if (length(dont.have) > 0 ) {
      stop('The following model(s) are not found: "', capture.output(dont.have), '".\n Check the spelling of inputs to model.names, and check that model(s) are present in the amml or model.list input you have provided.')
    }
  }
  
  # Trapezoid rule calculation from pracma so that we don't have to import the package.
  # http://svitsrv25.epfl.ch/R-doc/library/caTools/html/trapz.html
  #  This is a form of composite trapezoid rule as detailed here:
  #  http://www.aaronschlegel.com/the-trapezoidal-rule-of-numerical-integration-in-r/
  trapezoid.rule <- function(x,y){
    m <- length(x)
    xp <- c(x, x[m:1])
    yp <- c(numeric(m), y[m:1])
    n <- 2 * m
    p1 <- sum(xp[1:(n - 1)] * yp[2:n]) + xp[n] * yp[1]
    p2 <- sum(xp[2:n] * yp[1:(n - 1)]) + xp[1] * yp[n]
    return(0.5 * (p1 - p2))
  }
  
  # Set up list holding tank
  mods <- list()
  
  # If input object is an amml
  if (!(is.null(amml))) {
    
    # Get the names of all models in this amml
    if (is.null(model.names)) model.names <- names(modelMeta(amml))
    
    for (i in 1:length(model.names)) {
      mods[[i]] <- getAMModel(amml = amml, x = model.names[i])
    }
  }
  
  # If input object is merely a model list used for testing
  if (!(is.null(model.list))) {
    if (is.null(model.names)) model.names <- names(model.list)
    for (i in 1:length(model.names)) {
      # Into the mods object, extract only the model.names requested by the user
      mods[[i]] <- model.list[[model.names[i]]]
    }
  }
  names(mods) <- model.names
  
  # Initialize some colors for plotting
  col <- c('red','gray', 'blue','orange','green','lightslateblue',
           'firebrick4','goldenrod3','lightseagreen','deeppink4','black')
  
  # Find the label type and the true classes (should be same for all models input)
  label.type <- strsplit(x = names(mods)[[1]], split = '_', fixed = TRUE)[[1]][[3]]
  
  # Get the labels for training or testing data
  if (data.type == 'test') scoreIDs <- mods[[1]]$test.scoreID
  if (data.type == 'train') scoreIDs <- mods[[1]]$train.scoreID
  
  # Connect to the database:
  conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
  
  # Turn the SQLite foreign constraints on
  rs <- dbSendQuery(conn = conx, statement = "PRAGMA foreign_keys = ON;" )
  dbClearResult(rs)
  
  # Acquire scores for the desired label.type based on input templateID and score.threshold
  scrs <- data.table(
    dbGetQuery(
      conn = conx,
      statement = "SELECT scoreID, manualVerifyLibraryID, manualVerifySpeciesID
      FROM scores
      WHERE scoreID = $scoreID",
      params = list(scoreID = scoreIDs)))
  
  # Acquire classifications
  classif <- data.table(
    dbGetQuery(
      conn = conx,
      statement = 'SELECT * FROM classifications
      WHERE scoreID = $scoreID',
      params = list(scoreID = scoreIDs)
    )
  )
  
  # Find the true class labels for this label type, clean-up to be TS/FA labels
  ifelse(label.type == 'libraryID', colm <- 'manualVerifyLibraryID', colm <- 'manualVerifySpeciesID')
  if (data.type == 'test') {
    true.class <- scrs[,get(colm)]
    true.class[true.class == 0] <- 'FA'; true.class[true.class == 1] <- 'TS'
    true.class <- factor(x = true.class, levels = c('TS', 'FA'))
  }
  
  if (data.type == 'test') N <- length(mods) + length(ensembles)
  if (data.type == 'train') N <- length(mods)
  
  if (!is.null(ensembles)) {
    # Loop through ensembles
    ens.list <- list()
    for (e in 1:length(ensembles)) {
      ens.list[[e]] <- classifierEnsemble(db.path = db.path, 
                                          amml = amml,
                                          model.list = model.list,
                                          classificationsDT = classif,
                                          ensemble = ensembles[e])
    }
    names(ens.list) <- ensembles
  }
  
  aucs <- aucprs <- c()
  for (i in 1:N) {
    if (data.type == 'train') {
      # Gather probabilities of being a target signal
      probs <- mods[[i]]$training.fit$pred$TS
      
      # Gather true classes
      true.class <- mods[[i]]$training.fit$pred$obs
      
      # Order labels based on their TS probabilities
      actual.ord <- true.class[order(-probs)]
    }
    
    if (data.type == 'test') {
      # Order probabilities and classes
      if (i <= length(model.names)) {
        probs <- mods[[i]]$test.prediction$TS 
        actual.ord <- true.class[order(-probs)]
      } else {
        ens.list.i <- ens.list[[i - length(model.names)]]
        probs <- ens.list.i$modelProbability
        ens.list.i[,true.class := factor(NA)][
          class == 1, true.class := 'TS'][
            class == 0, true.class := 'FA']
        actual.ord <- ens.list.i$true.class[order(-probs)]
      }
    }
    
    # Sort probabilities
    sort.probs <- probs[order(-probs)]
    
    if (curve.type == 'roc') {
      pr.t <- sort.probs[actual.ord == 'TS']
      pr.f <- sort.probs[actual.ord == 'FA']
      
      tpr <- fpr <- 0
      for (p in 1:length(sort.probs)) {
        # Construct tpr & fpr by taking the mean of all probability values exceeding
        #  the value of the current step:
        tpr <- c(tpr, mean(pr.t >= sort.probs[p]))
        fpr <- c(fpr, mean(pr.f >= sort.probs[p]))
      }
      
      # Calculate AUC:
      # From: http://blog.revolutionanalytics.com/2016/11/calculating-auc.html
      dfpr <- c(diff(fpr), 0)
      dtpr <- c(diff(tpr), 0)
      aucs[i] <- sum(tpr * dfpr) + sum(dtpr*dfpr)/2
    }
    
    if (curve.type == 'pr') {
      
      # Calculate PR using the trapezoid rule
      if (data.type == 'train') {
        prrec <- pr(predictions = probs, labels = true.class)
      }
      if (data.type == 'test') {
        ifelse(i <= length(model.names), 
               lab <- true.class, 
               lab <- ens.list.i$true.class)
        prrec <- pr(predictions = probs, labels = lab)
      }
      
      aucprs[i] <- trapezoid.rule(prrec[,'Recall'], prrec[,'Precision'])
    }
    
    # Plot the curves
    if (i == 1) {
      if (curve.type == 'roc') {
        if (is.null(main)) main <- 'ROC Plot'
        plot(fpr, tpr, type = 'l', bty = 'l', col = col[i],
             xlab = 'False Positive Rate',
             ylab = 'True Positive Rate',...)
        title(main = main, line = 0)
        abline(0,1, lty = 2)
      } else {
        if (curve.type == 'pr') {
          if (is.null(main)) main <- 'Precision-Recall Plot'
          boundary.ratio <- sum(true.class == 'TS')/length(true.class)
          plot(x = prrec[,'Recall'], y = prrec[,'Precision'], type = 'l',
               xlim = c(0, 1), ylim = c(0, 1), col = col[i],
               xlab = 'Recall', ylab = 'Precision', ...)
          title(main = main, line = 0.5)
          abline(boundary.ratio, 0, lty = 2)} }
    } else {
      if (curve.type == 'roc') {
        lines(fpr, tpr, col = col[i], ...)
      } else {
        if (curve.type == 'pr') {
          lines(x = prrec[,'Recall'], y = prrec[,'Precision'],
                col = col[i], ...)
        }
      }
    }
    
    # Get classifier model.names
    leg <- unlist(lapply(strsplit(x = names(mods), split = '_', fixed = TRUE), '[[', 4))
    if (!is.null(ensembles)) leg <- c(leg, paste0('wt.by.', ensembles))
  }
  
  #ifelse(curve.type == 'roc', pos <- 'bottomright', pos <- 'bottomright')
  ifelse(curve.type == 'roc', msr <- aucs, msr <- aucprs)
  legend('bottomright', legend = paste0(leg,' [', round(msr,2),']'),
         pch = 20,
         bty = 'n',
         col = col,
         cex = 0.90)
  
  invisible(list(aucs = aucs, aucprs = aucprs))
  
}




#' @name plotVerifications
#' @title Plot side-by-side spectrograms of verified data
#' @description Generate plots of all verified target signals and false alarms for a given template and score.threshold combination, at the level of either the libraryID or the speciesID. All target signals are plotted side-by-side in a single plot window. All false alarms are plotted side-by-side in a separate plot window. 
#' @param db.path The file path that houses the SQLite database.
#' @param templateID Single character object indicating the templateID of detections to plot. All verified detections associated with this templateID are automatically plotted.
#' @param score.threshold Numeric score.threshold used to generate template-based detections. 
#' @param label.type Character string indicating whether to plot at the level of the speciesID or the libraryID. Options are c('speciesID', 'libraryID'). See \strong{Details} of \code{\link{scoresVerify}}.
#' @param plot.scoreID Default = FALSE. If TRUE, plots the scoreID over the detected event in red. Useful for identifying and correcting labeling errors.
#' @param new.window Default = TRUE. Logical value for whether to use \code{dev.new} to produce new plot windows for target signals and false alarms.
#' @param spec.col Default = gray.3(). The colors used to plot verification spectrograms. Spectrogram colors are adjustable, and users may create their own gradients for display. See \strong{Details}.
#' @param box.lwd Default = 1. Integer value for box line thickness.
#' @param plot.all Default = FALSE. Logical flag for whether to plot all verifications randomly together (target signals and false alarms). If TRUE, target signals and false alarms are plotted in random order, with no delineating boxes. 
#' @return Two plots; one of all verified target signals for this template, and the other of all verified false alarms for this template. \code{plotVerifications} also returns the row indices of each plotted verification in the scores data.table. These can be used to backtrack into the database and check or change labels that appear to be incorrect.
#' @details For plotting verifications, a few spectrogram color options are provided via the R package monitoR, including gray.1, gray.2, gray.3, rainbow.1, and topo.1, all of which are based on existing R colors. 
#' 
#' To produce two separate plot windows for target signals and false alarms, \code{plotVerifications} calls \code{dev.new} within the function when new.window = TRUE. Plot windows may not work the same on all operating systems. 
#' @family classifier
#' @seealso \href{placeholderurl}{AMMonitor Documentation Chapter 17: Classifications}, \code{\link{plotVerificationsAvg}}, \code{\link{scoresVerify}}
#' @export
#' @examples
#'\dontrun{
#'
#' # ------------------------------------------------------------
#' # Set up a demo AMMonitor database
#' # ------------------------------------------------------------
#' 
#' # Create the database (this will be deleted):
#' db.name <- 'demo.sqlite'
#' 
#' dbCreateSample(db.name = db.name, 
#'                file.path = paste0(getwd(),'/'), 
#'                tables = c('recordings', 'templates', 'locations', 
#'                           'equipment', 'people', 'accounts', 
#'                           'library', 'species', 'scores'))
#'                 
#' # Verify that the database exists in your current working directory:
#' file.exists(db.name)
#' 
#' # Connect to the db: 
#' db.path <- paste0(getwd(), '/', db.name)
#' conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
#' 
#' #----------------------------------------------------------------
#' # Add in verifications at the libraryID and speciesID level.
#' # Below, we manually input verification labels. (Note that in 
#' # practice, users should acquire verifications using scoresVerify(). 
#' # We update the verification labels manually here to easily 
#' # demonstrate the function.)
#' #----------------------------------------------------------------
#' 
#' # Update verifications in the manualVerifyLibraryID column
#' # of the scores table:
#' lib.verifications <- c(0,0,0,0,0,1,0,1,1,1,0,1,0,1,0,1,1,1,1,1,1,1)
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyLibraryID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = lib.verifications, scoreID = scoreIDs))
#'                      
#' # Update verifications in the manualVerifySpeciesID column
#' # of the scores table (remember, these will be different from the libraryID)
#' sp.verifications <- rep(1, length(scoreIDs))
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyspeciesID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = sp.verifications, scoreID = scoreIDs))
#'           
#' #----------------------------------------------------------------
#' # Plot verifications at the libraryID level in grayscale. 
#' # Set plot.scoreID = TRUE to plot scoreID values over each event. 
#' #----------------------------------------------------------------
#' plotVerifications(db.path = db.path, 
#'                   templateID = 'verd1', 
#'                   score.threshold = 0.2, 
#'                   label.type = 'libraryID', 
#'                   spec.col = gray.2(),
#'                   new.window = TRUE,
#'                   plot.scoreID = TRUE)
#'
#' #----------------------------------------------------------------
#' # Plot verifications at the species level in topo colors. 
#' # Set plot.scoreID = FALSE so that no scoreIDs are plotted. 
#' # (Note that there are no false alarms at the speciesID level, so
#' # the second plot will be empty, and the function returns a message 
#' # stating: "No false alarm verifications to display. Plotted target
#' # signals only.")
#' #----------------------------------------------------------------  
#' plotVerifications(db.path = db.path, 
#'                   templateID = 'verd1', 
#'                   score.threshold = 0.2, 
#'                   label.type = 'speciesID', 
#'                   spec.col = topo.1(),
#'                   new.window = TRUE, 
#'                   plot.scoreID = FALSE)
#'                   
#' # ------------------------------------------------------------
#' # Disconnect from demo database and delete from working directory
#' # ------------------------------------------------------------
#'
#' dbDisconnect(conx)
#' unlink(db.path)
#' 
#' }


plotVerifications <- function(db.path,
                              templateID,
                              score.threshold,
                              label.type,
                              plot.scoreID = FALSE,
                              new.window = TRUE, 
                              spec.col = gray.3(),
                              box.lwd = 1,
                              plot.all = FALSE)
{
  
  # This function currently makes a lot more sense for square templates only, since you 'verify'
  # based on a picture of a square (and not based within little squares)
  # So it might not be very useful for other types of templates in terms of interpretation 
  # Perform checks
  if (missing(label.type) | !(label.type == 'speciesID' | label.type == 'libraryID')) {
    stop("Please use the label.type argument to specify which verification type to plot. \n  Plot the call type ('libraryID') or the species ('speciesID').")
  }
  
  # Assign correct label.type
  ifelse(label.type == 'speciesID', type <- 'manualVerifySpeciesID', type <- 'manualVerifyLibraryID')
  
  # Connect to the database:
  conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
  
  # Turn the SQLite foreign constraints on
  rs <- dbSendQuery(conn = conx, statement = "PRAGMA foreign_keys = ON;" )
  dbClearResult(rs)
  
  # Acquire scores based on input scoreIDs and label.type
  qry <- paste("SELECT * FROM scores WHERE templateID = $templateID 
               AND scoreThreshold = $scoreThreshold AND", 
               type ,"is NOT NULL")
  scrs <- data.table(
    dbGetQuery(conn = conx,
               statement = qry,
               params = list(templateID = templateID,
                             scoreThreshold = score.threshold)))
  if (nrow(scrs) == 0) {
    stop('No verifications to plot for this templateID, score.threshold, and label.type.')
  }
  
  # Acquire and unserialize the template; access the underlying cor/binTemplate object
  template <- templatesUnserialize(db.path = db.path, templateID)@templates[[1]]  
  
  # Unserialize features associated with these scores
  feats <- lapply(scrs$features, 'unserialize')
  
  # If any failed to unserialize, try again: 
  check.ser <- sapply(feats, class)
  bad.indices <- which(check.ser == 'raw')
  if (length(bad.indices)) {
    for (i in bad.indices) {
      feats[[i]] <- unserialize(feats[[i]])
    }
  }
  
  # Find number of columnns and rows in this template
  ncols <- template@n.t.bins + 1
  nrows <- template@n.frq.bins + 1
  
  # Find indices of target signals and false alarms
  ts.inds <- which(scrs[,get(type)] == TRUE)
  fa.inds <- which(scrs[,get(type)] == FALSE)
  
  # Recreate spectrograms for each detection
  mats <- array(data = 0, dim = c(nrows, ncols, length(feats)))
  for (n in 1:length(feats)) {
    feat <- feats[[n]]
    amps <- as.numeric(feat[grep(pattern = 'amp', x = names(feat))])
    mat <- matrix(data = amps, nrow = nrows, ncol = ncols)
    mat[mat == 0] <- min(mat) - 10 # add color correction if encountering rect-select template
    mats[,,n] <- mat
  } # end for n
  
  # Find indices of target signals and false alarms; organize the matrices for plotting
  ts.mats <- mats[,,ts.inds]
  fa.mats <- mats[,,fa.inds]
  all.mats <- list(ts.mats = ts.mats, fa.mats = fa.mats)
  
  # If plotting target signals and false alarms in two separate windows:
  if (plot.all == FALSE) {
    for (i in 1:2) {
      if (new.window) dev.new()
      mats.i <- all.mats[[i]]
      
      ifelse(i == 1, 
             score.id.plot <- scrs[ts.inds, scoreID], 
             score.id.plot <- scrs[fa.inds, scoreID])
      
      if (length(mats.i) == 0 & i == 1) {
        message('No target signal verifications to display. Plotted false alarms only.')
        next
      }
      
      if (length(mats.i) == 0 & i == 2) {
        message('No false alarm verifications to display. Plotted target signals only.')
        return()
      }
      
      # Set up plotting dimensions
      sqrdim <- ceiling(sqrt(dim(mats.i)[3]))
      
      # If only one verification to plot, plot as 1x1
      if (is.na(sqrdim)) {
        dim1 <- dim2 <- 1 
      } else {
        # Else, calculate reasonable plotting dimensions 
        # for how to arrange the events in the plot window: 
        if (sqrdim^2 %% dim(mats.i)[3] > sqrdim) {
          dim1 <- sqrdim - 1; dim2 <- sqrdim
        } else {
          dim1 <- dim2 <- sqrdim
        }
      }
      
      # Set plot params according to these dimensions
      par(mar = rep(0,4), mfrow = c(dim1, dim2)) 
      
      # Loop through and add each verified event to the plot window
      for (pl in 1:length(score.id.plot)) {
        
        ifelse(length(score.id.plot) == 1, # deal with cases of only one verif.
               plot.this <- t(mats.i),
               plot.this <- t(mats.i[,,pl]))
        
        # Plot a single detected event
        image(x = 1:ncols, y = 1:nrows, plot.this,
              yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', col = spec.col)
        
        # Plot scoreIDs if desired
        if (plot.scoreID == TRUE) { 
          text(x = ncols/2, y = nrows/2, labels = score.id.plot[pl], col = 'red', cex = 2) 
        }
        ifelse(i == 1, box.col <- 'green3', box.col <- 'red')
        box(col = box.col, lwd = box.lwd)
      }
    }
    
  } else {
    # If plot.all == TRUE, then 
    # Plot TS and FA side by side in a jumble:
    if (new.window) dev.new()
    sqrdim <- ceiling(sqrt(dim(mats)[3]))
    if (sqrdim^2 %% dim(mats)[3] > sqrdim) {
      dim1 <- sqrdim - 1; dim2 <- sqrdim
    } else {
      dim1 <- dim2 <- sqrdim
    }
    
    # set pars
    par(mar = rep(0,4), mfrow = c(dim1, dim2)) 
    
    # Created a random sequence to plot them all in a jumble
    rand.seq <- sample(x = 1:dim(mats)[3], size = dim(mats)[3], replace = FALSE)
    for (pl in rand.seq) {
      image(x = 1:ncols, y = 1:nrows, t(mats[,,pl]),
            yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', col = spec.col)
    }
  }
  
  # If saved as an object, return score ids for target signals and false alarms
  #   (can be matched to the visual plot if plot.scoreID = TRUE, and later corrected in the database)
  ts.ids <- scrs[ts.inds, scoreID]
  fa.ids <- scrs[fa.inds, scoreID]
  score.ids <- list(target.signal.ids = ts.ids, false.alarm.ids = fa.ids)
  return(score.ids)
}





#' @name plotVerificationsAvg
#' @title Plot average template detections of verified events
#' @description Generate a four-panel plot that shows: 1. the template used, 2. a spectrogram showing the mean of all verified events (where each pixel in the spectrogram reflects the mean amplitude value at that pixel across all verified events), 3. a spectrogram of the mean target signal (where each pixel in the spectrogram reflects the mean amplitude value at that pixel across all events verified as target signals), and 4. a spectrogram of the mean false alarm (where each pixel in the spectrogram reflects the mean amplitude value at that pixel across all events verified as false alarms).
#' @param db.path The file path that houses the SQLite database.
#' @param templateID Single character object indicating the templateID of detections to plot. All verified detections associated with this templateID are automatically plotted.
#' @param score.threshold Numeric score.threshold used to generate detections. 
#' @param label.type Character string indicating whether to plot at the level of the speciesID or the libraryID. Options are c('speciesID', 'libraryID'). See \strong{Details} of \code{\link{scoresVerify}}.
#' @param spec.col Default = gray.2(). The colors used to plot verification spectrograms. Spectrogram colors are adjustable, and users may create their own gradients for display. See \strong{Details}.
#' @details For plotting verifications, a few spectrogram color options are provided via the R package monitoR, including gray.1, gray.2, gray.3, rainbow.1, and topo.1, all of which are based on existing R colors. 
#' @return Four-panel plot of the template, average verified detection, average verified target signal, and average verified false alarm.
#' @family classifier
#' @seealso \href{placeholderurl}{AMMonitor Documentation Chapter 17: Classifications}, \code{\link{plotVerifications}}, \code{\link{scoresVerify}}
#' @export
#' @examples
#' \dontrun{
#' 
#' # ------------------------------------------------------------
#' # Set up a demo AMMonitor database
#' # ------------------------------------------------------------
#' 
#' # Create the database (this will be deleted):
#' db.name <- 'demo.sqlite'
#' 
#' dbCreateSample(db.name = db.name, 
#'                file.path = paste0(getwd(),'/'), 
#'                tables = c('recordings', 'templates', 'locations', 
#'                           'equipment', 'people', 'accounts', 
#'                           'library', 'species', 'scores'))
#'                 
#' # Verify that the database exists in your current working directory:
#' file.exists(db.name)
#' 
#' # Connect to the db: 
#' db.path <- paste0(getwd(), '/', db.name)
#' conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
#' 
#' #----------------------------------------------------------------
#' # Add in verifications at the libraryID and speciesID level.
#' # Below, we manually input verification labels. (Note that in 
#' # practice, users should acquire verifications using scoresVerify(). 
#' # We update the verification labels manually here to easily 
#' # demonstrate the function.)
#' #----------------------------------------------------------------
#' 
#' # Update verifications in the manualVerifyLibraryID column
#' # of the scores table:
#' lib.verifications <- c(0,0,0,0,0,1,0,1,1,1,0,1,0,1,0,1,1,1,1,1,1,1)
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyLibraryID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = lib.verifications, scoreID = scoreIDs))
#'                      
#' # Update verifications in the manualVerifySpeciesID column
#' # of the scores table (remember, these will be different from the libraryID)
#' sp.verifications <- rep(1, length(scoreIDs))
#' scoreIDs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,33,34,35,36,37,38,39)
#' dbExecute(conn = conx, statement =  "UPDATE scores 
#'           SET manualVerifyspeciesID = $vers WHERE scoreID = $scoreID", 
#'           param = list(vers = sp.verifications, scoreID = scoreIDs))
#'
#' #----------------------------------------------------------------
#' # Plot average verified detections for the libraryID in grayscale
#' #----------------------------------------------------------------
#' # Specify a new plot window (may not work on non-Windows machines)
#' x11(height = 2, width = 7)   
#' 
#' # Plot average verifications   
#' plotVerificationsAvg(db.path = db.path, 
#'                      templateID = 'verd1',
#'                      score.threshold = 0.2,
#'                      label.type = 'libraryID',
#'                      spec.col = gray.2())
#'                   
#' #----------------------------------------------------------------
#' # Plot average verified detections for the speciesID in heat colors
#' # (note that there are no false alarms at the speciesID level for
#' # this example, so the fourth plot will be absent!)
#' #----------------------------------------------------------------     
#' # Specify a new plot window (may not work on non-Windows machines)
#' x11(height = 2, width = 7)        
#' 
#' # Plot average verifications
#' plotVerificationsAvg(db.path = db.path, 
#'                      templateID = 'verd1',
#'                      score.threshold = 0.2,
#'                      label.type = 'speciesID',
#'                      spec.col = rev(heat.colors(n = 100)))
#'                   
#' # ------------------------------------------------------------
#' # Disconnect from demo database and delete from working directory
#' # ------------------------------------------------------------
#'
#' dbDisconnect(conx)
#' unlink(db.path)
#' 
#' 
#' }

plotVerificationsAvg <- function(db.path,
                                 templateID,
                                 score.threshold,
                                 label.type,
                                 spec.col = gray.2())
{
  
  # This function currently makes a lot more sense for square templates only, since you 'verify'
  # based on a picture of a square (and not based within little squares)
  # So it might not be as useful for non-square templates in terms of interpretation 
  # Perform checks
  if (missing(label.type) | !(label.type == 'speciesID' | label.type == 'libraryID')) {
    stop("Please use the label.type argument to specify which verification type to plot. \n  Plot the call type ('libraryID') or the species ('speciesID').")
  }
  
  # Assign correct label.type
  ifelse(label.type == 'speciesID', type <- 'manualVerifySpeciesID', type <- 'manualVerifyLibraryID')
  
  # Connect to the database:
  conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
  
  # Turn the SQLite foreign constraints on
  rs <- dbSendQuery(conn = conx, statement = "PRAGMA foreign_keys = ON;" )
  dbClearResult(rs)
  
  # Acquire scores based on input scoreIDs and label.type
  qry <- paste("SELECT * FROM scores WHERE templateID = $templateID 
               AND scoreThreshold = $scoreThreshold AND", 
               type ,"is NOT NULL")
  scrs <- data.table(
    dbGetQuery(conn = conx,
               statement = qry,
               params = list(templateID = templateID,
                             scoreThreshold = score.threshold)))
  
  # Acquire and unserialize the template; access the underlying cor/binTemplate object
  template <- templatesUnserialize(db.path = db.path, templateID)@templates[[1]]  
  ncols <- template@n.t.bins + 1
  nrows <- template@n.frq.bins + 1
  
  # Unserialize features associated with these scores
  unserialized.feats <- lapply(scrs$features, 'unserialize')
  
  # If any failed to unserialize, try again: 
  check.ser <- sapply(unserialized.feats, class)
  bad.indices <- which(check.ser == 'raw')
  if (length(bad.indices)) {
    for (i in bad.indices) {
      unserialized.feats[[i]] <- unserialize(unserialized.feats[[i]])
    }
  }
  
  # Find indices of amplitude values (will be same for each event bc template is same)
  amp.inds <- grep(pattern = 'amp', names(unserialized.feats[[1]])) 
  
  # Arrange amplitude values from all detected events into an array 
  mats <- array(data = c(unlist(lapply(unserialized.feats, '[', amp.inds))), 
                dim = c(nrows, ncols, nrow(scrs)))
  
  # Take the mean of the array to get the overall average detected event:
  full.mn <- apply(mats, 1:2, mean, na.omit = TRUE)
  
  # Find Target Signals and False Alarms and find average for each type:
  ts.arr <- mats[ , , which(scrs[,get(type)] == 1)] # indices of target signals
  fa.arr <- mats[ , , which(scrs[,get(type)] == 0)] # indices of false alarms
  ts.mn <- apply(ts.arr, 1:2, mean) # mean target signal
  fa.mn <- apply(fa.arr, 1:2, mean) # mean false alarm
  
  if (class(template) == 'corTemplate') {
    
    tmp.pts <- template@pts
    
    # If using rect-select or cell-select
    if (nrow(tmp.pts) < nrows*ncols) {
      
      # Add color corrections
      dummy <- mats[,,1]
      dummy[dummy == 0] <- min(full.mn) - 10
      
      # Get original frequency range
      fr <- min(tmp.pts[,'frq']):max(tmp.pts[,'frq'])
      adj.frq.range <- 1:length(fr)
      
      # Adjust freq range indices
      orig.frq <- tmp.pts[,'frq']
      adj.frq <- adj.frq.range[as.factor(orig.frq)]
      
      # Input correct values
      for (m in 1:length(adj.frq)) {
        dummy[adj.frq[m], tmp.pts[m,'t']] <- tmp.pts[m,'amp']
      }
      
      tmp.mat <- matrix(data = dummy, nrow = nrows, ncol = ncols)
      
    }else{
      # If using a square template
      tmp.mat <- matrix(data = tmp.pts[ ,'amp'],
                        nrow = nrows, ncol = ncols)
    }
  }
  
  # For binTemplates, can only plot binary template (unless want to store or
  #  fetch the amplitude values associated with the on points of the template?)
  if (class(template) == 'binTemplate') {
    tmp.mat <- full.mn
    tmp.mat[tmp.mat != 0] <- 1
  }
  
  # Add color corrections
  full.mn[full.mn == 0] <- min(full.mn) - 10
  ts.mn[ts.mn == 0] <- min(ts.mn) - 10
  fa.mn[fa.mn == 0] <- min(fa.mn) - 10
  
  # Plot:
  
  # Set up the number of TS and FA signals
  ifelse(is.na(dim(ts.arr)[3]), ts.num <- 1, ts.num <- dim(ts.arr)[3])
  ifelse(is.na(dim(fa.arr)[3]), fa.num <- 1, fa.num <- dim(fa.arr)[3])
  
  par(mar = rep(2, 4), mfrow = c(1, 4))
  font <- 7 #TNR
  image(x = 1:ncols, y = 1:nrows, t(tmp.mat),
        yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', col = spec.col, frame = FALSE)
  title(paste0('Template Used: ', templateID), font.main = font)
  image(x = 1:ncols, y = 1:nrows, t(full.mn),
        yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', col = spec.col, frame = FALSE)
  title(paste0('Mean All (n = ', dim(mats)[3], ')'), font.main = font)
  image(x = 1:ncols, y = 1:nrows, t(ts.mn),
        yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', col = spec.col, frame = FALSE)
  title(paste0('Mean Target Signal (n = ', ts.num, ')'), font.main = font)
  image(x = 1:ncols, y = 1:nrows, t(fa.mn),
        yaxt = 'n', xaxt = 'n', xlab = '', ylab = '', col = spec.col, frame = FALSE)
  title(paste0('Mean False Alarm (n = ', fa.num, ')'), font.main = font)
  
  if (!is.na(dim(ts.arr)[3]) & !dim(ts.arr)[3])
    message('No Target Signals in this set, so none were plotted.')
  if (!is.na(dim(fa.arr)[3]) & !dim(fa.arr)[3])
    message('No False Alarms in this set, so none were plotted.')
}




#' @name pr
#' @title Calculate data points for a precision-recall plot
#' @description Internal function. In practice, please use \code{\link{plotROC}}.
#' @param predictions Probability predictions for the true positive class.
#' @param labels True class labels.
#' @return data.frame of precision and recall points to plot.
#' @seealso \code{\link{classifierModels}}
#' @references 
#' Code for this function is from: 
#' \href{https://github.com/kboyd/raucpr/blob/master/precision_recall.r}{https://github.com/kboyd/raucpr/blob/master/precision_recall.r}
#' @export


pr <- function(predictions, labels) {
  
  # code by kboyd: https://github.com/kboyd/raucpr/blob/master/precision_recall.r
  
  pos.values <- predictions[labels == 'TS']
  neg.values <- predictions[labels == 'FA']
  
  ## 0 for negative, 1 for positive
  l <- rep(c(0,1), c(length(neg.values), length(pos.values)))
  v <- c(neg.values, pos.values)
  
  ## sort by descending value
  indices <- order(v, decreasing = TRUE)
  labels <- l[indices]
  values <- v[indices]
  
  tp <- 0
  fp <- 0
  
  z <- 1
  precisions <- rep(0, length(indices) + 1)
  recalls <- rep(0, length(indices) + 1)
  
  for (i in 1:length(indices)) {
    if (labels[i] == 1) {
      tp <- tp + 1
    }
    else {
      fp <- fp + 1
    }
    
    ## make sure not between tied values
    if (i == length(indices) || values[i] > values[i + 1]) {
      ## can put a threshold
      
      ## if first, insert (r=0,p=1) point
      if (z == 1) {
        recalls[z] <- 0
        precisions[z] <- 1
        z <- z + 1
      }
      recalls[z] <- tp/(length(pos.values))
      precisions[z] <- tp/(tp + fp)
      z <- z + 1
    }
  }
  
  # Truncate recalls and precisions
  recalls <- recalls[1:(z - 1)]
  precisions <- precisions[1:(z - 1)]
  
  pr <- data.frame(cbind(precisions, recalls))
  names(pr) <- c('Precision', 'Recall')
  return(pr)
}




#' @name qryPkCheck
#' @title Check whether user-input field values are present in the database
#' @description Check whether user-input field values are present in the database (and therefore allowable as function input options). Though documented, this is generally only intended to be used as an internal function. 
#' @param conn DBI connection to an RSQLite database.
#' @param table Character string indicating the table to search within.
#' @param field Character string indicating the field name of content to search within for matches.
#' @param arg Default = field. Character string indicating the argument name (if different from field name). 
#' @param user.vals Character vector of user-given values.
#' @return If it encounters a mismatch between user-chosen values and options available in the database, \code{qryPkCheck} returns a \code{stop} error about which user inputs do not match values in the database. 
#' @family qry
#' @export
#' @examples
#'\dontrun{
#'
#' # ------------------------------------------------------------
#' # Set up a demo AMMonitor database
#' # ------------------------------------------------------------
#'
#' # Create the database (this will be deleted):
#' db.name <- 'demo.sqlite'
#'                 
#' dbCreateSample(db.name = db.name, 
#'                file.path = paste0(getwd(),'/'), 
#'                tables = c('species', 'library', 'people'))
#'                 
#' # Verify that the database exists in your current working directory:
#' file.exists(db.name)
#' 
#' # Connect to the db: 
#' db.path <- paste0(getwd(), '/', db.name)
#' conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
#'
#' # ------------------------------------------------------------
#' # Check primary key content for user.val options
#' # ------------------------------------------------------------
#'
#' # Returns a stop error due to 'stuff-not-in-the-library' 
#' # not being present in the database
#' qryPkCheck(conn = conx,
#'            table = 'library',
#'            field = 'libraryID',
#'            user.vals = c('verd_2notes', 
#'                          'stuff-not-in-the-library'))
#'                          
#' # Returns nothing (no error, no feedback) because both user.vals
#' # values are present in the database
#' qryPkCheck(conn = conx,
#'            table = 'library',
#'            field = 'libraryID',
#'            user.vals = c('verd_2notes', 
#'                          'verd_3notes'))
#'
#' # ------------------------------------------------------------
#' # Disconnect from demo database and delete from working directory
#' # ------------------------------------------------------------
#' 
#' dbDisconnect(conn = conx)
#' unlink(db.path)
#'
#'}

qryPkCheck <- function(conn, table, field, arg = field, user.vals) {
  # Check argument user-input values against options actually available in the table
  qrySelect <- paste0("SELECT ", field, " FROM ", table)
  opts <- unlist(dbGetQuery(conn, statement = qrySelect))
  bad.opts <- user.vals[which(!(user.vals %in% opts))]
  if (length(bad.opts) > 0 ) {
    stop(sprintf(paste0('The following ', field, 's in the ', arg, ' object do not exist in the ', table , ' table: \n %s'), 
                 paste(bad.opts, sep = '', collapse = ', ')))
  } 
}





#' @name qryDeployment
#' @title Query the database to return deployed devices at actively monitored locations
#' @description \code{qryDeployment} invokes SQL syntax to efficiently query an RSQLite-based AMMonitor database and return deployed equipment at actively monitored locations, as well as their associated Google accounts. Written primarily as an internal function to efficiently return equipment and account information required for automated routines, but may also be convenient for end users. See AMMonitor Documentation \href{placeholderurl}{Chapter 8: Equipment and Deployment} for details and context.
#' @param conn Class \code{SQLiteConnection} object, created by a call to \code{\link[DBI]{dbConnect}}. 
#' @param locationID Character vector of locationIDs for which to query the database for deployment information.
#' @return \code{qryDeployment} returns a data.frame of information about equipment deployed at active monitoring locations, composed of the following nine columns: equipmentID (character), locationID (character), accountID (character), dateDeployed (character), dateRetrieved (character), email (character), lat (numeric), long (numeric), tz (character).
#' @family qry
#' @export
#' @examples
#'\dontrun{
#'
#' # ------------------------------------------------------------
#' # Set up a demo AMMonitor database
#' # ------------------------------------------------------------
#' 
#' #' # Create the database (this will be deleted):
#' db.name <- 'demo.sqlite'
#'                 
#' dbCreateSample(db.name = db.name, 
#'                file.path = paste0(getwd(),'/'), 
#'                tables = c('locations', 'equipment', 'deployment',
#'                           'accounts', 'people'))
#'                 
#' # Verify that the database exists in your current working directory:
#' file.exists(db.name)
#' 
#' # Connect to the db: 
#' db.path <- paste0(getwd(), '/', db.name)
#' conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
#'
#'
#' # ------------------------------------------------------------
#' # Query database for deployment information about all actively
#' # monitored locations
#' # ------------------------------------------------------------
#'
#' all.deployment.info <- qryDeployment(conn = conx,
#'                                      locationID = 'all')
#' all.deployment.info
#'
#' # ------------------------------------------------------------
#' # Query database for deployment information about selected 
#' # locations
#' # ------------------------------------------------------------
#'
#' some.deployment.info <- qryDeployment(conn = conx,
#'                                       locationID = c('location@2', 
#'                                                     'location@10'))
#' some.deployment.info
#'
#'
#' # ------------------------------------------------------------
#' # Disconnect from demo database and delete from working directory
#' # ------------------------------------------------------------
#'
#' dbDisconnect(conx)
#' unlink(db.path)
#'
#'}

qryDeployment <- function(conn, locationID) {
  
  # Turn the SQLite foreign constraints on
  rs <- dbSendQuery(conn = conn, statement = "PRAGMA foreign_keys = ON;")
  dbClearResult(rs) 
  
  # Query to retrieve information only for actively deployed equipment 
  base.query <- "SELECT deployment.*, accounts.accountID, accounts.email, locations.lat, locations.long, locations.tz, deployment.dateRetrieved
  FROM (deployment INNER JOIN (accounts INNER JOIN equipment ON accounts.accountID = equipment.accountID) ON deployment.equipmentID = equipment.equipmentID) INNER JOIN locations ON deployment.locationID = locations.locationID"
  
  # Columns to return: 
  return.cols <- c('equipmentID', 'locationID', 'accountID', 'dateDeployed', 
                   'dateRetrieved', 'email', 'lat', 'long', 'tz')
  
  # Create and run queries depending on locationID argument
  if (length(locationID) == 1 && locationID == "all") {
    query.all <- paste0(base.query, " WHERE (((deployment.dateRetrieved) Is Null));")
    results <- dbGetQuery(conn, query.all)
  } else {
    query.loc <- paste0(base.query, " WHERE (((deployment.dateRetrieved) Is Null) 
                                      AND (deployment.locationID)= ?);")
    results <- dbGetQuery(conn, query.loc, param = list(locationID))
  }
  
  # Return desired columns as a df
  return(results[,return.cols])
  
}




#' @name scoresVerify
#' @title Manually verify detected events as target signals or false alarms
#' @description Interactively manually verify detected events as target signals or false alarms, either based on the speciesID or the libraryID. Verifications populate either the \strong{manualVerifySpeciesID} or \strong{manualVerifyLibraryID} column of the scores table, respectively, with 1s (to indicate target signals) or 0s (to indicate false alarms). See \strong{Details}. 
#' @details Verifications act as labelled training data, and pose the first step toward creating classifiers that assign a target signal probability to each detected event.
#' 
#' Verification can be performed either at the level of the species (choose \code{label.type = 'speciesID'}, which populates the manualVerifySpeciesID column of the species table) or at the level of the target signal (choose \code{label.type = 'libraryID'}, which populates the manualVerifyLibraryID column of the scores table). 
#' 
#' If we verify at the level of the speciesID, this means any detected sound produced by the target species will be labeled as a target signal. If we verify at the level of the libraryID, we are strictly seeking signals that match the libraryID associated with the template used to detect events. Sometimes it can be difficult to decide whether a detected event should count as a target signal -- generally, a research program should be careful to develop labeling standards consistent with their research objectives, and that reflect their knowledge of the target species. 
#' 
#'  When running \code{scoresVerify}, users are prompted to flag each detected event as a target signal or false alarm. To label target signals, users should input 'y' (yes) to the prompt. To label false alarms, users should input 'n' (no) to the prompt.
#'  
#' \strong{NOTE}: when verifying at the libraryID level, if a template has correctly found the call specified by the libraryID, then manualVerifySpeciesID will also automatically be verified and labeled as 1 by default.
#' @param db.path The file path that houses the SQLite database.
#' @param date.range The first option for how to specify which scores should be verified. \code{date.range} is an optional length 2 character vector of date ranges (inclusive) over which to run \code{scoresVerify}. Dates should be given in YYYY-mm-dd format. e.g. c('2016-03-04', '2016-03-12'). The alternative to using \code{date.range} is \code{recordingID}; do not use both arguments in one function call. 
#' @param recordingID The second option for how to specify which scores should be verified. \code{recordingID} is an optional character vector of recordingIDs over which to run \code{scoresVerify}. If scores should be verified for all recordings, set \code{recordingID = 'all'}. The alternative to using \code{recordingID} is \code{date.range}; do not use both arguments in one function call. 
#' @param templateID Character vector of the templateID whose scores should be verified. Only one templateID may be entered at a time. 
#' @param label.type Character string indicating whether to manually verify at the level of the speciesID or the libraryID. Options are c('speciesID', 'libraryID'). See \strong{Details}.
#' @param directory If not providing 'token.path', 'directory' is a character string indicating the local directory where recordings are stored. If providing 'token.path', 'directory' should provide the Dropbox Web folder with which to interact. In either case, 'directory' may be a single folder (e.g. 'recordings') or a nested folder (e.g. 'myfolder/recordings') so long as there are no slashes at the beginning or end of the directory string. 
#' @param token.path If using a 'directory' on Dropbox, provide a character string of the file path to an RDS file containing a Dropbox token generated with the OAuth2 connection (see \strong{Details} of \code{\link{dropboxMoveBatch}} or AMMonitor Documentation \href{placeholderurl}{Chapter 11: Recordings} for further context).
#' @param db.insert Default = FALSE. Logical flag for whether to insert verifications generated by this function into the database. If TRUE, verifications are added into the manualVerifyLibraryID and/or manualVerifyScoresID column(s) of the scores table in the the database. If FALSE, events are returned as a data.table for examination by the user, but not added to the database.
#' @param overwrite Default = FALSE. If TRUE and db.insert = TRUE, the user will overwrite any existing verifications in the manualVerifyLibraryID and/or manualVerifySpeciesID column(s) of the scores table in the database. If FALSE, no overwriting occurs, and the user only verifies events that are currently unlabelled.
#' @param fd.rat Default = 4. Ratio of plot frame (time duration of plots) to template duration. In essence, this argument specifies how much visual "context" to place around the detected event to be verified. Using smaller values (e.g., fd.rat = 1) will cause the detected event to take up nearly the entire plot (users will not see much of the surrounding spectrogram). Larger values (e.g. fd.rat = 8) will cause the detected event to be much smaller relative to the rest of the spectrogram, and users will be able to see much more of the surrounding spectrogram. 
#' @param f.lim Default = c(0, 12). Length-two numeric vector specifying frequency limits (y axis limits) for the spectrogram. 
#' @param spec.col Default = gray.3(). A vector of colors for the spectrogram. 
#' @param box Default = TRUE. Logical flag for whether to plot a box around detected events in the spectrogram. Box boundaries are based on template duration and frequency limits. Can also be set to "template" to see the template points plotted over the detection.
#' @param on.col Default = "#FFA50050". Colors for the on points of a binary point matching template, if \code{box = 'template'}. Uses #RRGGBBAA color notation, where the last two digits, AA, specify transparency.
#' @param off.col Default = "#0000FF50". Colors for the off points of a binary point matching template, if \code{box = 'template'}. Uses #RRGGBBAA color notation, where the last two digits, AA, specify transparency.
#' @param off.col Default = "#80008050". Colors for the points of a spectrogram cross correlation template, if \code{box = 'template'}. Uses #RRGGBBAA color notation, where the last two digits, AA, specify transparency.
#' @return A data.table with columns equivalent to the scores table of an AMMonitor database, with verifications added to the manualVerifyLibraryID and/or manualVerifySpeciesID column(s).
#' @family classifier
#' @seealso \href{placeholderurl}{AMMonitor Documentation Chapter 17: Classifications}, \code{\link{plotVerifications}}, \code{\link{plotVerificationsAvg}}, \code{\link{scoresDetect}}
#' @export
#' @examples
#' \dontrun{
#' 
#' # ------------------------------------------------------------
#' # Set up a demo AMMonitor database
#' # ------------------------------------------------------------
#' 
#' # Create the database (this will be deleted):
#' db.name <- 'demo.sqlite'
#' 
#' dbCreateSample(db.name = db.name, 
#'                file.path = paste0(getwd(),'/'), 
#'                tables = c('recordings', 'templates', 'locations', 
#'                           'equipment', 'people', 'accounts', 
#'                           'library', 'species', 'scores'))
#'                 
#' # Verify that the database exists in your current working directory:
#' file.exists(db.name)
#' 
#' # Connect to the db: 
#' db.path <- paste0(getwd(), '/', db.name)
#' conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
#'
#' #----------------------------------------------------------------
#' # Read in wave files and write to wd 
#' #----------------------------------------------------------------
#' # Read in sampleRecordings data
#' data(sampleRecordings)
#' 
#' # Write four waves to wd
#' writeWave(object = sampleRecordings[[1]], 
#'           filename = "midEarth3_2016-03-12_07-00-00.wav")
#' writeWave(object = sampleRecordings[[2]], 
#'           filename = "midEarth4_2016-03-04_06-00-00.wav")
#' writeWave(object = sampleRecordings[[3]], 
#'           filename = "midEarth4_2016-03-26_07-00-00.wav")
#' writeWave(object = sampleRecordings[[4]], 
#'           filename = "midEarth5_2016-03-21_07-30-00.wav")
#'
#'                        
#' #----------------------------------------------------------------
#' # Verify detected events in the scores table using the date.range
#' # argument, at the level of the libraryID (label.type = 'libraryID),
#' # without inserting records to the database
#' #----------------------------------------------------------------
#' 
#' test <- scoresVerify(db.path = db.path,
#'                      date.range = c('2016-03-04', '2016-03-21'),
#'                      templateID = 'verd1', 
#'                      label.type = 'libraryID', 
#'                      directory = getwd(), 
#'                      token.path = NULL,
#'                      db.insert = FALSE, 
#'                      fd.rat = 7)
#'                      
#' #----------------------------------------------------------------
#' # Verify detected events in the scores table using the recordingID
#' # argument, at the level of the speciesID (label.type = 'speciesID'),
#' # inserting verifications to the database by setting db.insert = TRUE
#' #----------------------------------------------------------------
#' 
#' scoresVerify(db.path = db.path,
#'              recordingID = c('midEarth3_2016-03-12_07-00-00.wav',
#'                              'midEarth4_2016-03-04_06-00-00.wav'),
#'              templateID = 'verd1', 
#'              label.type = 'speciesID', 
#'              directory = getwd(), 
#'              token.path = NULL,
#'              db.insert = TRUE)
#'              
#' # View scores table of the database to ensure verifications were added
#' # for these recordings: 
#' dbGetQuery(conx, 'SELECT * FROM scores WHERE templateID = "verd1" ')            
#'
#' # ------------------------------------------------------------
#' # Disconnect from demo database and delete from working directory
#' # ------------------------------------------------------------
#'
#' dbDisconnect(conx)
#' unlink(db.path)
#' 
#' }
#'


scoresVerify <- function(db.path,
                         date.range, 
                         recordingID,
                         templateID,
                         label.type, 
                         directory,
                         token.path = NULL, 
                         db.insert = FALSE,
                         overwrite = FALSE,
                         
                         # Args that customize the verification experience
                         fd.rat = 4,
                         f.lim = c(0, 12),
                         spec.col = gray.3(),
                         box = TRUE,
                         on.col = "#FFA50050",
                         off.col = "#0000FF50",
                         pt.col = "#80008050")
{
  # Perform argument checks
  if (missing(date.range) & missing(recordingID)) {
    stop('Please use either the date.range OR recordingID argument to specify which recordings should be verified.')
  }
  
  if (missing(label.type) | !(label.type == 'speciesID' | label.type == 'libraryID')) {
    stop("Please use the label.type argument to specify how you would like to verify these recordings. \n  Verify based on the call type ('libraryID') or the species ('speciesID').")
  }
  
  # Assign correct label.type
  ifelse(label.type == 'speciesID', type <- 'manualVerifySpeciesID', type <- 'manualVerifyLibraryID')
  
  # Connect to the database:
  conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
  
  # Turn the SQLite foreign constraints on
  rs <- dbSendQuery(conn = conx, statement = "PRAGMA foreign_keys = ON;" )
  dbClearResult(rs)
  
  # Save existing working directory to reset it after
  #  function is done running:
  owd <- setwd(getwd())
  on.exit(setwd(owd))
  
  # If a token is provided, we are using dropbox, otherwise use local directory
  ifelse(!is.null(token.path), dropbox <- TRUE, dropbox <- FALSE)
  
  # Add a forward slash to recordings directory 
  user.dir <- directory
  if (grepl("\\/$", directory) == FALSE) {
    directory <- paste0(directory, '/')
  }
  
  # Acquire and unserialize the template
  tmp <- templatesUnserialize(db.path = db.path, templateID)
  
  # Set up verification window frame
  frame <- fd.rat*tmp@templates[[templateID]]@duration
  template <- tmp@templates[[1]] #[[1]] only keep first if duplicates
  
  # IF USING DATE RANGE
  if (!missing(date.range)) {
    # If a date.range is present, select relevant recordings
    dates <- list(first = date.range[1], second = date.range[2])
    recordingID <- unlist(
      dbGetQuery(conn = conx,
                 statement = "SELECT recordingID FROM recordings 
                 WHERE startDate BETWEEN $first AND $second",
                 params = dates)
    )
  }
  
  # Acquire all scores associated with the input recordings and template
  combos <- expand.grid(templateID = templateID, 
                        recordingID = recordingID, 
                        stringsAsFactors = FALSE)
  scrs <-
    data.table(
      dbGetQuery(conn = conx,
                 statement = "SELECT * FROM scores 
                 WHERE templateID = $templateID 
                 AND recordingID = $recordingID",
                 params = as.list(combos)
      )
    )
  
  if (scrs[,.N] == 0) {
    return(message("Nothing to verify for this template at this score threshold."))
  }
  
  # If not overwriting existing verifications, only keep scrs events that still require verifying
  if (overwrite == FALSE) {
    cat('You have verified',
        length(which(!(is.na(scrs[,get(type)])))),
        'out of', scrs[,.N], 'events. \nSince overwrite == FALSE, only the remaining',
        scrs[,.N] - length(which(!(is.na(scrs[,get(type)])))),
        'events will be verified.')
    scrs <- scrs[which(is.na(scrs[,type, with = FALSE])), ]
    if (scrs[,.N] == 0) return(message('\n Overwrite == FALSE and you have already verified all records for this template at this score threshold.'))
  }
  
  # Find unique recordings in this dataset
  unique.recs <- unique(scrs$recordingID)
  counter <- 0
  
  # For each recording:
  for (r in 1:length(unique.recs)) {
    
    # Assign recording
    rec.string <- unique.recs[r]
    
    # Only look at scores for this recording 
    pks <- scrs[recordingID == rec.string]
    
    if (dropbox == TRUE) {
      
      cat('\nGetting next recording from Dropbox...')
      
      # Try catch for recordings that may not be present
      catch.error <- tryCatch(
        dropboxGetOneFile(file = rec.string, 
                          directory = directory,
                          token.path = token.path),
        error = function(e) e
      ) 
      
      if (inherits(catch.error, "error")) {
        message(paste0('\n Can\'t find ', rec.string, ' in directory "',
                       user.dir,
                       '". \nSkipping to next recording.'))
        next
      } # End tryCatch 
      
      # Set working directory as local directory
      setwd(getwd())
      file.name <- file.path(getwd(), rec.string)
      
    } else {
      file.name <- file.path(user.dir, rec.string)
    }
    
    # Read in wave file, extract helpful acoustic parameters:
    wav <- readWave(file.name)
    reclen <- length(wav@left)/wav@samp.rate
    fft.data <- monitoR:::spectro(wave = wav, wl = template@wl,
                                  wn = template@wn, ovlp = template@ovlp)
    trec <- fft.data$time
    frec <- fft.data$freq
    arec <- fft.data$amp
    which.frq.bins <- which(frec >= f.lim[1] & frec <= f.lim[2])
    frec <- frec[which.frq.bins]
    arec <- arec[which.frq.bins, ]
    verC <- NULL
    
    # For plotting, save object to help label hh:mm:ss on x-axis
    trec.times <- as.ITime(trec) 
    time.lab <- 'Time (hh:mm:ss)'
    
    # If recording is less than a minute long, keep display time in s
    if (reclen < 60) {
      trec.times <- trec
      time.lab <- 'Time (s)'
    }
    
    # Not sure what ask means here? legacy code from viewSpec interactive
    ask <- FALSE
    oldask <- par(ask = par("ask"))
    on.exit(par(oldask))
    
    # Skip this recording if no detections.
    if (pks[,.N] == 0) {
      message("No events detected within this recording. Skipping to next recording.")
      next
    }
    
    for (i in 1:pks[,.N]) {
      
      # Plot the clip to be verified:
      x <- "x"
      t.start <- max(pks$time[i] - frame/2, 0)
      t.end <- min(pks$time[i] + frame/2, reclen) #
      rec.clip <- cutWave(wave = wav, from = t.start, to = t.end)
      time.inds <- which(trec %in% trec[trec >= t.start & trec <= t.end])
      times <- trec[time.inds]
      amp.clip <- arec[, trec %in% times]
      par(mfrow = c(1,1), mar = c(3,3,2,1), mgp = c(2,1,0))
      image(x = times, y = frec, z = t(amp.clip), col = spec.col,
            xlab = time.lab, ylab = "Frequency (kHz)", xaxt = "n",
            bty = 'n', axes = FALSE,
            main = paste("Peak",i, 'of',nrow(pks),
                         'in this recording. Score = ',
                         round(pks$score[i], 2)))
      if (reclen > 60)
        axis(1, at = pretty(times), labels = as.ITime(pretty(times)))
      else 
        axis(1, at = pretty(times), labels = pretty(times))
      
      # Draw a box around the verification area:
      if (box == TRUE & nrow(pks) > 0) {
        xleft <- pks$time[i] - template@duration/2
        xright <- pks$time[i] + template@duration/2
        ylwr <- template@frq.lim[1]
        yupr <- template@frq.lim[2]
        polygon(x = c(xleft, xleft, xright, xright),
                y = c(ylwr, yupr, yupr, ylwr), border = "blue")
      }
      
      # If no box
      else if (tolower(box) == "template" & nrow(pks) > 0) {
        xleft <- pks$time[i] - template@duration/2
        ylwr <- template@frq.lim[1]
        if (class(template) == "binTemplate") {
          pt.on <- template@pt.on
          pt.off <- template@pt.off
          pt.on[, "t"] <- pt.on[, "t"] + ((xleft - min(times))/template@t.step)
          pt.off[, "t"] <- pt.off[, "t"] + ((xleft - min(times))/template@t.step)
          pt.on[, "frq"] <- pt.on[, "frq"] + ylwr - 1
          pt.off[, "frq"] <- pt.off[, "frq"] + ylwr - 1
          bin.amp <- 0 * amp.clip
          bin.amp[pt.on[, c(2, 1)]] <- 1
          bin.amp[pt.off[, c(2, 1)]] <- 2
          image(x = times, y = frec, t(bin.amp),
                zlim = c(0,2),
                col = c("transparent", on.col, off.col),
                add = TRUE)
        }
        else if (class(template) == "corTemplate") {
          pts <- template@pts
          pts[, "t"] <- pts[, "t"] + ((xleft - min(times))/template@t.step)
          pts[, "frq"] <- pts[, "frq"] + ylwr - 1
          bin.amp <- 0 * amp.clip
          bin.amp[pts[, c(2, 1)]] <- 1
          image(x = times, y = frec, t(bin.amp),
                zlim = c(0, 1), col = c("transparent", pt.col),
                add = TRUE)
        }
        else stop("Template list class not recognized: ",
                  class(template))
      }
      
      while (length(x) == 0 || !x %in% c("y", "n", NA)) {
        
        cat(paste0("\n This is recording ", r, " out of ", 
                   length(unique.recs), ": ", rec.string, "\n",
                   " The scoreID for this verification is: ", pks[i,scoreID], "\n",
                   " This is verification ", counter + i, 
                   " out of ", scrs[,.N]), "\n")
        cat(paste0("\n", i, ". True detection for this ", 
                   label.type,
                   "?\n Enter y for yes, n for no, NA for NA, or q to exit (q will exit and save any verifications you have already completed): "))
        x <- tolower(readLines(n = 1)[1])
        
        if (length(x) == 0) {
          cat("\nYou didn't enter a response.\n")
          next
        }
        if (!is.na(x) && x == "na")
          x <- NA
        if (is.na(x)) {
          cat("NA\n")
          break
        }
        cat(switch(x, n = FALSE, 
                   y = TRUE, 
                   q = "Exiting and saving what you have already verified.", "Value not recognized. Enter y, n, NA, or q."),
            "\n")
        if (!x %in% c("y", "n", "q"))
          next
        if (x == "q")
          break
      } # end while x %in% y, n, NA
      
      # If q, break out of the peaks loop
      if (!is.na(x) & x == 'q') break
      
      if (is.na(x) || x != "r")
        verC[i] <- x
      par(ask = ask)
      if (!is.na(x) && x == "r")
        i <- i - 1
      else i <- i + 1
      if (i < 1)
        i <- 1
    } # end for i in 1:pks[,.N]
    
    counter <- counter + pks[,.N]
    
    cat("\n")
    
    # Update verification labels in scores table:
    verC[verC == 'y'] <- 1
    verC[verC == 'n'] <- 0
    vers <- as.integer(verC)
    update.these <- pks[,scoreID]
    # If we quit before the end of a recording, we only 
    # update labels for the ones we actually verified
    if (x == 'q') {
      update.these <- pks[,scoreID][seq_along(vers)]
    }
    
    # If there is anything to update
    if (length(update.these) > 0) {
      scrs[scoreID %in% update.these, (type) := vers]
      
      # Update verifications directly in database if desired
      if (db.insert) {
        params <- list(scoreID = update.these, vers = vers)
        if (type == 'manualVerifySpeciesID') {
          dbExecute(conn = conx, statement = 
                      "UPDATE scores 
                    SET manualVerifySpeciesID = $vers
                    WHERE scoreID = $scoreID",
                    param = params)
        }
        
        if (type == 'manualVerifyLibraryID') {
          dbExecute(conn = conx, statement = 
                      "UPDATE scores 
                    SET manualVerifyLibraryID = $vers
                    WHERE scoreID = $scoreID",
                    param = params)
          
        }
      }
    } # end if length(update.these)
    
    cat("Finished verifying", label.type, "for this recording.\n")
    
    # Remove temporarily downloaded file if using web dropbox
    if (dropbox == TRUE) {
      cat(paste0('\nRemoving Dropbox Web wave file from local directory.'))
      file.remove(file.name)
    }
    
  } # end recordings loop
  
  # If libraryID is labeled TRUE, so is manualVerifySpeciesID
  if (db.insert == TRUE) {
    if (length(update.these) > 0) {
      dbExecute(conn = conx, statement = 
                  "UPDATE scores 
                SET manualVerifySpeciesID = 1
                WHERE manualVerifyLibraryID = 1
                AND scoreID = $scoreID",
                param = list(scoreID = update.these))
      cat('\nVerifications saved to database.')
    }
  }
  scrs[manualVerifyLibraryID == 1, manualVerifySpeciesID := 1]
  
  # Return scores table updated with verifications
  return(scrs)
} 



#' @name templatesUnserialize
#' @title Extract templates from an AMMonitor SQLite database
#' @description Extract and unserialize templates from a SQLite database for subsequent use in R. SQLite databases cannot store specialized R S4 objects like templates without first converting them into BLOBs (binary large objects used to store any kind of data of unlimited size). \code{templatesUnserialize} converts template BLOBs back into templateLists for use in R.
#' @param db.path The file path that houses the SQLite database.
#' @param templateID Character vector of templateIDs to convert back into binTemplateList or corTemplateList objects. Must match templateIDs present in the database. Additionally, bin and cor template types may not be mixed in one function call.
#' @return S4 object of class \code{binTemplateList} or \code{corTemplateList}
#' @family templates
#' @seealso \code{\link{templatesInsert}}, \href{placeholderurl}{AMMonitor Documentation Chapter 15: Templates}
#' @export
#' @examples
#' \dontrun{ 
#'

#' # ------------------------------------------------------------
#' # Set up a demo AMMonitor database
#' # ------------------------------------------------------------
#' 
#' # Create the database (this will be deleted):
#' db.name <- 'demo.sqlite'
#'                 
#' dbCreateSample(db.name = db.name, 
#'                file.path = paste0(getwd(),'/'), 
#'                tables = c('library', 'species', 'people'))
#'                 
#' # Verify that the database exists in your current working directory:
#' file.exists(db.name)
#' 
#' # Connect to the db: 
#' db.path <- paste0(getwd(), '/', db.name)
#' 
#' conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
#'
#' #----------------------------------------------------------------
#' # Read in a wave file 
#' #----------------------------------------------------------------
#' # Read in sampleRecordings data
#' data(sampleRecordings)
#'  
#' # Write the fourth recording to the working directory
#' writeWave(object = sampleRecordings[[4]], 
#'           filename = "midEarth5_2016-03-21_07-30-00.wav")
#' 
#' 
#' #----------------------------------------------------------------
#' # Create some cor and bin templates
#' #----------------------------------------------------------------
#' verd1 <- makeCorTemplate(clip = "midEarth5_2016-03-21_07-30-00.wav", 
#'                          t.lim = c(4.45, 4.95), 
#'                          frq.lim = c(3.8,6), 
#'                          select = "auto", 
#'                          name = 'verd1')
#'                          
#' verd2 <- makeCorTemplate(clip = "midEarth5_2016-03-21_07-30-00.wav", 
#'                          t.lim = c(8.8, 9.32), 
#' 
#'                          lim = c(3.8,6),   
#'                          select = "auto", 
#'                          name = 'verd2')
#'                          
#' verd3 <- makeBinTemplate(clip = "midEarth5_2016-03-21_07-30-00.wav", 
#'                          t.lim = c(8.8, 9.32), 
#'                          frq.lim = c(3.8,6),   
#'                          amp.cutoff = -28,
#'                          select = "auto", 
#'                          name = 'verd3')
#'                          
#' #----------------------------------------------------------------
#' # Insert a list of corTemplates
#' #----------------------------------------------------------------
#' templatesInsert(db.path = db.path, 
#'                 template.list = combineCorTemplates(verd1, verd2), 
#'                 libraryID = c('verd_2notes', 'verd_other'), 
#'                 personID = c('fbaggins', 'fbaggins'))
#' 
#' #----------------------------------------------------------------
#' # Insert a single binTemplate
#' #----------------------------------------------------------------
#' templatesInsert(db.path = db.path, 
#'                 template.list = verd3, 
#'                 libraryID = 'verd_2notes',
#'                 personID = 'bbaggins')
#' 
#' #----------------------------------------------------------------
#' # Return a corTemplateList from the templates table blobs
#' #----------------------------------------------------------------
#' ctmps <- templatesUnserialize(db.path = db.path,
#'                               templateID = c('verd1', 'verd2'))
#' 
#' #----------------------------------------------------------------
#' # Return a binTemplateList from the templates table blobs
#' #----------------------------------------------------------------
#' btmps <- templatesUnserialize(db.path = db.path,
#'                               templateID = 'verd3')
#' 
#' # ------------------------------------------------------------
#' # Disconnect from demo database and delete from working directory
#' # ------------------------------------------------------------
#'
#' dbDisconnect(conx)
#' unlink(db.path)
#' 
#'}

templatesUnserialize <- function(db.path, 
                                 templateID){
  
  # Connect to the database:
  conx <- dbConnect(drv = dbDriver('SQLite'), dbname = db.path)
  
  # Turn the SQLite foreign constraints on
  rs <- dbSendQuery(conn = conx, statement = "PRAGMA foreign_keys = ON;" )
  dbClearResult(rs)
  
  # Check user-supplied templateIDs against those in the templates table (stop if wrong)
  qryPkCheck(conn = conx, table = 'templates', field = 'templateID', user.vals = templateID) 
  
  # Gather desired records from templates table
  get.templates <- dbGetQuery(
    conn = conx, 
    statement = 'SELECT templateID, class, template FROM  templates where templateID = ?', 
    param = list(templateID)
  )
  
  # Ensure all templates are of the same class
  cl <- unique(get.templates$class)
  if (length(cl) > 1) {
    message('Class mismatch in the templateIDs you have chosen:')
    print(get.templates)
    stop('Please ensure all templateIDs in the templateID argument point to templates of the same class \n(i.e., either all binTemplateList or all corTemplateList.')
  }
  
  # Unserialize the templates
  unser <- lapply(get.templates$template, 'unserialize')
  
  # Convert unserialized templates back into bin or corTemplateList object
  names(unser) <- templateID # must have a name slot
  template.list <- new(Class = cl, templates = unser)
  
  # Return a bin or corTemplateList of the unserialized templates
  return(template.list)
}





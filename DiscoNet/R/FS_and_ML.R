#' Internal function: Check for validitity of ML functions
#'
#' Return an error and stop the function if arguments are not how they should be
#' @param df A numerical dataframe
#' @param class A boolean vector representing the class (TRUE/FALSE)
#' @param nbr_of_boot A natural number higher than 1, the number of bootstraps done.
#' @param seed seed for reproductability
#' @param nfolds number of folds during cross validation
#' @param nCores number of cores for parallel computing
#' @param verbose logical, wether to print or not internal calculations.
#' @return Return an error, or let the function continue if everything is fine.
check_for_validity2 = function(df, class, nbr_of_boot = NULL ,seed = NULL, nfolds = NULL, nCores = NULL, verbose = NULL){
  if(!is.data.frame(df)){
    stop("df is not a dataframe!")
  }
  if(any(is.na(df))){
    stop("df contains missing values")
  }
  if(!all(unlist(lapply(df, is.finite)))){
    stop("df contains not finite values")
  }
  if(any(c("character","logical","factor") %in% sapply(df, class))){
    stop("df cannot contains character, logical, or factor variables. Please convert the class type of  your features!")
  }
  if(!is.logical(class)){
    stop("class should be a logical vector of TRUE and FALSE!")
  }
  if(nrow(df) != length(class)){
    stop("nrow(df) is not equal to length(class)!")
  }
  if(!is.null(nbr_of_boot)){
    if(!is.numeric(nbr_of_boot)){
      stop("nbr_of_boot is not recognised")
    }
    if((nbr_of_boot < 1)){
      stop("nbr_of_boot should be higher than 1. Try 10!")
    }
  }
  if(!is.null(seed)){
    if(!is.numeric(seed)){
      stop("seed is not recognised")
    }
  }
  if(!is.null(nfolds)){
    if(!is.numeric(nfolds)){
      stop("nbr_of_boot is not recognised")
    }
    if((nfolds < 1)){
      stop("nbr_of_boot should be higher than 1. Try 10!")
    }
  }
  if(!is.null(verbose)){
    if(!is.logical(verbose)){
      stop("verbose is not recognised")
    }
  }
  if(!is.null(nCores)){
    if(!is.numeric(nCores)){
      stop("nCores is not recognised")
    }
    if((nCores < 1)){
      stop("nCores should be higher than 1. Try 1, or more for parallel computing!")
    }
  }
}



#' Calculate information gain
#'
#' Calculate the information gain of each features (column) in a dataframe
#' @param df A numerical dataframe or matrix
#' @param class A boolean vector representing the class (TRUE/FALSE)
#' @param verbose Default = TRUE, print internal calculations.
#' @return A dataframe of two columns: InfGain (information gain value) and feature (name of the feature).
#' @examples
#' if(requireNamespace("mlbench", quietly = TRUE)) {
#' data("BreastCancer", package = "mlbench")
#' }
#'
#' # Remove rows with missing values
#' BreastCancer = BreastCancer[rowSums(is.na(BreastCancer)) == 0, ]
#'
#' # Remove the ID column and the class column.
#' df = BreastCancer[, -c(1, 11)]
#'
#' # Convert each column to numeric class
#' df = as.data.frame(apply(df, 2, as.numeric))
#'
#' # Retrieve the class as a logical vector
#' class = BreastCancer$Class == "malignant"
#'
#' # Calculate information gain for each features
#' BreastCancer_InfoGain = InformationGain(df, class)
#' @export
#' @import mlbench
#' @importFrom magrittr %>%
InformationGain = function(df, class, verbose = TRUE){
  if(verbose){message("InformationGain is running...")}
  check_for_validity2(df, class)
  feature_names = colnames(df)
  df = data.frame(class = as.factor(class), df)
  if(ncol(df) < 20){
    Results = utils::stack(RWeka::InfoGainAttributeEval(class ~ ., df, na.action = "na.omit"))
  } else {
    chunksy = chunks(2:ncol(df), round((ncol(df)-1)/10))
    Results = list()
    ni=0
    for(chunk in chunksy){
      ni = ni+nbr(chunk)
      add(Results, utils::stack(RWeka::InfoGainAttributeEval(class ~ ., df[, c(1,chunk)], na.action = "na.omit")))
      if(verbose){cat('\r', paste0(ni,'/',ncol(df), " features"))}
    }
    Results = do.call(rbind, Results)
  }
  best_InfGain = InformationGain_perfect(class = as.factor(class))
  Results = Results %>%
    `colnames<-`(c("InfGain", "feature")) %>%
    dplyr::mutate(feature = feature_names) %>%
    dplyr::arrange(dplyr::desc(InfGain)) %>%
    dplyr::mutate(InfGain = (InfGain - min(InfGain)) / (best_InfGain - min(InfGain))) %>%
    dplyr::mutate(InfGain = signif(InfGain, 5))
  if(verbose){message("\nDone.")}
  return(Results)
}


if(requireNamespace("mlbench", quietly = TRUE)) {
  data("BreastCancer", package = "mlbench")
}

#' Calculate information gain through bootstraping
#'
#' Calculate the information gain of each features (column) in a dataframe using multiple bootstraping
#' @param df A numerical dataframe or matrix
#' @param class A boolean vector representing the class (TRUE/FALSE)
#' @param nbr_of_boot Number of bootstraps
#' @param seed Optional. For reproducability.
#' @param verbose Default = TRUE, print internal calculations.
#' @return A dataframe of with for each feature, the information gain calculated for each bootstraps
#' @examples
#' if(requireNamespace("mlbench", quietly = TRUE)) {
#' data("BreastCancer", package = "mlbench")
#' }
#'
#' # Remove rows with missing values
#' BreastCancer = BreastCancer[rowSums(is.na(BreastCancer)) == 0, ]
#'
#' # Remove the ID column and the class column.
#' df = BreastCancer[, -c(1, 11)]
#'
#' # Convert each column to numeric class
#' df = as.data.frame(apply(df, 2, as.numeric))
#'
#' # Retrieve the class as a logical vector
#' class = BreastCancer$Class == "malignant"
#'
#' # Calculate information gain for each features
#' BreastCancer_InfoGain = InformationGain_Bootstrap(df, class, nbr_of_boot = 5)
#' @export
#' @import mlbench
#' @importFrom magrittr %>%
InformationGain_Bootstrap = function(df, class, nbr_of_boot = 2, seed = NULL, verbose = TRUE){
  if(verbose){message("InformationGain_Bootstrap is running...")}
  check_for_validity2(df, class, nbr_of_boot = nbr_of_boot, seed = seed)
  get_bootstrap_index = function(class, nbr_of_boot){
    list_index = list()
    for(i in 1:nbr_of_boot){
      rows = 1:length(class)
      classes = unique(class)
      row_class1 = rows[class == classes[1]]
      row_class2 = rows[class == classes[2]]

      new_row_class1 = sample(x = row_class1, size = length(row_class1), replace = TRUE)
      new_row_class2 = sample(x = row_class2, size = length(row_class2), replace = TRUE)

      new_rows = c(new_row_class1, new_row_class2)
      add(list_index, new_rows)
    }
    return(list_index)
  }
  set.seed(seed)
  Boot_index = get_bootstrap_index(class, nbr_of_boot)
  results_infgain = list()
  i = 0
  for(index in Boot_index){
    i = i+1
    add(results_infgain, InformationGain(df = df[index, ], class = class[index], verbose = FALSE))
    if(verbose){cat('\r', paste0(i,'/',nbr(Boot_index), " bootstraps"))}
  }
  names(results_infgain) = paste0("Boot_", 1:length(results_infgain))
  Summary_infgain_rep = dplyr::bind_rows(results_infgain, .id = "DataFrame")
  Summary_infgain_rep = tidyr::pivot_wider(Summary_infgain_rep, names_from = DataFrame, values_from = InfGain)
  Summary_infgain_rep = Summary_infgain_rep %>%
    dplyr::mutate(RowMean = rowMeans(dplyr::select(., tidyselect::starts_with("Boot_")), na.rm = TRUE)) %>%
    dplyr::mutate(RowSd = {apply(dplyr::select(., tidyselect::starts_with("Boot_")), 1, function(x) stats::sd(x,na.rm = TRUE))}) %>%
    dplyr::select(feature, RowMean,RowSd, tidyselect::everything()) %>%
    dplyr::arrange(dplyr::desc(RowMean))
  if(verbose){message("\nDone.")}
  return(Summary_infgain_rep)
}



#' Internal function: calculate the maximum information gain for a feature
#'
#' #' Calculate the maximum information gain for a dataset to later normalize the values obtained so that they go to 1.
#' @param class A boolean vector representing the class (TRUE/FALSE)
#' @return A single value representing the maximum information gain.
InformationGain_perfect = function(class){
  df = data.frame(class = as.factor(class), perfect_feature = as.factor(class))
  return(RWeka::InfoGainAttributeEval(class ~ ., df, na.action = "na.omit"))
}



#' Internal function: Do Adaptive LASSO for feature selection
#'
#' Return the list of features chosen through Adaptive LASSO during cross validation
#' @param index A numerical vector indicating which rows of the dataframe to work on
#' @param Dataset A numerical dataframe with features as columns and observations as rows
#' @param Label A boolean vector representing the class (TRUE/FALSE)
#' @param nfolds The number of folds during cross validation for tuning lambda
#' @param i Optional. Used during bootstraping
#' @param verbose Default to TRUE.Printing internal calculations.
#' @return A string vector returning the selected features.
run_lasso = function(index, Dataset, Label, nfolds, i = NULL, verbose = TRUE){
  # Launch models
  label = Label[index]
  weights = ifelse(label == 1, sum(label == 0)/length(label), sum(label == 1)/length(label))
  if(is.null(i) & verbose){
    message("Calculating adaptive Ridge weight")
  } else {
    if(verbose){message(paste0("Calculating adaptive Ridge weights, CV: ", i))
    }
  }
  cv_ridge = glmnet::cv.glmnet(x = as.matrix(Dataset[index, ]), y = label, weights = weights, alpha = 0, nfolds = nfolds, family = "binomial", type.measure = "auc")
  ada_weights = 1/abs(matrix(coef(cv_ridge, s=cv_ridge$lambda.min)[, 1][-1]))
  if(is.null(i) & verbose){
    message("Running LASSO")
  } else {
    if(verbose){message(paste0("Running LASSO, CV : ", i))
    }
  }
  cv.lasso = glmnet::cv.glmnet(x = as.matrix(Dataset[index, ]), y = label, weights = weights, alpha = 1, nfolds = nfolds, family = "binomial", type.measure = "auc", penalty.factor = ada_weights)
  coef <- coef(cv.lasso, s='lambda.1se')
  features = rownames(coef)[-1][coef[, 1][-1] != 0]
  return(features)
}



#' Adaptive LASSO for feature selection
#'
#' Return the list of features chosen through Adaptive LASSO during cross validation
#' @param df A numerical dataframe with features as columns and observations as rows
#' @param class A boolean vector representing the class (TRUE/FALSE)
#' @param nfolds The number of folds during cross validation for tuning lambda
#' @param verbose Default to TRUE.Printing internal calculations.
#' @return A string vector returning the selected features.
#' @examples
#' if(requireNamespace("mlbench", quietly = TRUE)) {
#' data("BreastCancer", package = "mlbench")
#' }
#'
#' # Remove rows with missing values
#' BreastCancer = BreastCancer[rowSums(is.na(BreastCancer)) == 0, ]
#'
#' # Remove the ID column and the class column.
#' df = BreastCancer[, -c(1, 11)]
#'
#' # Convert each column to numeric class
#' df = as.data.frame(apply(df, 2, as.numeric))
#'
#' # Retrieve the class as a logical vector
#' class = BreastCancer$Class == "malignant"
#'
#' # Calculate selected features
#' BreastCancer_AdaLASSO = AdaLASSO(df, class, nfolds = 10)
#' @export
#' @import mlbench
#' @importFrom magrittr %>%
AdaLASSO = function(df, class, nfolds = 10, verbose = TRUE){
  if(verbose){message("select_features_ADA is running...")}
  check_for_validity2(df, class, nfolds = nfolds, verbose = verbose)
  class = as.numeric(class)
  Results = run_lasso(index = 1:length(class), Dataset = df, Label = class, nfolds = nfolds, verbose=verbose)
  if(verbose){message("Done.")}
  return(Results)
}



#' Adaptive LASSO for feature selection with bootstraping
#'
#' Return the list of features chosen through Adaptive LASSO during cross validation
#' @param df A numerical dataframe with features as columns and observations as rows
#' @param class A boolean vector representing the class (TRUE/FALSE)
#' @param nfolds The number of folds during cross validation for tuning lambda
#' @param nbr_of_boot The number of bootstraps to perform
#' @param seed Optional. For reproducibility
#' @param nCores Default to 1. For parallel computing on each bootstrap.
#' @param verbose Default to TRUE.Printing internal calculations.
#' @return A string vector returning the selected features.
#' @examples
#' if(requireNamespace("mlbench", quietly = TRUE)) {
#' data("BreastCancer", package = "mlbench")
#' }
#'
#' # Remove rows with missing values
#' BreastCancer = BreastCancer[rowSums(is.na(BreastCancer)) == 0, ]
#'
#' # Remove the ID column and the class column.
#' df = BreastCancer[, -c(1, 11)]
#'
#' # Convert each column to numeric class
#' df = as.data.frame(apply(df, 2, as.numeric))
#'
#' # Retrieve the class as a logical vector
#' class = BreastCancer$Class == "malignant"
#'
#' # Calculate selected features
#' BreastCancer_AdaLASSO = AdaLASSO_Bootstrap(df, class, nbr_of_boot = 7)
#' @export
#' @import mlbench
#' @importFrom magrittr %>%
AdaLASSO_Bootstrap = function(df, class, nfolds = 10, nbr_of_boot = 10, seed = NULL, nCores = 1, verbose = TRUE){
  if(verbose){message("select_features_ADA is running...")}
  check_for_validity2(df, class, nfolds = nfolds, nbr_of_boot = nbr_of_boot, seed = seed, nCores = nCores, verbose = verbose)
  get_bootstrap_index = function(class, nbr_of_boot){
    list_index = list()
    for(i in 1:nbr_of_boot){
      rows = 1:length(class)
      classes = unique(class)
      row_class1 = rows[class == classes[1]]
      row_class2 = rows[class == classes[2]]

      new_row_class1 = sample(x = row_class1, size = length(row_class1), replace = TRUE)
      new_row_class2 = sample(x = row_class2, size = length(row_class2), replace = TRUE)

      new_rows = c(new_row_class1, new_row_class2)
      add(list_index, new_rows)
    }
    return(list_index)
  }
  run_lasso = function(index, Dataset, Label, nfolds, i){
    # Launch models
    label = Label[index]
    weights = ifelse(label == 1, sum(label == 0)/length(label), sum(label == 1)/length(label))
    if(verbose){message(paste0("Calculating adaptive Ridge weights, Boot: ", i))}
    cv_ridge = glmnet::cv.glmnet(x = as.matrix(Dataset[index, ]), y = label, weights = weights, alpha = 0, nfolds = nfolds, family = "binomial", type.measure = "auc")
    ada_weights = 1/abs(matrix(coef(cv_ridge, s=cv_ridge$lambda.min)[, 1][-1]))
    if(verbose){message(paste0("Running LASSO, Boot : ", i))}
    cv.lasso = glmnet::cv.glmnet(x = as.matrix(Dataset[index, ]), y = label, weights = weights, alpha = 1, nfolds = nfolds, family = "binomial", type.measure = "auc", penalty.factor = ada_weights)

    coef <- coef(cv.lasso, s='lambda.1se')
    features = rownames(coef)[-1][coef[, 1][-1] != 0]
    return(features)
  }
  class = as.logical(class)
  if(!all(unique(class) %in% c("TRUE", "FALSE"))){
    stop("class not recognised, should be a vector of logicals.")
  }
  class = as.numeric(class)
  set.seed(seed)
  Boot_index = get_bootstrap_index(class, nbr_of_boot)
  Results = parallel::mclapply(seq_along(Boot_index), function(i) {print(i) ; run_lasso(index = Boot_index[[i]], Dataset = df, Label = class, nfolds = nfolds, i = i)}, mc.cores = nCores)
  names(Results) = paste0("Boot_", 1:length(Boot_index))
  Results = data.frame(Features = colnames(df), lapply(Results, function(x) as.numeric(colnames(df) %in% x)))
  if(verbose){message("Done.")}
  return(Results)
}




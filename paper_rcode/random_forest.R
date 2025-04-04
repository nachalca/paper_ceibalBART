library(tidymodels)

#Train a RF model 
random_forest <- function(train, test, default = F){

  cei_train_rf <- train |> mutate(aprobado = factor(aprobado, levels=c(1,0)))
  cei_test_rf <- test |> mutate(aprobado = factor(aprobado, levels=c(1,0)))
  
  i <- max(unlist(Filter(Negate(is.na), lapply(colnames(cei_train_rf), function(x){
    suppressWarnings(as.integer(strsplit(x,'_')[[1]][2]))
  }))))
  
  rf_best_fit <- c()
  
  model_name <- ifelse(default, paste0('models/rf/rf_default_',i), paste0('models/rf/rf_',i))
  
  if(!file.exists(model_name)){
    rf_recipe <- recipe(aprobado ~ ., data = cei_train_rf) |> 
      update_role(id_persona, new_role = "ID") |>
      update_role(c(puntos,nivel), new_role = "dependiente")
    
    rf_model <-
      rand_forest(
        mtry = tune(), 
        trees = tune(),
        min_n = tune()
      ) %>%
      set_mode("classification") %>%
      set_engine("ranger")
    
    rf_workflow <- workflow() |>
      add_recipe(rf_recipe) |>
      add_model(rf_model)

    rf_best <- NULL
    
    if(!default){
      message(paste0('Training RF for month ', i))
      
      rf_folds <- vfold_cv(cei_train_rf, v = 5)
      
      rf_param <- rf_workflow |> 
        extract_parameter_set_dials() |>
        update( mtry = mtry(c(1, ncol(cei_train_rf) - 4)), #Subtract 4: id_persona, aprobado, nivel, puntos.
                min_n = min_n(),
                trees = trees())
      
      #Uso hipercubos latinos para armar la grilla.
      rf_grid <- rf_param |> grid_latin_hypercube(size = 8)
      
      f_meas_beta <- function(data, truth, estimate, na_rm = TRUE, ...) {
        f_meas(
          data = data,
          truth = !!rlang::enquo(truth),
          estimate = !!rlang::enquo(estimate),
          # set bias = TRUE
          beta = 1.6,
          na_rm = na_rm,
          ...
        ) |> mutate(.metric = 'f_meas_beta')
      }
      
      # Use `new_numeric_metric()` to formalize this new metric function
      f_meas_beta <- new_class_metric(f_meas_beta, "maximize")
      

      #Ajusta todos los modelos de la grilla usando todos los folds.
      initial_grid <- tune_grid(
        rf_workflow,
        resamples = rf_folds,
        grid = rf_grid,
        metrics = metric_set(accuracy),  # Pass the custom function
        control = control_resamples(save_pred = TRUE, save_workflow = TRUE)
      )
      
      ctrl <- control_bayes(verbose = TRUE, uncertain = 10, no_improve = 50)
      
      tune_res <- tune_bayes(rf_workflow,
                             resamples = rf_folds,
                             initial = initial_grid,
                             param_info = rf_param,
                             iter = 200,
                             metrics = metric_set(accuracy),
                             control = ctrl)
      
      rf_best <- select_best(tune_res, metric = "accuracy")
    }    
    else {
      message(paste0('Training RF with the default hyperparamaters for month ', i))
      # https://arxiv.org/pdf/1804.03515
      rf_best <- tibble(trees = 500, mtry = round(sqrt(ncol(cei_train_rf) -4)), min_n = 1L) 
    }

    rf_workflow <- rf_workflow |> finalize_workflow(rf_best)
    rf_best_fit <- fit(rf_workflow, cei_train_rf)
    
    saveRDS(rf_best_fit, file = model_name)
    
  }
  else{
    message(paste0('loading model for month ', i, ' since the model was already trained'))
    rf_best_fit <- readRDS(file = model_name)
  }
  
  rf_best_predicted <- predict(rf_best_fit, cei_test_rf, type = 'prob') |> select(.pred_1)
}

random_forest_roc <- function(test, month){
  test <- test |> mutate(aprobado = factor(aprobado, levels=c(1,0)))
  
  rf_fit <- readRDS(file = paste0('models/rf/rf_',month))
  
  predicted_class_prob <- predict(rf_fit, test, type = "prob")
  
  res_class_prob <- test |>
    select(aprobado, puntos) |>
    bind_cols(predicted_class_prob)

  roc_curve(res_class_prob, aprobado, .pred_1)
}


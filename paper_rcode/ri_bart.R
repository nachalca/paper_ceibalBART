library(tidymodels)
library(data.table)
library(dbarts)
library(ParBayesianOptimization)


ri_bart <- function(train, test, default = F){
  
  cei_train_bart <- train |> select(-c(puntos,nivel, id_persona)) |> mutate(id_centro = as.factor(id_centro))
  cei_test_bart <- test |> select(-c(puntos,nivel, id_persona)) |> mutate(id_centro = as.factor(id_centro))
  
  i <- max(unlist(Filter(Negate(is.na), lapply(colnames(cei_train_bart), function(x){
    suppressWarnings(as.integer(strsplit(x,'_')[[1]][2]))
  }))))
  
  best_ri_bart <- NA
  model_name <- ifelse(default, paste0("models/riBART/riBART_default_",i), paste0("models/riBART/riBART_",i))
  
  if(!file.exists(model_name)){
    
      if(!default){
        message(paste0('Training riBART for month ', i))
        
        #Define folds for CV
        n_folds <- 5
        folds <- vfold_cv(cei_train_bart, n_folds, strata = aprobado)
      
        # score <- function(true, pred){
        #   spec <- sum(if_else(true == 0 & pred == 0, 1, 0))/sum(if_else(true == 0, 1, 0))
        #   sens <- sum(if_else(true == 1 & pred == 1, 1, 0))/sum(if_else(true == 1, 1, 0))
        #   (1+1.6^2)*(spec*sens)/(1.6*spec + sens)
        # }
        
        score <- function(true, pred){
          sum(if_else(true == pred, 1, 0))/length(true)
        }
        
        #Search space
        bounds <- list( 
          power = c(1, 3),
          base = c(0.0000001, 0.9999999), 
          k = c(0.1, 5),
          n.trees = c(50L, 500L)
        )
        
        scoringFunction <- function(power, base, k, n.trees) {
          scores <- sapply((1:n_folds), function(x) { #Do the CV
            
            #Extract train / validation set from the folds
            train_bart <- training(folds[[1]][[x]])
            validation_bart <- testing(folds[[1]][[x]])

            rbartFit <- rbart_vi(aprobado ~ . - id_centro, data = train_bart, 
                                 test = validation_bart, group.by = id_centro, group.by.test = id_centro,
                                 power = power, base = base,
                                 k = k,
                                 n.samples = 200, n.burn = 1000,
                                 keepTrees = TRUE, keepCall = TRUE, keepTestFits = TRUE,
                                 n.trees = 200,  n.chains = 4, n.thin = 1, seed = 1784)    
            
            pred <- if_else(rbartFit$yhat.test.mean < 0, 0, 1)
            score(validation_bart$aprobado, pred)
          })
          
          return(list(Score = mean(scores)))
        }  
        
        optObj <- bayesOpt(
          FUN = scoringFunction, 
          bounds = bounds, 
          initPoints = 5, 
          iters.n = 100, 
          iters.k = 1
        )
        
        best_params <- getBestPars(optObj)
      }
      else {
        message(paste0('Training riBART model with the default hyperparamaters for month ', i))
        best_params <- list(power = 2.0, base = 0.95, k = 2, n.trees = 200)
      }
      
      best_ri_bart <- rbart_vi(aprobado ~ . - id_centro, 
                               data = cei_train_bart,test = cei_test_bart, 
                               group.by = id_centro, group.by.test = id_centro,
                               power = best_params$power, base = best_params$base, #Hyperparameters
                               k = best_params$k,  n.trees = best_params$n.tree, #Hyperparameters
                               keepTrees = TRUE, keepCall = TRUE, keepTestFits = TRUE,
                               seed = 1784)
      
      saveRDS(best_params, file =  paste0("models/hyperparameters/riBART_",i)) # Save hyperparameters
      saveRDS(best_ri_bart, file = model_name) # Save model
  }
  else{
    message(paste0('loading model for month ', i, ' since the model was already trained')) 
    best_ri_bart <- readRDS(file = model_name) 
  }
  
  pred <- pnorm(best_ri_bart$yhat.test.mean)
}

ribart_roc <- function(test, month){
  ribart_fit <- readRDS(file =  paste0("models/riBART/riBART_",month))
 
  bart_predict <- data.frame(truth = as.factor(test$aprobado), 
                             estimated = pnorm(ribart_fit$yhat.test.mean))
  
  roc_curve(bart_predict, truth, estimated, event_level = 'second')
}

plot_random_effects <- function(data, month){
  rbartFit <- readRDS(file =  paste0("models/riBART/riBART_",month))
    
  random_effects_means <- rbartFit$ranef.mean
  
  random_effects_sd <- unlist(lapply(1:n_distinct(cei$id_centro), function(x) sd(unlist(as.list(rbartFit$ranef[,,x])))))
  
  random_effects <- tibble(mean = random_effects_means,  
                           sd = random_effects_sd) |> 
                    mutate(is_significant = as.factor(if_else((mean - sd > 0), 1, if_else((mean + sd < 0),-1,0))))
  
  positive_effect_qty <- nrow(random_effects |> filter(is_significant == 1))
  negative_effect_qty <-  nrow(random_effects |> filter(is_significant == -1))
    
  random_effects$names <- names(random_effects_means)
  
  centros <- cei |> 
    group_by(id_centro) |> 
    summarise(contexto_sociocultural = max(contexto_sociocultural)) |> 
    mutate(id_centro = as.character(id_centro)) |>
    inner_join(random_effects, by = join_by(id_centro == names))
  
  #ggplot(centros, aes(x = reorder(id_centro, mean), y = mean)) +
  ggplot(centros, aes(x = rank(mean), y = mean)) +  
    geom_point(aes(color = is_significant, fill = is_significant), size = 2, shape = 21) +
    # Add error bars for standard errors
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd,  color = as.factor(is_significant)) ) +
    geom_hline(yintercept = 0) +
    labs(y = "Random Effect", x= 'Centers',
         #caption = paste0(positive_effect_qty, " schools with positive effect, ", negative_effect_qty, " schools with negative effect")
         ) +
    scale_color_manual("is_significant", breaks=c(-1,0,1),values=c("#d63031", "grey75", "#00b894")) +
    scale_fill_manual("is_significant", breaks=c(-1,0,1),values=c("#d63031", "grey75", "#00b894"))+
    facet_wrap( ~ contexto_sociocultural, 
                nrow = 2,
                labeller = labeller(contexto_sociocultural = c("Quintil 1" = "Quintile 1",
                                                               "Quintil 2" = "Quintile 2",
                                                               "Quintil 3" = "Quintile 3",
                                                               "Quintil 4" = "Quintile 4",
                                                               "Quintil 5" = "Quintile 5"
                                                              ))) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          #axis.text = element_text(size = 16), 
          #axis.title = element_text(size = 18),
          #plot.caption = element_text(size = 13)
    )
    
}

prob_quantile_activity <- function(train, test, month){
  syntethic_students <- test |>
              select(-id_persona, -starts_with("class_lb_")) |>
              summarise(across(
                where(is.numeric), 
                list(
                  mean = ~ifelse(is.integer(.), round(quantile(., 0.5)), quantile(., 0.5)),
                  d90  = ~ifelse(is.integer(.), round(quantile(., 0.9)), quantile(., 0.9)),
                  d10  = ~ifelse(is.integer(.), round(quantile(., 0.1)), quantile(., 0.1))
                )
              )) |>
              pivot_longer(
                everything(), 
                names_to = c("variable", ".value"), 
                names_pattern = "(.*)_(.*)"
              ) |> 
              column_to_rownames(var = "variable") |> 
              t() |> 
              as_tibble()
  
  syntethic_students$id_persona <- c(1, 2, 3)
  syntethic_students$departamento <- rep("MONTEVIDEO", 3)
  syntethic_students$zona <- rep("Urbana", 3)
  syntethic_students$id_centro <-  c("1261", "1261", "1261")
  syntethic_students$aprobado <- c(0,0,0) #Dummy values
  syntethic_students$nivel <- rep('A2.1', 3) #Dummy values
  syntethic_students[paste0("class_lb_", 3:month)] <- lapply(3:month, function(x) rep(1016, 3))
  
  syntethic_students <- syntethic_students |> 
              mutate_at(vars(id_persona, 
                             aprobado, 
                             contains("acumu"), 
                             starts_with("class_lb_"),
              ), as.integer)
  
  prlist <- readRDS(file = paste0("models/hyperparameters/riBART_",month))
  
  prob_por_contexto <- function(ctxt, best_params){
    syntethic_students$contexto_sociocultural <- rep(ctxt, 3)
    
    syntethic_students <- syntethic_students[, colnames(test)]  #Reorder syntethic_students to match test
    test <- rbind(test, syntethic_students)
    
    cei_train_bart <- train |> select(-c(puntos,nivel, id_persona)) |> mutate(id_centro = as.factor(id_centro))
    cei_test_bart  <- test |> select(-c(puntos,nivel, id_persona)) |> mutate(id_centro = as.factor(id_centro))
    
    rbartFit_2 <-  rbart_vi(aprobado ~ . - id_centro, 
                           data = cei_train_bart,test = cei_test_bart, 
                           group.by = id_centro, group.by.test = id_centro,
                           power = best_params$power, base = best_params$base,
                           k = best_params$k,  n.trees = best_params$n.tree,
                           keepTrees = TRUE, keepCall = TRUE, keepTestFits = TRUE,
                           seed = 1784)
    
    desvios <- c(sd(pnorm(unlist(as.list(rbartFit_2$yhat.test[,,568])))), 
                 sd(pnorm(unlist(as.list(rbartFit_2$yhat.test[,,569])))), 
                 sd(pnorm(unlist(as.list(rbartFit_2$yhat.test[,,570])))))
    
    c <- data.frame(decile = c("d50", "d90", "d10"), contexto = rep(ctxt, 3), 
                    prob = pnorm(rbartFit_2$yhat.test.mean[c(568, 569, 570)]),
                    desv = desvios)
    return(c)
  }
  
  c5 <- prob_por_contexto("Quintil 5", prlist)
  c4 <- prob_por_contexto("Quintil 4", prlist)
  c3 <- prob_por_contexto("Quintil 3", prlist) 
  c2 <- prob_por_contexto("Quintil 2", prlist)
  c1 <- prob_por_contexto("Quintil 1", prlist)
  
  prob_quantiles <- rbind(c5,c4,c3,c2,c1) |> arrange(desc(decile), desc(contexto))
  
  ggplot(prob_quantiles, aes(x = contexto, y = prob, color = decile, group = decile)) + 
    geom_point() + geom_line() +
    geom_errorbar(aes(ymin = prob - desv, ymax = prob + desv),
                  width = 0.2, alpha = 0.7, linewidth = 1) +
    labs(x = "Socioeconomic context",
         y = "Prob. to reach A2.1",
         color = "LB_use") +
    scale_x_discrete(labels = c("Quintil 1" = "Quintile 1",
                                "Quintil 2" = "Quintile 2",
                                "Quintil 3" = "Quintile 3",
                                "Quintil 4" = "Quintile 4",
                                "Quintil 5" = "Quintile 5")) +  
    theme_minimal() +
    theme(aspect.ratio = 1,
          legend.position = "top",
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 17),
          legend.title = element_text(size =18), 
          plot.caption = element_text(size = 13, hjust = 0.5))
}

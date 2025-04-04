library(tidymodels)
library(data.table)
library(dbarts)
library(ParBayesianOptimization)

source('Rcode_paper/random_forest.R')
source('Rcode_paper/ri_bart.R')

set.seed(1780)

cei <- fread('data_output/datos_6to_bruno.csv')
cei_split <- initial_split(cei, prop = 0.80, strata = aprobado)
cei_train <- training(cei_split) 
cei_test <- testing(cei_split) 

dfs <- c()

for(i in 3:11){
  j <- i+1
  regex <- paste0(paste(j:12, collapse = "$|"),"$")
  dfs[[i-2]] <- list(cei_train |> select(!matches(regex)),
                     cei_test |> select(!matches(regex)))
}

col_names <- sapply(3:11, function(x) paste0(x))
met <- metric_set(accuracy, sens, spec)

getResults <- function(model){
  f <- NULL
  if(model == 'rf') 
    f <- function(train, test) random_forest(train, test)
  else if(model == 'rf_default')
    f <- function(train, test) random_forest(train, test, default = T)
  else if(model == 'riBART')
    f <- function(train, test) ri_bart(train, test)
  else if(model == 'riBART_default')
    f <- function(train, test) ri_bart(train, test, default = T)
  else
    warning("No model matches")

  results <- lapply(dfs, function(x){
    f(train = x[[1]], test = x[[2]])
  })
  
  results <- as.data.frame(results)
  colnames(results) <- col_names
  
  results <- cei_test |>
    select(aprobado) |>
    cbind(results) |>
    mutate(across(everything(), function(x) factor(x,  levels=c(1,0))))
  
  result_metrics <- lapply(col_names, function(x){
    results |> 
      met(truth = aprobado, estimate = x) |>
      mutate(month = c(x), model = model) |>
      select(-c(.estimator)) |>
      pivot_wider(names_from = .metric, values_from = .estimate) 
  })
  
  result_metrics <- map_dfr(result_metrics, bind_rows) 
}

rf <- getResults('rf')
rf_default <- getResults('rf_default')
riBART <- getResults('riBART')
riBART_default <- getResults('riBART_default')

######### SHOW RESULTS ###########
combined_metrics <- rbind(rf, rf_default, riBART, riBART_default) |> 
                      pivot_longer(!c(month, model), names_to = 'metric', values_to = 'value') |>
                      mutate(month = as.integer(month))
  

ggplot(combined_metrics, aes(x = month, y = value, color = metric, linetype = metric)) +
  geom_line() +
  geom_point() +
  labs(x = "Month", y = "Value", color = "Model", linetype = "Metric") +
  scale_x_continuous(breaks= seq(3, 11, by=1)) +
  facet_wrap( ~ model, nrow = 2) 

######## Work only with month 7 #############
month <- 3

## ROC Curve
roc_rf <- random_forest_roc(test = dfs[[month - 2]][[2]], month)
roc_bart <- ribart_roc(test = dfs[[month - 2]][[2]], month)

roc_combinadas <- data.frame(
  spec_1 = c(1 - roc_rf$specificity, 1 - roc_bart$specificity),
  sens = c(roc_rf$sensitivity, roc_bart$sensitivity),
  Model = c(rep("Random Forest", length(roc_rf$sensitivity)), rep("riBART", length(roc_bart$sensitivity)))
)

ggplot(roc_combinadas, aes(x = spec_1, y = sens, color = Model)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "1 - Specificity",
       y = "Sensitivity",
       color = "Model") +
  scale_color_manual(values = c("Random Forest" = "blue", "riBART" = "red")) +
  theme(
    axis.text = element_text(size = 16), 
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 13),
    legend.title = element_text(size =14),
    legend.position = "top"
  )

## Random effects
plot_random_effects(cei, month)


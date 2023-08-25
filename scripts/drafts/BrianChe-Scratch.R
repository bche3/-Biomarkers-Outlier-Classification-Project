library(tidyverse)
library(ggplot2)
library(dlookr)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(discrim)

# Data Preprocessing

# get names
var_names <- read_csv('data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# function for trimming outliers (good idea??)
trim <- function(x, .at){
  x[abs(x) > .at] <- sign(x[abs(x) > .at])*.at
  return(x)
}

# read in data
biomarker_clean_trimmed <- read_csv('data/biomarker-raw.csv', 
                                    skip = 2,
                                    col_select = -2L,
                                    col_names = c('group', 
                                                  'empty',
                                                  pull(var_names, abbreviation),
                                                  'ados'),
                                    na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # log transform, center and scale, and trim
  mutate(across(.cols = -c(group, ados), 
                ~ trim(scale(log10(.x))[, 1], .at = 3))) %>%
  # reorder columns
  select(group, ados, everything())
# read in data

biomarker_clean_trimmed


# Question 2: Temporarily remove the outlier trimming from preprocessing and do some exploratory analysis of outlying values. Are there specific subjects (not values) 
# that seem to be outliers? If so, are outliers more frequent in one group or the other? (Hint: consider tabluating the number of outlying values per subject.)

# Remove outlier trimming
biomarker_clean <- read_csv('data/biomarker-raw.csv', 
                            skip = 2,
                            col_select = -2L,
                            col_names = c('group', 
                                          'empty',
                                          pull(var_names, abbreviation),
                                          'ados'),
                            na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # log transform, center and scale, and trim
  mutate(across(.cols = -c(group, ados), 
                ~ scale(log10(.x))[, 1], .at = 3)) %>%
  # reorder columns
  select(group, ados, everything())

biomarker_clean %>%
  plot_normality(CHIP, CEBPB)




biomarker_transformed <- read_csv('data/biomarker-raw.csv', 
                                  skip = 2,
                                  col_select = -2L,
                                  col_names = c('group', 
                                                'empty',
                                                pull(var_names, abbreviation),
                                                'ados'),
                                  na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  mutate(across(.cols = -c(group, ados), log10)) %>%
  mutate(across(.cols = -c(group, ados), scale)) %>%
  
  select(group, ados, everything())

biomarker_transformed
biomarker_clean

# EDA of outliers

dlookr::diagnose_outlier(biomarker_clean)


# Question 4: Use any method to find either:
# - a simpler panel that achieves comparable classification accuracy
# - an alternative panel that achieves improved classification accuracy
# Benchmark your results against the in-class analysis.


#### Using linear discriminant analysis model (LDA) for classification


s_star <- c("DERM", "RELT", "IgD", "PTN", "FSTL1")
biomarker_clean_2 <- biomarker_transformed %>%
  select(group, any_of(s_star)) %>%
  # convert group (chr) to binary (lgl)
  mutate(class = (group == 'ASD')) %>%
  select(-group)


# split the data into training and testing sets

set.seed(101422)
biomarker_split <- initial_split(biomarker_clean_2, prop = 0.8)
biomarker_train <- training(biomarker_split)
biomarker_test <- testing(biomarker_split)

# Convert the outcome variable "class" to factor from logical
biomarker_train$class <- as.factor(biomarker_train$class)
biomarker_test$class <- as.factor(biomarker_test$class)
class(biomarker_train$class)
class(biomarker_test$class)

# creating recipe
biomarker_recipe <- recipe(class ~ ., data = biomarker_train)

lda_mod = discrim_linear() %>%
  set_mode("classification") %>%
  set_engine("MASS")

lda_wkflow = workflow() %>%
  add_model(lda_mod) %>%
  add_recipe(biomarker_recipe)

lda_fit = fit(lda_wkflow, biomarker_train)

predict(lda_fit, new_data = biomarker_train, type = "prob")

augment(lda_fit, new_data = biomarker_train) %>%
  conf_mat(truth = class, estimate = .pred_class) 

lda_acc_train <- augment(lda_fit, new_data = biomarker_train) %>%
  accuracy(truth = class, estimate = .pred_class)

lda_acc_train

#### Using Quadratic Discriminant Analysis (QDA) for classification on training set
qda_mod <- discrim_quad() %>% 
  set_mode("classification") %>% 
  set_engine("MASS")

qda_wkflow <- workflow() %>% 
  add_model(qda_mod) %>% 
  add_recipe(biomarker_recipe)

qda_fit = fit(qda_wkflow, data= biomarker_train)

predict(qda_fit, new_data = biomarker_train, type = "prob")

augment(qda_fit, new_data = biomarker_train) %>%
  conf_mat(truth = class, estimate = .pred_class) 

qda_acc_train <- augment(qda_fit, new_data = biomarker_train) %>%
  accuracy(truth = class, estimate = .pred_class)


#### Using Naive Bayes classification on training set

nb_mod <- naive_Bayes() %>% 
  set_mode("classification") %>% 
  set_engine("klaR")# %>% 
# set_args(usekernel = FALSE) 

nb_wkflow <- workflow() %>% 
  add_model(nb_mod) %>% 
  add_recipe(biomarker_recipe)

nb_fit <- fit(nb_wkflow, biomarker_train)

predict(nb_fit, new_data = biomarker_train, type = "prob")

augment(nb_fit, new_data = biomarker_train) %>%
  conf_mat(truth = class, estimate = .pred_class) 

nb_acc_train <- augment(nb_fit, new_data = biomarker_train) %>%
  accuracy(truth = class, estimate = .pred_class)

nb_acc_train



#### Fitting to test
### Using LDA on testing set
lda_fit_test <- fit(lda_wkflow, biomarker_test)

predict(lda_fit_test, new_data = biomarker_test, type = "prob")

augment(lda_fit_test, new_data = biomarker_test) %>%
  conf_mat(truth = class, estimate = .pred_class) 

multi_metric <- metric_set(accuracy, sensitivity, specificity)

lda_acc_test <- augment(lda_fit_test, new_data = biomarker_test) %>%
  multi_metric(truth = class, estimate = .pred_class)
lda_acc_test

### Using QDA on testing set
qda_fit_test <- fit(qda_wkflow, biomarker_test)

predict(qda_fit_test, new_data = biomarker_test, type = "prob")

augment(qda_fit_test, new_data = biomarker_test) %>%
  conf_mat(truth = class, estimate = .pred_class) 

qda_acc_test <- augment(qda_fit_test, new_data = biomarker_test) %>%
  accuracy(truth = class, estimate = .pred_class)

multi_metric <- metric_set(accuracy, sensitivity, specificity)

qda_acc_test <- augment(qda_fit_test, new_data = biomarker_test) %>%
  multi_metric(truth = class, estimate = .pred_class)
qda_acc_test


### Using Naive Bayes classification on testing set
nb_fit_test <- fit(nb_wkflow, biomarker_test)

predict(nb_fit_test, new_data = biomarker_test, type = "prob")

multi_metric <- metric_set(accuracy, sensitivity, specificity)

nb_acc_test <- augment(nb_fit_test, new_data = biomarker_test) %>%
  multi_metric(truth = class, estimate = .pred_class)
nb_acc_test

lda_acc_test
nb_acc_test

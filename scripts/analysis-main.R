### Question 1 ###
library(tidyverse)


# get names
var_names <- read_csv("/Users/shannon/Documents/PSTAT197/biomarkers-group-1/data/biomarker-raw.csv", 
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
biomarker_clean <- read_csv("/Users/shannon/Documents/PSTAT197/biomarkers-group-1/data/biomarker-raw.csv", 
                            skip = 2,
                            col_select = -2L,
                            col_names = c('group', 
                                          'empty',
                                          pull(var_names, abbreviation),
                                          'ados'),
                            na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # reorder columns
  select(group, ados, everything())

ggplot(biomarker_clean) + geom_histogram(aes(x = CHIP))
ggplot(biomarker_clean) + geom_histogram(aes(x = CEBPB))
ggplot(biomarker_clean) + geom_histogram(aes(x = NSE))
ggplot(biomarker_clean) + geom_histogram(aes(x = PIAS4))


library(tidyverse)


# get names
var_names <- read_csv("/Users/shannon/Documents/PSTAT197/biomarkers-group-1/data/biomarker-raw.csv", 
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
biomarker_clean <- read_csv("/Users/shannon/Documents/PSTAT197/biomarkers-group-1/data/biomarker-raw.csv", 
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

ggplot(biomarker_clean) + geom_histogram(aes(x = CHIP))
ggplot(biomarker_clean) + geom_histogram(aes(x = CEBPB))
ggplot(biomarker_clean) + geom_histogram(aes(x = NSE))
ggplot(biomarker_clean) + geom_histogram(aes(x = PIAS4))

### Question 2 ###
#Reading in data
var_names <- read_csv('C:/Users/longw/OneDrive/Desktop/biomarkers-group-1/data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()


#Pre-processing data but keeping outliers
biomarker_outliers <- read_csv("C:/Users/longw/OneDrive/Desktop/biomarkers-group-1/data/biomarker-raw.csv", 
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


#Adding variable to track number of outlying protein values per subject
biomarker_outliers$outliers <- 0
biomarker_outliers <- biomarker_outliers %>%
  select(group, ados, outliers, everything())

#Nested loop to count outliers
#Execution takes a while, maybe more efficient method?
#X = observation number
#Y is variable index
#Y starts at 4, because first 3 vars = (group, ados, outliers)
#Condition for outlier is the same as specified in the original methodology 
for (x in 1:154) {
  for (y in 4:1320) {
    if (abs(biomarker_outliers[x,y]) > 3) {
      biomarker_outliers[x,3] = biomarker_outliers[x,3] + 1
    }  
  }
}

#First try to identify if any particular subjects are outliers, in terms of having an extreme number of outlier protein values

summary(biomarker_outliers$outliers)
#Mean is 15.45, but median is only 8.5
#Some subjects have such a high amount of extreme outlier variables that it is skewing our distribution.

hist(biomarker_outliers$outliers, breaks = seq(0, 160, 10))
#Vast majority of subjects do not have more than 20 outlier proteins in their blood sample

#Knowing that, we'll set our arbitrary outlier value cutoff at 20. Any subject exceeding 20 protein outliers will be deemed an outlier subject. 

#Adding new variable to classify subject as outlier or not
biomarker_outliers$outlier.subject <- FALSE
biomarker_outliers <- biomarker_outliers %>%
  select(group, ados, outliers, outlier.subject, everything())

biomarker_outliers <- biomarker_outliers %>%
  mutate(outlier.subject = ifelse(outliers > 20, TRUE, FALSE))

summary(biomarker_outliers$outlier.subject)
#21 outlier subjects

tapply(biomarker_outliers$outlier.subject, biomarker_outliers$group, summary)


#A few more outliers in the TD group than the ASD group, but it doesn't seem statistically significant

#Viewing outlier subjects
outlier_subjects <- biomarker_outliers[biomarker_outliers$outlier.subject == TRUE,]

outlier_subjects






### Question 3 ###

library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
load('/Users/shannon/Documents/PSTAT197/biomarkers-group-1/data/biomarker-clean.RData')
## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################
# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)
# This gives us 4 proteins: DERM, RELT, IgD, and FSTL1


biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit <- glm(class ~ ., 
           data = training(biomarker_split), 
           family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

testing(biomarker_split) %>%
  add_predictions(fit, type = 'response') %>%
  class_metrics(estimate = factor(pred > 0.5),
                truth = factor(class), pred,
                event_level = 'second')

load('/Users/shannon/Documents/PSTAT197/biomarkers-group-1/data/biomarker-clean.RData')
## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

## RANDOM FOREST
##################

# store predictors and response separately
predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

## LOGISTIC REGRESSION
#######################

set.seed(27)
proteins <- data.frame(sets = c("proteins_s1", "proteins_s2"),
                       elements = c(proteins_s1, proteins_s2),
                       fuzzy = runif(2))
fuzzy <- tidySet(proteins)
proteins_sstar_f <- intersection(fuzzy, sets = c("proteins_s1", "proteins_s2"))
# This gives us 4 proteins: DERM, RELT, IgD, and FSTL1

proteins_sstar_f <- data.frame(proteins_sstar_f)
proteins_sstar_f$elements

biomarker_sstar_2 <- biomarker_clean %>%
  select(group, any_of(proteins_sstar_f$elements)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split_2 <- biomarker_sstar_2 %>%
  initial_split(prop = 0.8)

# fit logistic regression model to training set
fit_2 <- glm(class ~ ., 
             data = training(biomarker_split_2), 
             family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

testing(biomarker_split_2) %>%
  add_predictions(fit_2, type = 'response') %>%
  class_metrics(estimate = factor(pred > 0.5),
                truth = factor(class), pred,
                event_level = 'second')

### QUESTION 4 ###

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



### Using Naive Bayes classification on testing set
nb_fit_test <- fit(nb_wkflow, biomarker_test)

predict(nb_fit_test, new_data = biomarker_test, type = "prob")

multi_metric <- metric_set(accuracy, sensitivity, specificity)

nb_acc_test <- augment(nb_fit_test, new_data = biomarker_test) %>%
  multi_metric(truth = class, estimate = .pred_class)

lda_acc_test
nb_acc_test

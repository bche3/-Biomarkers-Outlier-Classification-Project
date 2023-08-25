# Measuring Classification Accuracy

# load packages
library(tidyverse)
library(tidymodels)
library(modelr)
library(rsample)
library(yardstick)

# read data
url <- 'https://raw.githubusercontent.com/pstat197/pstat197a/main/materials/labs/lab4-logistic/data/biomarker_clean.csv'

s_star <- c("DERM", "RELT", "IgD", "PTN", "FSTL1")
biomarker <- read_csv(url) %>%
  # subset to proteins of interest and group
  select(group, any_of(s_star)) %>%
  # convert group (chr) to binary (lgl)
  mutate(class = (group == 'ASD')) %>%
  select(-group)


### Data Partitioning ###
# setting aside a random subset of observations 
# that will be used only to assess predictive accuracy and not to fit any models

# for reproducibility
set.seed(102022)
# must run set.seed with whatever other code we want to run

# partition data
# this is completely random
partitions <- biomarker %>%
  initial_split(prop = 0.8) # we want 80% of data in training set

# examine
partitions

# To retrieve the data partitions
# training set
training(partitions) %>% head(4)

# testing set
testing(partitions) %>% head(4)

### Model Fitting ###
# fit glm
fit <- glm(class ~ ., 
           data = biomarker, 
           family = binomial(link = "logit"))

# glm function can fit many kinds of generalized linear models (such as logistic)

#parameter estimates
tidy(fit)

# what is the estimated change in P(ASD) associated with a +1SD change in log protein level?
# If we were to increase IgD by +1SD, we will see a 0.662 increase in P(ASD)


# modelr makes it easy to compute predictions for a wide range of model objects
# compute predictions on the test set
testing(partitions) %>%
  add_predictions(fit)


# manually transform to probabilities
testing(partitions) %>%
  add_predictions(fit) %>%
  mutate(probs = 1/(1 + exp(-pred))) %>%
  # we are transforming the data
  select(class, pred, probs) %>%
  head(5)


# predict on scale of response
testing(partitions) %>%
  add_predictions(fit, type = 'response') %>%
  select(class, pred) %>%
  head(5)

# predict classes
testing(partitions) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(pred.class = (pred > 0.5)) %>%
  select(class, pred, pred.class) %>%
  head(5)

### Accuracy Measures ###
# tabulate
testing(partitions) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(pred.class = (pred > 0.5)) %>%
  select(class, pred.class) %>%
  table()

# store predictions as factors
pred_df <- testing(partitions) %>%
  add_predictions(fit, type = 'response') %>%
  mutate(pred.class = (pred > 0.5),
         group = factor(class, labels = c('TD', 'ASD')),
         pred.group = factor(pred.class, labels = c('TD', 'ASD'))) 

# check order of factor levels
pred_df %>% pull(group) %>% levels()

# compute specificity
pred_df %>%
  specificity(truth = group, 
              estimate = pred.group,
              event_level = 'second')

# sensitivity
pred_df %>%
  sensitivity(truth = group,
              estimate = pred.group,
              event_level = 'second')

## Action ##
# compute accuracy
pred_df %>%
  accuracy(truth = group,
              estimate = pred.group,
              event_level = 'second')

# define panel (arguments must be yardstick metric function names)
panel_fn <- metric_set(sensitivity, specificity)

# compute
pred_df %>%
  panel_fn(truth = group, #observed values
           estimate = pred.group, #predicted values
           event_level = 'second')

## Action ##
metrics <- metric_set(precision, recall, f_meas)
pred_df %>% 
  metrics(truth = group,
          estimate = pred.group,
          event_level = 'second')

pred_df %>% conf_mat(truth = group, estimate = pred.group)

library(tidyverse)
# install.packages('infer') # execute once then comment out
# iteration strategies for the ASD data set
# data location
url <- 'https://raw.githubusercontent.com/pstat197/pstat197a/main/materials/labs/lab3-iteration/data/biomarker-clean.csv'

# function for outlier trimming
trim_fn <- function(x){
  x[x > 3] <- 3
  x[x < -3] <- -3
  
  return(x)
}

# read in and preprocess data
asd <- read_csv(url) %>%
  select(-ados) %>%
  # log transform
  mutate(across(.cols = -group, log10)) %>%
  # center and scale
  mutate(across(.cols = -group, ~ scale(.x)[, 1])) %>%
  # trim outliers
  mutate(across(.cols = -group, trim_fn))


# Check your understanding
rslt <- rep(NA, 4)
input_vals <- rnorm(n = 3)
for(i in 1:4){
  rslt[i] <- 2*input_vals[i]
}

rslt

# The reason why this prints out NA is because there are only 3 input_vals but we need it to be length 4

# Multiple Testing with Loops
x <- asd %>% filter(group == 'ASD') %>% pull(CHIP)
y <- asd %>% filter(group == 'TD') %>% pull(CHIP)
t.test(x, y, var.equal = F)

t.test(x, y) %>% str()


#p-value
t.test(x, y, var.equal = F)$p.value

n_tests <- 100
p_vals <- rep(NA, n_tests)
for(i in 1:n_tests){
  x <- asd %>% filter(group == 'ASD') %>% pull(i + 1)
  y <- asd %>% filter(group == 'TD') %>% pull(i + 1)
  p_vals[i] <- t.test(x, y, var.equal = F)$p.value
}

tibble(protein = colnames(asd)[2:(n_tests + 1)],
       p = p_vals)

n_tests <- 100
rslt <- tibble(protein = colnames(asd)[2:(n_tests + 1)],
               p = NA)
for(i in 1:n_tests){
  x <- asd %>% filter(group == 'ASD') %>% pull(i + 1)
  y <- asd %>% filter(group == 'TD') %>% pull(i + 1)
  rslt$p[i] <- t.test(x, y, var.equal = F)$p.value
}

# write a loop that stores both the -values and the estimated differences for the first 50 proteins
for(i in 1:n_tests){x <- asd %>% filter(group == 'ASD') %>% pull(i + 1)
y <- asd %>% filter(group == 'TD') %>% pull(i + 1)
rslt$p[i] <- t.test(x, y, var.equal = F)$p.value
}

# The Apply Family
vals <- rnorm(n = 4)
simple_fn <- function(x){2*x}
lapply(vals, simple_fn)

sapply(vals, simple_fn)

# apply a function to an index set
simple_fn_ix <- function(i){2*vals[i]}
rslt_apply <- sapply(1:length(vals), simple_fn_ix)

# equivalent for loop
rslt_loop <- rep(NA, length(vals))
for(i in 1:length(vals)){
  rslt_loop[i] <- 2*vals[i]
}

# compare
rbind(rslt_loop, rslt_apply)




# number of tests to perform
n_tests <- 100

# convert to a list
asd_list <- asd %>% 
  select(1:(n_tests + 1)) %>%
  pivot_longer(cols = -group,
               names_to = 'protein',
               values_to = 'level') %>%
  group_by(protein) %>%
  group_split()

# first entry in list
asd_list[[1]]

t.test(level ~ group, data = asd_list[[1]])

# p value for ith protein
tt_fn <- function(i){
  t.test(level ~ group, data = asd_list[[i]])$p.value
}

# check
tt_fn(1)

# to preform  many tests
sapply(1:n_tests, tt_fn)


# to time our code (loop vs apply function)
# apply time
start <- Sys.time()
rslt <- sapply(1:n_tests, tt_fn)
end <- Sys.time()

end - start

# loop time
start <- Sys.time()
n_tests <- 100
rslt <- tibble(protein = colnames(asd)[2:(n_tests + 1)],
               p = NA)
for(i in 1:n_tests){
  x <- asd %>% filter(group == 'ASD') %>% pull(i + 1)
  y <- asd %>% filter(group == 'TD') %>% pull(i + 1)
  rslt$p[i] <- t.test(x, y, var.equal = F)$p.value
}
end <- Sys.time()

end - start
# apply can return a p-value and test statistic
tt_fn <- function(i){
  test_rslt <- t.test(level ~ group, data = asd_list[[i]])
  out <- c(pval = test_rslt$p.value, 
           tstat = test_rslt$statistic)
  out
}

tt_fn(1)

# will return a matrix
sapply(1:5, tt_fn) %>% t() %>% as_tibble()


##Use sapply to obtain the estimated differences and standard errors for the groupwise comparisons for the first 50 proteins.
##Arrange the result in a data frame with a column indicating the protein, a column indicating the estimated group difference, and a column indicating the standard error.


# Nesting
asd_nested <- asd %>%
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  nest(data = c(level, group))

asd_nested %>% head(5)

asd_nested %>%
  slice(1L) %>%
  pull(data)


# The Map function
tt_fn <- function(.df){
  t.test(level ~ group, data = .df)
}

rslt <- asd_nested %>%
  slice(1:10) %>%
  mutate(ttest.out = map(data, tt_fn))

rslt

rslt %>% slice(1L) %>% pull(ttest.out)

asd_nested %>% 
  slice(1L) %>% 
  unnest(cols = data) %>% 
  infer::t_test(formula = level ~ group,
                order = c('ASD', 'TD'),
                alternative = 'two-sided',
                var.equal = F)

# wrapper around infer::t_test
tt_fn <- function(.df){
  infer::t_test(.df, 
                formula = level ~ group,
                order = c('ASD', 'TD'),
                alternative = 'two-sided',
                var.equal = F)
}

# compute test results
tt_out <- asd_nested %>%
  slice(1:n_tests) %>%
  mutate(ttest = map(data, tt_fn))

# preview
tt_out %>% head(4)

tt_out %>% 
  unnest(ttest) %>%
  head(4)

# time it
start <- Sys.time()
tt_out <- asd_nested %>%
  slice(1:n_tests) %>%
  mutate(ttest = map(data, tt_fn))
end <- Sys.time()

end - start

# bonferroni correction
# 
tt_out %>% 
  unnest(ttest) %>%
  mutate(p_adj = p_value*n_tests) %>%
  select(protein, p_value, p_adj) %>%
  arrange(p_adj) %>%
  head(4)





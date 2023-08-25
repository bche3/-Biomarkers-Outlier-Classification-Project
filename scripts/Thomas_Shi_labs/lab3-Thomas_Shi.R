library(tidyverse)

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

for(i in 1:4){
  print(2*i)
}

flag_vals <- c(1, 2, 3, 4)
for(i in flag_vals){
  out <- 2*i
  print(out)
}

rslt <- rep(NA, 4)
for(i in 1:4){
  rslt[i] <- 2*i
}
rslt

slt <- rep(NA, 4)
input_vals <- c(15, 27, 3, 12.6)
for(i in 1:4){
  rslt[i] <- 2*input_vals[i]
}
rslt

rslt <- rep(NA, 4)
input_vals <- rnorm(n = 3)
for(i in 1:4){
  rslt[i] <- 2*input_vals[i]
}

rslt

x <- asd %>% filter(group == 'ASD') %>% pull(CHIP)
y <- asd %>% filter(group == 'TD') %>% pull(CHIP)
t.test(x, y, var.equal = F)

t.test(x, y) %>% str()

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

p_vals <- rep(NA, 50)
diff <- rep(NA, 50)
for(i in 1:50){
  x <- asd %>% filter(group == 'ASD') %>% pull(i + 1)
  y <- asd %>% filter(group == 'TD') %>% pull(i + 1)
  p_vals[i] <- t.test(x, y, var.equal = F)$p.value
  diff[i] <- t.test(x,y,var.equal = F)$estimate[2]-
    t.test(x,y,var.equal = F)$estimate[1]
}
tibble(protein = colnames(asd)[2:(50 + 1)],
       p = p_vals, diff = diff)



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
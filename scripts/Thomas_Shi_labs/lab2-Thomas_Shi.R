library(tidyverse)

# retrieve class survey data
url <- 'https://raw.githubusercontent.com/pstat197/pstat197a/main/materials/labs/lab2-tidyverse/data/'

background <- paste(url, 'background-clean.csv', sep = '') %>%
  read_csv()

interest <- paste(url, 'interest-clean.csv', sep = '') %>%
  read_csv()

metadata <- paste(url, 'survey-metadata.csv', sep = '') %>%
  read_csv()

# filter rows
background %>%
  filter(math.comf > 3)

# select a column
background %>%
  select(math.comf)

# pull a column
background %>%
  pull(rsrch)

# define a new variable
background %>%
  mutate(avg.comf = (math.comf + prog.comf + stat.comf)/3)

# sequence of verbs
background %>%
  filter(stat.prof == 'Adv') %>%
  mutate(avg.comf = (math.comf + prog.comf + stat.comf)/3) %>%
  select(avg.comf, rsrch)

background %>%
  filter(rsrch == TRUE) %>%
  filter(updv.num == '6-8') %>%
  select(contains('prof'))

background %>%
  filter(rsrch == FALSE) %>%
  filter(updv.num == '6-8') %>%
  select(contains('prof'))

# I don't see any obvious differences

# a summary
background %>%
  filter(stat.prof == 'Adv') %>%
  mutate(avg.comf = (math.comf + prog.comf + stat.comf)/3) %>%
  select(avg.comf, rsrch) %>%
  summarize(prop.rsrch = mean(rsrch))

# equivalent
background %>%
  filter(stat.prof == 'Adv') %>%
  mutate(avg.comf = (math.comf + prog.comf + stat.comf)/3) %>%
  select(avg.comf, rsrch) %>%
  pull(rsrch) %>%
  mean()

background %>%
  filter(stat.prof == 'Adv') %>%
  mutate(avg.comf = (math.comf + prog.comf + stat.comf)/3) %>%
  select(avg.comf, rsrch) %>%
  summarize(prop.rsrch = mean(rsrch),
            med.comf = median(avg.comf))

# average comfort levels across all students
background %>%
  select(contains('comf')) %>%
  summarise_all(.funs = mean)

# create a grouping
background %>%
  group_by(stat.prof)

# count observations
background %>%
  group_by(stat.prof) %>%
  count()

# a grouped summary
background %>%
  group_by(stat.prof) %>%
  select(contains('.comf')) %>%
  summarize_all(.funs = mean)



background %>%
  select(contains('.comf')) %>%
  summarise_all(.funs = median)

background %>%
  group_by(updv.num) %>%
  select(contains('.comf')) %>%
  summarise_all(.funs = median)
#Level of comfort of programming and statistic doesn't change with number of upper division class taken
#Level of comfort of Math fluctuate with number of upper courses taken.

# many variables, many summaries
comf_sum <- background %>%
  select(contains('comf')) %>%
  summarise_all(.funs = list(mean = mean, 
                             median = median,
                             min = min, 
                             max = max))

comf_sum

# gather columns into long format
comf_sum %>% gather(stat, val) 

# separate into rows and columns
comf_sum %>%
  gather(stat, val) %>%
  separate(stat, into = c('variable', 'stat'), sep = '_') 

# spread into table
comf_sum %>%
  gather(stat, val) %>%
  separate(stat, into = c('variable', 'stat'), sep = '_') %>%
  spread(stat, val)

# summary of classes taken
classes <- background %>%
  select(11:29) %>%
  mutate_all(~factor(.x, levels = c('no', 'yes'))) %>%
  mutate_all(~as.numeric(.x) - 1) %>%
  summarize_all(mean) %>%
  gather(class, proportion)

classes

classes %>% arrange(desc(proportion))

# plot it
classes %>%
  ggplot(aes(x = proportion, y = class)) +
  geom_point()

fig <- classes %>%
  ggplot(aes(x = proportion, y = reorder(class, proportion))) +
  geom_point()

fig

# adjust labels
fig + labs(x = 'proportion of class', y = '')


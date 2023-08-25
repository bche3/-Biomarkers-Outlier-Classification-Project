url <- 'https://raw.githubusercontent.com/pstat197/pstat197a/main/materials/scripts/package-installs.R'

source(url)

# setup
library(tidyverse)
library(tidytext)
library(tokenizers)
library(textstem)
library(stopwords)
url <- 'https://raw.githubusercontent.com/pstat197/pstat197a/main/materials/labs/lab5-text/data/drseuss.txt'

# read data
seuss_lines <- read_lines(url, skip_empty_rows = T)

seuss_lines %>% head()

# flag lines with a document id
seuss_lines_df <- tibble(line_lag = c(seuss_lines, NA)) %>%
  mutate(flag = str_detect(line_lag, 'Dr. Seuss'),
         line = lag(line_lag, n = 1),
         doc = cumsum(flag)) %>% 
  select(doc, line) %>%
  slice(-1) %>%
  fill(doc)

seuss_lines_df

# grab titles
titles <- seuss_lines_df %>% 
  group_by(doc) %>%
  slice_head() %>%
  pull(line) %>%
  tolower()

# label docs
seuss_lines_df <- seuss_lines_df %>%
  mutate(doc = factor(doc, labels = titles))

seuss_lines_df

# remove header lines (title/author)
seuss_lines_clean <- seuss_lines_df %>%
  group_by(doc) %>%
  mutate(line_num = row_number() - 2) %>%
  filter(line_num > 0)

seuss_lines_df

seuss_lines_clean %>%
  group_by(doc) %>% 
  summarize(num_lines = max(line_num))
            
seuss_lines_clean %>%
  group_by(doc) %>% 
  summarize(num_bump = sum(str_detect(doc, 'bump')))


# collapse lines into one long string
seuss_text <- seuss_lines_clean %>% 
  summarize(text = str_c(line, collapse = ' '))

cat_in_hat <- seuss_text %>% slice(1) %>% pull(text)

cat_in_hat %>%
  str_remove_all('[[:punct:]]') %>%
  tolower()

clean_fn <- function(.text){
  str_remove_all(.text, '[[:punct:]]') %>% tolower()
}

seuss_text_clean <- seuss_text %>%
  mutate(text = clean_fn(text))


stpwrd <- stop_words %>%
  pull(word) %>%
  str_remove_all('[[:punct:]]')

seuss_tokens_long <- seuss_text_clean %>%
  unnest_tokens(output = token, # specifies new column name
                input = text, # specifies column containing text
                token = 'words', # how to tokenize
                stopwords = stpwrd) %>% # optional stopword removal
  mutate(token = lemmatize_words(token)) 


seuss_tokens_long %>%
  group_by(doc) %>%
  count(token) %>%
  summarise(max_token = token[which(n==max(n))])

seuss_tokens_long %>%
  group_by(doc) %>%
  summarise(max_token = token[which(n==max(n))])

seuss_tfidf <- seuss_tokens_long %>%
  count(doc, token) %>%
  bind_tf_idf(term = token,
              document = doc,
              n = n) 

seuss_df <- seuss_tfidf %>%
  pivot_wider(id_cols = doc, 
              names_from = token,
              values_from = tf_idf,
              values_fill = 0)

seuss_df

seuss_tfidf %>%
  group_by(doc) %>%
  slice_max(tf_idf, n = 2)

seuss_tfidf %>%
  group_by(doc) %>%
  slice_max(tf, n = 2)


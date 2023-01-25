# Finding a consensuns for epithelial markers based on my references
library(tidyverse)

markers <- read.csv('../../resources/cell_type_markers_small_intestine_epithelial.csv', na.strings = '')

# Selecting markers that appear in at least 2 datasets
consensus <- markers %>% 
  unite(col = 'all.markers', 2:5, sep = ', ', na.rm = TRUE) %>%
  separate_rows(2, sep = ', ') %>%
  group_by(Cell.Type, all.markers) %>%
  summarise(count = n()) %>%
  filter(count >= 2) %>%
  ungroup() %>%
  select(Cell.Type, all.markers) %>%
  group_by(Cell.Type) %>%
  summarise(consensuns = str_c(unlist(cur_data()), collapse=', '))


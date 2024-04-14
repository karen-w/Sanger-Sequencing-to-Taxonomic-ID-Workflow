#library(tidyverse)
#library(data.table)
#
file_move <- function(table){
  tbl <- read.table(table, sep = "\t", header = TRUE)
  tbl$new_folder <- map(tbl$new_folder, function(x) ifelse(str_sub(x, -1, -1) == "/", x, paste0(x, "/"))) %>% unlist()
  file_names <- tbl$rename %>% str_split("\\/") %>% map(function(x) x[length(x)]) %>% unlist()
  new_folder_files <- map2(tbl$new_folder, file_names, function(x, y) paste0(x, y)) %>% unlist()
  walk(tbl$new_folder %>% unique(), function(x) dir.create(path = x, showWarnings = FALSE))
  walk2(tbl$rename, new_folder_files, file.rename)
}
#file_move("D:/MSL Projects/ECF Seaslug/BGI_Batch1_20240326/Rename2/BGI_file_rename.tsv")

#library(tidyverse)
#library(data.table)
#
rename_files <- function(input_dir, output_dir, incl_primer = TRUE, output_tbl = FALSE){
  #input_dir <- "D:/MSL Projects/ECF Seaslug/BGI_Batch1_20240326/Success";output_dir <- "D:/MSL Projects/ECF Seaslug/BGI_Batch1_20240326/rename"
  if(str_sub(output_dir,-1,-1) != "/"){output_dir <- paste0(output_dir,"/")}
  files <- c(list.files(input_dir, full.names = TRUE) %>% .[which(str_detect(., "ab1"))],
             list.files(input_dir, full.names = TRUE) %>% .[which(str_detect(., "seq"))])
  files0 <- c(list.files(input_dir) %>% .[which(str_detect(., "ab1"))],
             list.files(input_dir) %>% .[which(str_detect(., "seq"))])
  #
  tbl <- data.frame(file = files,
                    file_type = str_sub(files0, -4, -1),
                    smp = str_split_i(files0, "__", 1),
                    primer = str_split_i(files0, "__", 2) %>% str_split_i("\\)", 2)) %>%
    mutate(rename = paste0(output_dir, smp,"_", primer, file_type))
  #
  if(incl_primer == FALSE){
    tbl <- data.frame(file = files,
                      file_type = str_sub(files0, -4, -1),
                      smp = str_split_i(files0, "__", 1)) %>%
      mutate(rename = paste0(output_dir, smp, file_type))
  }
  #
  dir.create(output_dir)
  walk2(tbl$file, tbl$rename, file.copy)
  #
  if (output_tbl == TRUE){
    setwd(output_dir)
    write.table(tbl, "BGI_file_rename.tsv", sep = "\t", row.name = FALSE, quote = FALSE)
  }
}


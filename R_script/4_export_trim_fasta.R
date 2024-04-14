#library(tidyverse)
#library(sangerseqR)
#library(Biostrings)
#
export_fa_fm_tbl <- function(trim_tbl, file_name = NULL){
  fa_vec <- c()
  for (i in 1:nrow(trim_tbl)){
    file_abi <- trim_tbl$file[i]
    basecall_ratio <- trim_tbl$basecall_ratio[i]
    chrom <- sangerseqR::readsangerseq(file_abi)
    chrom_bs <- makeBaseCalls(chrom, ratio = basecall_ratio)
    trim5 <- trim_tbl$trim5[i]
    trim3 <- trim_tbl$trim3[i]
    fa_vec <- c(fa_vec, subseq(chrom_bs@primarySeq, start = trim5+1, end = nchar(chrom_bs@primarySeq)-trim3))
  }
  fa_set <- fa_vec %>% DNAStringSet()
  names(fa_set) <-   map(trim_tbl$file, function(x) str_split(x, "\\/", simplify = TRUE)[length(str_split(x, "\\/", simplify = TRUE))]) %>% 
    unlist() %>%
    str_split_i("\\.", i = 1)
  #
  if(!is.null(file_name)){
    writeXStringSet(fa_set, filepath = file_name)
  }
  #
  return(fa_set)
}
#fa_test <- export_fa_fm_tbl(trim_df, file_name = "NUD_Batch1_16S.fasta") #test
#library(tidyverse)
#library(sangerseqR)
#
gen_trim_tbl <- function(file_abi, basecall_ratio = 0.33,
                         trimmable_5 = 0.1,
                         trimmable_3 = 0.05,
                         low_peak_pct = 0.1, output_chrom_name = NULL){
  chrom <- sangerseqR::readsangerseq(file_abi)
  chrom_bs <- makeBaseCalls(chrom, ratio = basecall_ratio)
  chrom_bs_traceMat <- chrom_bs@traceMatrix
  chrom_bs_posMat <- chrom_bs@peakPosMatrix
  chrom_bs_ampMat <- chrom_bs@peakAmpMatrix
  # 
  seq_len <- chrom@primarySeq %>% length()
  trim_allow_region <- list(trim5 = 1:floor(seq_len*trimmable_5),
                            trim3 = ceiling(seq_len*(1-trimmable_3)):seq_len)
  #
  peak_max <- map(1:nrow(chrom_bs_ampMat), function(x) chrom_bs_ampMat[x, ] %>% max()) %>% unlist()
  low_peak <- which(peak_max < quantile(peak_max, low_peak_pct))  
  #
  pwa <- Biostrings::pairwiseAlignment(chrom@primarySeq, chrom_bs@primarySeq)
  mismatch <- pwa@subject@mismatch[[1]]
  #
  trim5_vec <- map(trim_allow_region[[1]], function(x) sum(x %in% low_peak, x %in% mismatch)) %>% unlist()
  trim3_vec <- map(trim_allow_region[[2]], function(x) sum(x %in% low_peak, x %in% mismatch)) %>% unlist()
  #
  if(sum(trim5_vec == 2) > 0){
    trim5_pos <- max(which(trim5_vec == 2))
    if (trim5_pos < seq_len*0.03) {
      trim5_pos <- max(which(trim5_vec == 1))
    }
  } else if (sum(trim5_vec == 2) == 0 & sum(trim5_vec == 1) > 0){trim5_pos <- max(which(trim5_vec == 1))
  } else if (sum(trim5_vec == 2) == 0 & sum(trim5_vec == 1) == 0){trim5_pos <- 0}
  #
  if(sum(trim3_vec == 2) > 0){trim3_pos <- min(which(trim3_vec == 2))
  } else if (sum(trim3_vec == 2) == 0 & sum(trim3_vec == 1) > 0){trim3_pos <- min(which(trim3_vec == 1))
  } else if (sum(trim3_vec == 2) == 0 & sum(trim3_vec == 1) == 0){trim3_pos <- 0}
  #
  out <- data.frame(
    file = file_abi,
    trim5 = trim5_pos,
    trim3 = trim3_pos,
    basecall_ratio = basecall_ratio,
    low_peak_pct = low_peak_pct,
    trimmable_5 = trimmable_5,
    trimmable_3 = trimmable_3,
    score = pwa@score,
    trimmable_5_pos = floor(seq_len*trimmable_5),
    trimmable_3_pos = ceiling(seq_len*(1-trimmable_3)),
    mismatch_pos = paste(pwa@subject@mismatch[[1]], collapse = ";"),
    low_peak_pos = paste(low_peak, collapse = ";")
  )
  #
  if (!is.null(output_chrom_name)){
    chromatogram(chrom_bs, trim5 = trim5_pos, trim3 = trim3_pos, 
                 width = 100, height = 2, showcalls = "primary", showtrim = TRUE, 
                 filename = output_chrom_name)
  }
  #
  return(out)
}
#file_names <- list.files("D:/MSL Projects/ECF Seaslug/BGI_Batch1_20240326/Rename/16S/", full.names = TRUE) %>% .[str_detect(., ".ab1")] %>% .[1:10]
#trim_df <- map(file_names, gen_trim_tbl) %>% rbindlist()
#save("trim_df", file = "D:/MSL Miscellaneous/20240408 BGI Pipeline Github/example_trim_df.Rdata")
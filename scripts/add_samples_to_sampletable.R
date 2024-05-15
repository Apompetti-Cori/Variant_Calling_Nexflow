library(tidyverse)
init_here <- function() {
  `%noin%` = Negate(`%in%`)
  files <- dir( all.files = T )
  while ( ".here" %noin% files & getwd()!="/" ) {
    setwd("..")
    files <- dir( all.files = T )
  }
  i_am(".here")
}

init_here()

fastq <- list.files(here::here("fastq"), full.names = T)
sample_ids <- stringr::str_match(basename(fastq), pattern = "(.*)_CKDL*")[,2] %>% unique()
sample_table_template <- readr::read_csv(here::here("sample_table_template.csv"))
sample_table_template[1:length(sample_ids),] <- NA
sample_table_template$sample <- sample_ids
sample_table_template <- as.data.frame(sample_table_template)
rownames(sample_table_template) <- sample_table_template$sample

for(sample in sample_table_template$sample){
  hit <- print(grep(pattern = sample, x = fastq, value = TRUE))
  sample_table_template[sample,]$r1_L1 <- hit[1]
  sample_table_template[sample,]$r2_L1 <- hit[2]
}

sample_table_template[is.na(sample_table_template)] <- ""

write.csv(sample_table_template, file = here::here("sample_table.csv"), quote = F, row.names = F)
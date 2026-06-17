#===================================================================
# Script: Metabolite Annotation & KEGG Name Fetching
#===================================================================

#---- 0. Load packages ----
library(tidyverse)

#---- 1. Set paths ----
base_path <- "/Users/wangjie/Desktop/LAB/1_Data/5.1_ (Re-analyze) Metabolomic/0_Raw Data"

input_path <- file.path(base_path, "1_Data(reanalysis)_Celldeath_Sup_Each_Filtered")
output_path <- file.path(base_path, "5_Metabolite Annotated")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

#---- 2. Import Final Filtered Measurement Data ----
NEG_measure_final <- read_csv(
  file.path(input_path, "17_NEG_measurement_with_KEGG_mapped_only_15ppm_top10.csv"),
  show_col_types = FALSE
)

POS_measure_final <- read_csv(
  file.path(input_path, "17_POS_measurement_with_KEGG_mapped_only_15ppm_top10.csv"),
  show_col_types = FALSE
)

#---- 3. Fetch KEGG Compound Names ----
kegg_compound_raw <- readLines("https://rest.kegg.jp/list/compound", warn = FALSE)

kegg_compound_names <- tibble(raw = kegg_compound_raw) %>%
  separate(raw, into = c("KEGG_ID_raw", "Metabolite_Name_Full"), sep = "\t", extra = "merge") %>%
  mutate(
    KEGG_ID = str_replace(KEGG_ID_raw, "^cpd:", ""),
    # 只提取第一格 Primary name
    Primary_Name = str_squish(str_split(Metabolite_Name_Full, ";", simplify = TRUE)[, 1])
  ) %>%
  select(KEGG_ID, Primary_Name)

#---- 4. Expand, Map, and Collapse Function ----
# 邏輯：把 "C05984;C05988" 拆成兩列，分別去查字典，查完再重新合併成 "Name1;Name2"
map_multiple_ids <- function(df, dict) {
  
  # 給每一列一個獨立 Row_ID
  df_temp <- df %>% mutate(Row_ID = row_number())
  
  # 展開多個 KEGG ID 去對照
  mapped_names <- df_temp %>%
    select(Row_ID, KEGG_ID) %>%
    filter(!is.na(KEGG_ID)) %>%
    separate_rows(KEGG_ID, sep = ";") %>%
    mutate(KEGG_ID = str_squish(KEGG_ID)) %>%
    left_join(dict, by = "KEGG_ID") %>%
    group_by(Row_ID) %>%
    summarise(
      # 將查到的結果用分號合併，找不到的標記為 Unknown
      Primary_Name = paste(replace_na(Primary_Name, "Unknown"), collapse = ";"),
      .groups = "drop"
    )
  
  # 將對照結果合併回原表
  df_temp %>%
    left_join(mapped_names, by = "Row_ID") %>%
    relocate(Primary_Name, .after = KEGG_ID) %>%
    select(-Row_ID)
}

#---- 5. Map Names & Export ----
NEG_final_annotated <- map_multiple_ids(NEG_measure_final, kegg_compound_names)
POS_final_annotated <- map_multiple_ids(POS_measure_final, kegg_compound_names)

write_csv(
  NEG_final_annotated, 
  file.path(output_path, "1_NEG_Final_Measurement_with_KEGG_Names.csv")
)

write_csv(
  POS_final_annotated, 
  file.path(output_path, "1_POS_Final_Measurement_with_KEGG_Names.csv")
)

write_csv(
  kegg_compound_names, 
  file.path(output_path, "0_KEGG_Compound_Dictionary.csv")
)
#==============================================================
# KEGG pathway enrichment analysis
# Cell-death metabolites vs DMSO
#
# Analysis version:
#   Single_KEGG + Multiple_KEGG
#   But each Feature_ID keeps only one best KEGG_ID
#
# Purpose:
#   Avoid signal inflation from expanding Multiple_KEGG features.
#
# Input:
#   17_NEG_measurement_with_KEGG_mapped_only_15ppm_top10.csv
#   17_POS_measurement_with_KEGG_mapped_only_15ppm_top10.csv
#
# Output:
#   Pathway enrichment tables with:
#   Condition, Pathway_ID, Pathway, Total, Hit, P_value, FDR, Hit_KEGG_IDs
#==============================================================


#----0. load packages----
# 這一步在做什麼？
# 載入資料整理、統計分析、KEGG REST 讀取需要的 packages。

library(tidyverse)


#----1. set paths----
# 這一步在做什麼？
# 設定 input / output 資料夾。
#
# raw_data_path:
#   你目前 reanalysis 的主資料夾。
#
# input_path:
#   前面整理好的 KEGG-mapped measurement data 所在位置。
#
# output_path:
#   這次 Best_KEGG per feature enrichment 結果會存在這裡。

raw_data_path <- "/Users/wangjie/Desktop/LAB/1_Data/5.1_ (Re-analyze) Metabolomic/0_Raw Data"

input_path <- file.path(
  raw_data_path,
  "1_Data(reanalysis)_Celldeath_Sup_Each_Filtered"
)

output_path <- file.path(
  raw_data_path,
  "4_KEGG_Enrichment_BestKEGG_15ppm_top10"
)

dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

input_path
output_path


#----2. import KEGG-mapped measurement data----
# 這一步在做什麼？
# 讀入已經整理好的主分析資料。
#
# 這一步得到的 Data 代表什麼？
# NEG_measure_KEGG:
#   NEG mode 中，有 intensity 且有 KEGG_ID 的 features。
#
# POS_measure_KEGG:
#   POS mode 中，有 intensity 且有 KEGG_ID 的 features。
#
# 注意：
# 這兩張表裡可能同一個 feature 有 Multiple_KEGG，
# KEGG_ID 可能是 Cxxxxx;Cyyyyy 的形式。

NEG_measure_KEGG <- read_csv(
  file.path(input_path, "17_NEG_measurement_with_KEGG_mapped_only_15ppm_top10.csv"),
  show_col_types = FALSE
)

POS_measure_KEGG <- read_csv(
  file.path(input_path, "17_POS_measurement_with_KEGG_mapped_only_15ppm_top10.csv"),
  show_col_types = FALSE
)


#----3. combine POS and NEG measurement data----
# 這一步在做什麼？
# 把 POS / NEG 合併，並建立 Feature_ID。
#
# 這一步得到的 Data 代表什麼？
# measure_KEGG_all:
#   所有 KEGG-mapped features 的 measurement table。
#
# 注意：
# Feature_ID = Mode + Compound。
# 這樣 POS 和 NEG 不會因為 Compound 名稱相同而混在一起。

measure_KEGG_all <- bind_rows(
  NEG_measure_KEGG %>% mutate(Mode = "NEG"),
  POS_measure_KEGG %>% mutate(Mode = "POS")
) %>%
  mutate(
    Feature_ID = paste(Mode, Compound, sep = "_")
  ) %>%
  distinct(Feature_ID, .keep_all = TRUE)

write_csv(
  measure_KEGG_all,
  file.path(output_path, "1_All_KEGG_mapped_measurement_POS_NEG.csv")
)


#----4. choose one best KEGG ID per feature----
# 這一步在做什麼？
# 1. Single_KEGG features 本來就只有一個 KEGG_ID，直接保留。
# 2. Multiple_KEGG features 可能有多個 KEGG_ID，用 ; 分隔。
# 3. 這一步把 KEGG_ID 拆開後，每個 Feature_ID 只保留一個 best KEGG_ID。
#
# 這一步得到的 Data 代表什麼？
# measure_best_KEGG:
#   每個 feature 只保留一個 KEGG_ID 的主分析資料。
#
# 為什麼要這樣？
# 如果把 Multiple_KEGG 全部展開，會把同一個 intensity 訊號重複計算，
# 造成 pathway enrichment 假性膨脹。
#
# best KEGG 的選擇邏輯：
# 這裡先使用 KEGG_ID 出現順序的第一個。
# 因為在目前 measurement-level table 中，已經沒有 ID-level source ranking。
#
# 風險：
# Multiple_KEGG 的 best ID 不一定是真實身分。
# 所以這版是 expanded coverage / exploratory main，
# 最保守版本仍然是 Single_KEGG only。

measure_best_KEGG <- measure_KEGG_all %>%
  filter(!is.na(KEGG_ID)) %>%
  mutate(
    KEGG_ID_original = as.character(KEGG_ID)
  ) %>%
  separate_rows(KEGG_ID, sep = ";") %>%
  mutate(
    KEGG_ID = str_squish(KEGG_ID),
    KEGG_ID = str_replace(KEGG_ID, "^cpd:", ""),
    KEGG_ID = na_if(KEGG_ID, "")
  ) %>%
  filter(!is.na(KEGG_ID)) %>%
  group_by(Feature_ID) %>%
  arrange(
    KEGG_mapping_status,
    Best_abs_mass_error_ppm,
    desc(Best_fragmentation_score),
    desc(Best_score),
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    KEGG_analysis_strategy = case_when(
      KEGG_mapping_status == "Single_KEGG" ~ "Single_KEGG",
      KEGG_mapping_status == "Multiple_KEGG" ~ "Best_from_Multiple_KEGG",
      TRUE ~ "Other"
    )
  )

write_csv(
  measure_best_KEGG,
  file.path(output_path, "1_Best_KEGG_per_feature_measurement_for_enrichment.csv")
)

measure_best_KEGG_summary <- measure_best_KEGG %>%
  count(Mode, KEGG_analysis_strategy)

write_csv(
  measure_best_KEGG_summary,
  file.path(output_path, "1_Best_KEGG_per_feature_measurement_Summary.csv")
)

measure_best_KEGG_summary


#----5. define sample columns----
# 這一步在做什麼？
# 指定 DMSO 和各死亡條件的 replicate 欄位。
#
# 這一步得到的 Data 代表什麼？
# dmso_cols:
#   control samples。
#
# condition_groups:
#   每一種死亡條件的 sample columns。

dmso_cols <- c("DMSO-1", "DMSO-2", "DMSO-3")

condition_groups <- list(
  STS = c("STS-1", "STS-2", "STS-3"),
  Fas = c("Fas-1", "Fas-2", "Fas-3"),
  RSL = c("RSL-1", "RSL-2", "RSL-3"),
  Nec = c("Nec-1", "Nec-2", "Nec-3")
)

sample_cols <- c(
  dmso_cols,
  unlist(condition_groups)
)

stopifnot(all(sample_cols %in% colnames(measure_best_KEGG)))


#----6. calculate differential metabolites: each condition vs DMSO----
# 這一步在做什麼？
# 對每個 metabolite feature，計算每個死亡條件相對 DMSO 的：
# 1. DMSO mean
# 2. Condition mean
# 3. log2FC
# 4. p value
# 5. FDR
#
# 這一步得到的 Data 代表什麼？
# differential_metabolite_table:
#   每個 feature 在每個死亡條件 vs DMSO 的統計結果。

test_condition_vs_dmso <- function(df, condition_name, condition_cols) {
  
  df %>%
    rowwise() %>%
    mutate(
      DMSO_mean = mean(c_across(all_of(dmso_cols)), na.rm = TRUE),
      Condition_mean = mean(c_across(all_of(condition_cols)), na.rm = TRUE),
      log2FC = log2((Condition_mean + 1) / (DMSO_mean + 1)),
      p_value = tryCatch(
        t.test(
          x = c_across(all_of(condition_cols)),
          y = c_across(all_of(dmso_cols))
        )$p.value,
        error = function(e) NA_real_
      )
    ) %>%
    ungroup() %>%
    mutate(
      Condition = condition_name
    )
}

differential_metabolite_table <- bind_rows(
  lapply(
    names(condition_groups),
    function(cond) {
      test_condition_vs_dmso(
        df = measure_best_KEGG,
        condition_name = cond,
        condition_cols = condition_groups[[cond]]
      )
    }
  )
) %>%
  group_by(Condition) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH")
  ) %>%
  ungroup() %>%
  select(
    Condition,
    Mode,
    Feature_ID,
    Compound,
    KEGG_ID,
    KEGG_ID_original,
    KEGG_mapping_status,
    KEGG_analysis_strategy,
    DMSO_mean,
    Condition_mean,
    log2FC,
    p_value,
    FDR,
    Best_abs_mass_error_ppm,
    Best_score,
    Best_fragmentation_score,
    everything()
  )

write_csv(
  differential_metabolite_table,
  file.path(output_path, "2_Differential_metabolites_Best_KEGG_vs_DMSO.csv")
)


#----7. define high-up metabolite thresholds----
# 這一步在做什麼？
# 設定 high-up metabolite 的篩選門檻。
#
# p005 version:
#   log2FC >= 1 且 p_value < 0.05
#   這是探索性版本，比較容易有 pathway 結果。
#
# FDR01 version:
#   log2FC >= 1 且 FDR < 0.1
#   這是比較嚴格版本，可能結果很少。

log2FC_cutoff <- 1
pvalue_cutoff <- 0.05
FDR_cutoff <- 0.1


#----8. high-up metabolites: p-value version----
# 這一步在做什麼？
# 篩選每個死亡條件中，相對 DMSO 上升的 metabolites。
#
# 這一步得到的 Data 代表什麼？
# high_up_metabolites_p005:
#   log2FC >= 1 且 p_value < 0.05 的 features。

high_up_metabolites_p005 <- differential_metabolite_table %>%
  filter(
    log2FC >= log2FC_cutoff,
    p_value < pvalue_cutoff
  ) %>%
  arrange(Condition, desc(log2FC))

write_csv(
  high_up_metabolites_p005,
  file.path(output_path, "3_High_up_metabolites_log2FC1_p005_Best_KEGG.csv")
)

high_up_metabolites_p005_summary <- high_up_metabolites_p005 %>%
  count(Condition, Mode, KEGG_analysis_strategy)

write_csv(
  high_up_metabolites_p005_summary,
  file.path(output_path, "3_High_up_metabolites_log2FC1_p005_Summary.csv")
)

high_up_metabolites_p005_summary


#----9. high-up metabolites: FDR version----
# 這一步在做什麼？
# 篩選更嚴格的 high-up metabolites。
#
# 這一步得到的 Data 代表什麼？
# high_up_metabolites_FDR01:
#   log2FC >= 1 且 FDR < 0.1 的 features。

high_up_metabolites_FDR01 <- differential_metabolite_table %>%
  filter(
    log2FC >= log2FC_cutoff,
    FDR < FDR_cutoff
  ) %>%
  arrange(Condition, desc(log2FC))

write_csv(
  high_up_metabolites_FDR01,
  file.path(output_path, "3_High_up_metabolites_log2FC1_FDR01_Best_KEGG.csv")
)

high_up_metabolites_FDR01_summary <- high_up_metabolites_FDR01 %>%
  count(Condition, Mode, KEGG_analysis_strategy)

write_csv(
  high_up_metabolites_FDR01_summary,
  file.path(output_path, "3_High_up_metabolites_log2FC1_FDR01_Summary.csv")
)

high_up_metabolites_FDR01_summary


#----10. prepare KEGG background and input lists----
# 這一步在做什麼？
# 建立 enrichment analysis 需要的 background 和 input。
#
# background:
#   所有 Best_KEGG-per-feature 的 KEGG IDs。
#
# input:
#   每個死亡條件 high-up metabolites 的 KEGG IDs。
#
# 為什麼 background 不能用整個 KEGG？
# 因為你的實驗只偵測到其中一部分 metabolites。
# enrichment 的背景應該是「本實驗中有機會被選到的 KEGG metabolites」。

kegg_background <- measure_best_KEGG %>%
  select(Feature_ID, KEGG_ID) %>%
  distinct() %>%
  filter(!is.na(KEGG_ID))

background_KEGG_IDs <- unique(kegg_background$KEGG_ID)

kegg_input_by_condition_p005 <- high_up_metabolites_p005 %>%
  select(Condition, Feature_ID, KEGG_ID) %>%
  distinct() %>%
  filter(!is.na(KEGG_ID))

kegg_input_by_condition_FDR01 <- high_up_metabolites_FDR01 %>%
  select(Condition, Feature_ID, KEGG_ID) %>%
  distinct() %>%
  filter(!is.na(KEGG_ID))

write_csv(
  kegg_background,
  file.path(output_path, "4_KEGG_background_Best_KEGG.csv")
)

write_csv(
  kegg_input_by_condition_p005,
  file.path(output_path, "4_KEGG_input_by_condition_log2FC1_p005_Best_KEGG.csv")
)

write_csv(
  kegg_input_by_condition_FDR01,
  file.path(output_path, "4_KEGG_input_by_condition_log2FC1_FDR01_Best_KEGG.csv")
)


#----11. download KEGG compound-pathway mapping----
# 這一步在做什麼？
# 從 KEGG REST 下載 compound 和 pathway 的對應表。
#
# 重要：
# KEGG link/pathway/compound 回傳的是 reference pathway：
#   path:map00230
# 不是：
#   path:hsa00230
#
# 所以內部 enrichment 用 path:mapXXXXX。
# 最後輸出時才轉成 path:hsaXXXXX 顯示。

kegg_link_url <- "https://rest.kegg.jp/link/pathway/compound"

kegg_compound_pathway_raw <- readLines(kegg_link_url, warn = FALSE)

kegg_compound_pathway <- tibble(raw = kegg_compound_pathway_raw) %>%
  separate(raw, into = c("KEGG_raw", "Pathway_raw"), sep = "\t") %>%
  mutate(
    KEGG_ID = str_replace(KEGG_raw, "^cpd:", ""),
    Pathway_number = str_extract(Pathway_raw, "\\d{5}"),
    Pathway_ID_map = paste0("path:map", Pathway_number),
    Pathway_ID_display = paste0("path:hsa", Pathway_number)
  ) %>%
  select(
    KEGG_ID,
    Pathway_ID_map,
    Pathway_ID_display
  ) %>%
  filter(
    !is.na(KEGG_ID),
    !is.na(Pathway_ID_map)
  ) %>%
  distinct()

write_csv(
  kegg_compound_pathway,
  file.path(output_path, "5_KEGG_compound_to_map_pathway_mapping.csv")
)

kegg_compound_pathway_summary <- kegg_compound_pathway %>%
  summarise(
    Mapping_rows = n(),
    Unique_KEGG_IDs = n_distinct(KEGG_ID),
    Unique_pathways = n_distinct(Pathway_ID_map)
  )

write_csv(
  kegg_compound_pathway_summary,
  file.path(output_path, "5_KEGG_compound_to_map_pathway_mapping_Summary.csv")
)

kegg_compound_pathway_summary


#----12. download KEGG pathway names----
# 這一步在做什麼？
# 下載 KEGG reference pathway 的名稱。
#
# 重要：
# 這裡使用 list/pathway/map。
# 內部 join 使用 path:mapXXXXX。
# 輸出顯示使用 path:hsaXXXXX。
#
# 這一步得到的 Data 代表什麼？
# kegg_pathway_names:
#   Pathway_ID_map = path:mapXXXXX
#   Pathway_ID_display = path:hsaXXXXX
#   Pathway = pathway name

kegg_pathway_name_raw <- readLines(
  "https://rest.kegg.jp/list/pathway/map",
  warn = FALSE
)

kegg_pathway_names <- tibble(raw = kegg_pathway_name_raw) %>%
  separate(raw, into = c("Pathway_raw", "Pathway"), sep = "\t") %>%
  mutate(
    Pathway_number = str_extract(Pathway_raw, "\\d{5}"),
    Pathway_ID_map = paste0("path:map", Pathway_number),
    Pathway_ID_display = paste0("path:hsa", Pathway_number),
    Pathway = str_replace(Pathway, " - Reference pathway$", "")
  ) %>%
  select(
    Pathway_ID_map,
    Pathway_ID_display,
    Pathway
  ) %>%
  filter(
    !is.na(Pathway_ID_map),
    !is.na(Pathway)
  ) %>%
  distinct()

write_csv(
  kegg_pathway_names,
  file.path(output_path, "5_KEGG_map_pathway_names.csv")
)

kegg_pathway_names_summary <- kegg_pathway_names %>%
  summarise(
    Pathway_name_rows = n(),
    Unique_pathways = n_distinct(Pathway_ID_map)
  )

write_csv(
  kegg_pathway_names_summary,
  file.path(output_path, "5_KEGG_map_pathway_names_Summary.csv")
)

kegg_pathway_names_summary


#----13. check KEGG mapping and input before enrichment----
# 這一步在做什麼？
# 在正式 enrichment 前先確認：
# 1. KEGG compound-pathway mapping 不是空的
# 2. background KEGG IDs 有 pathway annotation
# 3. 各 condition input KEGG IDs 有 pathway annotation

kegg_precheck_background <- tibble(
  Background_KEGG_IDs = length(background_KEGG_IDs),
  Background_KEGG_IDs_with_pathway = length(intersect(
    background_KEGG_IDs,
    unique(kegg_compound_pathway$KEGG_ID)
  ))
)

kegg_precheck_input_p005 <- kegg_input_by_condition_p005 %>%
  group_by(Condition) %>%
  summarise(
    Input_KEGG_IDs = n_distinct(KEGG_ID),
    Input_KEGG_IDs_with_pathway = length(intersect(
      unique(KEGG_ID),
      unique(kegg_compound_pathway$KEGG_ID)
    )),
    .groups = "drop"
  )

kegg_precheck_input_FDR01 <- kegg_input_by_condition_FDR01 %>%
  group_by(Condition) %>%
  summarise(
    Input_KEGG_IDs = n_distinct(KEGG_ID),
    Input_KEGG_IDs_with_pathway = length(intersect(
      unique(KEGG_ID),
      unique(kegg_compound_pathway$KEGG_ID)
    )),
    .groups = "drop"
  )

write_csv(
  kegg_precheck_background,
  file.path(output_path, "5_Precheck_background_KEGG_pathway_mapping.csv")
)

write_csv(
  kegg_precheck_input_p005,
  file.path(output_path, "5_Precheck_input_KEGG_pathway_mapping_p005.csv")
)

write_csv(
  kegg_precheck_input_FDR01,
  file.path(output_path, "5_Precheck_input_KEGG_pathway_mapping_FDR01.csv")
)

kegg_precheck_background
kegg_precheck_input_p005
kegg_precheck_input_FDR01


#----14. define pathway enrichment function----
# 這一步在做什麼？
# 建立一個 function，用 hypergeometric test 做 pathway over-representation analysis。
#
# 每個 pathway 的欄位意義：
# Total:
#   background 中屬於這個 pathway 的 KEGG IDs 數量。
#
# Hit:
#   high-up input 中屬於這個 pathway 的 KEGG IDs 數量。
#
# P_value:
#   hypergeometric enrichment p value。
#
# FDR:
#   BH-adjusted p value。
#
# Hit_KEGG_IDs:
#   hit 到這個 pathway 的 KEGG compound IDs。
#
# 注意：
# 內部使用 Pathway_ID_map = path:mapXXXXX。
# 輸出使用 Pathway_ID = path:hsaXXXXX。

run_kegg_enrichment <- function(input_ids, background_ids, compound_pathway_df, pathway_name_df) {
  
  input_ids <- unique(na.omit(input_ids))
  background_ids <- unique(na.omit(background_ids))
  
  input_ids <- intersect(input_ids, background_ids)
  
  pathway_background <- compound_pathway_df %>%
    filter(KEGG_ID %in% background_ids)
  
  pathway_input <- compound_pathway_df %>%
    filter(KEGG_ID %in% input_ids)
  
  N <- length(background_ids)
  n <- length(input_ids)
  
  if (N == 0 | n == 0 | nrow(pathway_background) == 0) {
    return(tibble(
      Pathway_ID = character(),
      Pathway = character(),
      Total = integer(),
      Hit = integer(),
      P_value = numeric(),
      FDR = numeric(),
      Hit_KEGG_IDs = character()
    ))
  }
  
  enrichment_result <- pathway_background %>%
    group_by(Pathway_ID_map, Pathway_ID_display) %>%
    summarise(
      Total = n_distinct(KEGG_ID),
      .groups = "drop"
    ) %>%
    left_join(
      pathway_input %>%
        group_by(Pathway_ID_map, Pathway_ID_display) %>%
        summarise(
          Hit = n_distinct(KEGG_ID),
          Hit_KEGG_IDs = paste(sort(unique(KEGG_ID)), collapse = ";"),
          .groups = "drop"
        ),
      by = c("Pathway_ID_map", "Pathway_ID_display")
    ) %>%
    mutate(
      Hit = replace_na(Hit, 0),
      Hit_KEGG_IDs = replace_na(Hit_KEGG_IDs, ""),
      P_value = phyper(
        q = Hit - 1,
        m = Total,
        n = N - Total,
        k = n,
        lower.tail = FALSE
      )
    ) %>%
    filter(Hit > 0) %>%
    mutate(
      FDR = p.adjust(P_value, method = "BH")
    ) %>%
    left_join(
      pathway_name_df,
      by = c("Pathway_ID_map", "Pathway_ID_display")
    ) %>%
    mutate(
      Pathway = if_else(
        is.na(Pathway),
        "Pathway name not found",
        Pathway
      )
    ) %>%
    select(
      Pathway_ID = Pathway_ID_display,
      Pathway,
      Total,
      Hit,
      P_value,
      FDR,
      Hit_KEGG_IDs
    ) %>%
    arrange(P_value)
  
  return(enrichment_result)
}


#----15. run KEGG enrichment function by condition----
# 這一步在做什麼？
# 建立一個 wrapper function，讓 p005 和 FDR01 兩種 input 都可以重複使用。

run_enrichment_by_condition <- function(kegg_input_by_condition, output_prefix) {
  
  if (nrow(kegg_input_by_condition) == 0) {
    
    empty_result <- tibble(
      Condition = character(),
      Pathway_ID = character(),
      Pathway = character(),
      Total = integer(),
      Hit = integer(),
      P_value = numeric(),
      FDR = numeric(),
      Input_KEGG_ID_count = integer(),
      Hit_KEGG_IDs = character()
    )
    
    write_csv(
      empty_result,
      file.path(output_path, paste0("6_KEGG_enrichment_", output_prefix, "_Best_KEGG.csv"))
    )
    
    return(empty_result)
  }
  
  kegg_enrichment_result <- bind_rows(
    lapply(
      unique(kegg_input_by_condition$Condition),
      function(cond) {
        
        input_ids <- kegg_input_by_condition %>%
          filter(Condition == cond) %>%
          pull(KEGG_ID)
        
        run_kegg_enrichment(
          input_ids = input_ids,
          background_ids = background_KEGG_IDs,
          compound_pathway_df = kegg_compound_pathway,
          pathway_name_df = kegg_pathway_names
        ) %>%
          mutate(
            Condition = cond,
            Input_KEGG_ID_count = length(unique(input_ids))
          ) %>%
          select(
            Condition,
            Pathway_ID,
            Pathway,
            Total,
            Hit,
            P_value,
            FDR,
            Input_KEGG_ID_count,
            Hit_KEGG_IDs
          )
      }
    )
  )
  
  write_csv(
    kegg_enrichment_result,
    file.path(output_path, paste0("6_KEGG_enrichment_", output_prefix, "_Best_KEGG.csv"))
  )
  
  top20_kegg_pathways <- kegg_enrichment_result %>%
    group_by(Condition) %>%
    arrange(P_value, .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    ungroup()
  
  write_csv(
    top20_kegg_pathways,
    file.path(output_path, paste0("7_Top20_KEGG_pathways_by_condition_", output_prefix, "_Best_KEGG.csv"))
  )
  
  concise_kegg_pathway_table <- kegg_enrichment_result %>%
    select(
      Condition,
      Pathway_ID,
      Pathway,
      Total,
      Hit,
      P_value,
      FDR
    ) %>%
    arrange(Condition, P_value)
  
  write_csv(
    concise_kegg_pathway_table,
    file.path(output_path, paste0("8_Concise_KEGG_pathway_table_", output_prefix, "_Best_KEGG.csv"))
  )
  
  return(kegg_enrichment_result)
}


#----16. run KEGG enrichment: log2FC >= 1 and p < 0.05----
# 這一步在做什麼？
# 使用探索性 high-up metabolite list 進行 KEGG pathway enrichment。

kegg_enrichment_p005 <- run_enrichment_by_condition(
  kegg_input_by_condition = kegg_input_by_condition_p005,
  output_prefix = "log2FC1_p005"
)

kegg_enrichment_p005 %>%
  group_by(Condition) %>%
  slice_min(P_value, n = 10, with_ties = FALSE) %>%
  ungroup()


#----17. run KEGG enrichment: log2FC >= 1 and FDR < 0.1----
# 這一步在做什麼？
# 使用更嚴格的 high-up metabolite list 進行 KEGG pathway enrichment。

kegg_enrichment_FDR01 <- run_enrichment_by_condition(
  kegg_input_by_condition = kegg_input_by_condition_FDR01,
  output_prefix = "log2FC1_FDR01"
)

kegg_enrichment_FDR01 %>%
  group_by(Condition) %>%
  slice_min(P_value, n = 10, with_ties = FALSE) %>%
  ungroup()


#----18. final check----
# 這一步在做什麼？
# 檢查 pathway name 是否成功接上。
#
# Missing_pathway_name 應該要是 0。
# 如果不是 0，代表 Pathway_ID 和 pathway name 對照表格式還沒接上。

final_check_p005 <- kegg_enrichment_p005 %>%
  summarise(
    Rows = n(),
    Missing_pathway_name = sum(is.na(Pathway) | Pathway == "Pathway name not found")
  )

final_check_FDR01 <- kegg_enrichment_FDR01 %>%
  summarise(
    Rows = n(),
    Missing_pathway_name = sum(is.na(Pathway) | Pathway == "Pathway name not found")
  )

write_csv(
  final_check_p005,
  file.path(output_path, "9_Final_check_p005.csv")
)

write_csv(
  final_check_FDR01,
  file.path(output_path, "9_Final_check_FDR01.csv")
)

final_check_p005
final_check_FDR01
#===================================================================
# Script: Separate Heatmaps, 2 Names Kept, Pro-Simplification, No Brackets
#===================================================================

#---- 0. Load packages ----
library(tidyverse)
library(pheatmap)

#---- 1. Set paths ----
base_path <- "/Users/wangjie/Desktop/LAB/1_Data/5.1_ (Re-analyze) Metabolomic/0_Raw Data"

data_path <- file.path(base_path, "5_Metabolite Annotated")
output_path <- file.path(base_path, "6_HeatMap")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

#---- 2. Custom Function: 專業代謝物縮寫大師 (Pro Version + 2 Names + 無括號) ----
simplify_metabolite_names_pro <- function(names_vec) {
  
  # 1. 提取前 2 個同義詞，並用 " / " 合併
  primary_names <- sapply(names_vec, function(x) {
    parts <- str_squish(str_split(x, ";")[[1]])
    # 取前兩個，如果只有一個就不會有斜線
    paste(head(parts, 2), collapse = " / ")
  })
  
  # 2. 建立生化標準縮寫字典進行全取代
  clean_names <- primary_names %>%
    # 【暴力清除括號】：把所有 ( ) 以及裡面的東西全部刪除
    str_replace_all("\\(.*?\\)", "") %>%
    
    # 核苷酸與能量分子
    str_replace_all("(?i)\\bAdenosine 5'-triphosphate\\b", "ATP") %>%
    str_replace_all("(?i)\\bAdenosine 5'-diphosphate\\b", "ADP") %>%
    str_replace_all("(?i)\\bAdenosine 5'-monophosphate\\b", "AMP") %>%
    str_replace_all("(?i)\\bUridine 5'-triphosphate\\b", "UTP") %>%
    str_replace_all("(?i)\\bUridine 5'-diphosphate\\b", "UDP") %>%
    str_replace_all("(?i)\\bUridine 5'-monophosphate\\b", "UMP") %>%
    str_replace_all("(?i)\\bGuanosine 5'-triphosphate\\b", "GTP") %>%
    str_replace_all("(?i)\\bGuanosine 5'-diphosphate\\b", "GDP") %>%
    str_replace_all("(?i)\\bGuanosine 5'-monophosphate\\b", "GMP") %>%
    str_replace_all("(?i)\\bCytidine 5'-triphosphate\\b", "CTP") %>%
    str_replace_all("(?i)\\bCytidine 5'-diphosphate\\b", "CDP") %>%
    str_replace_all("(?i)\\bCytidine 5'-monophosphate\\b", "CMP") %>%
    
    # 輔酶群
    str_replace_all("(?i)\\bNicotinamide adenine dinucleotide phosphate\\b", "NADP+") %>%
    str_replace_all("(?i)\\bNicotinamide adenine dinucleotide\\b", "NAD+") %>%
    str_replace_all("(?i)\\bCoenzyme A\\b", "CoA") %>%
    str_replace_all("(?i)Acetyl-CoA", "Ac-CoA") %>%
    
    # 醣類與衍生物 (新增 UDP-Glc)
    str_replace_all("(?i)\\bUDP-glucose\\b", "UDP-Glc") %>%
    str_replace_all("(?i)UDP-N-acetyl-glucosamine", "UDP-GlcNAc") %>%
    str_replace_all("(?i)UDP-N-acetyl-galactosamine", "UDP-GalNAc") %>%
    str_replace_all("(?i)N-Acetyl-glucosamine", "GlcNAc") %>%
    str_replace_all("(?i)N-Acetyl-galactosamine", "GalNAc") %>%
    str_replace_all("(?i)N-Acetyl-neuraminate", "Neu5Ac") %>%
    
    # 脂質常見縮寫
    str_replace_all("(?i)\\bPhosphatidylcholine\\b", "PC") %>%
    str_replace_all("(?i)\\bPhosphatidylethanolamine\\b", "PE") %>%
    str_replace_all("(?i)\\bPhosphatidylserine\\b", "PS") %>%
    str_replace_all("(?i)\\bPhosphatidylinositol\\b", "PI") %>%
    str_replace_all("(?i)\\bSphingomyelin\\b", "SM") %>%
    str_replace_all("(?i)\\bCeramide\\b", "Cer") %>%
    
    # 常見官能基與後綴簡化
    str_replace_all("(?i)\\btriphosphate\\b", "TP") %>%
    str_replace_all("(?i)\\bdiphosphate\\b", "DP") %>%
    str_replace_all("(?i)\\bphosphate\\b", "P") %>%
    str_replace_all("(?i)N-Acetyl-", "NAc-") %>%
    str_replace_all("(?i)\\bacid\\b", "") %>%
    
    # 刪除無意義立體結構與前綴 (alpha, beta, D, L, sn-, etc.)
    str_replace_all("(?i)alpha-D-|beta-D-|alpha-|beta-|D-|L-|sn-", "") %>%
    
    # 最後收尾：去除多餘空白與符號
    str_replace_all("\\s+", " ") %>%
    str_replace_all("-$", "") %>%
    str_replace_all(" / $", "") %>% # 如果斜線在最後面，把它修掉
    str_squish()
  
  return(unname(clean_names))
}

#---- 3. 設定要處理的文件列表 ----
file_list <- list(
  POS = "1_POS_Final_Measurement_with_KEGG_Names.csv",
  NEG = "1_NEG_Final_Measurement_with_KEGG_Names.csv"
)

#---- 4. 主循環：逐一處理 POS 和 NEG 資料 ----
for (mode_name in names(file_list)) {
  
  file_name <- file_list[[mode_name]]
  full_data <- read_csv(file.path(data_path, file_name), show_col_types = FALSE)
  
  # 步驟 A: 專業命名、改名、全 0 過濾
  processed_df <- full_data %>%
    filter(!is.na(Primary_Name) & Primary_Name != "Unknown") %>%
    mutate(Short_Name = simplify_metabolite_names_pro(Primary_Name)) %>%
    rename(
      `Ferroptosis-1` = `RSL-1`,
      `Ferroptosis-2` = `RSL-2`,
      `Ferroptosis-3` = `RSL-3`
    ) %>%
    # 雙 0 過濾：兩組(共6樣本)加總為0則刪除
    filter((`DMSO-1` + `DMSO-2` + `DMSO-3` + `Ferroptosis-1` + `Ferroptosis-2` + `Ferroptosis-3`) > 0)
  
  # 步驟 B: 計算 Fold Change 和 T-test P-value
  stats_df <- processed_df %>%
    mutate(across(c(`DMSO-1`, `DMSO-2`, `DMSO-3`, `Ferroptosis-1`, `Ferroptosis-2`, `Ferroptosis-3`), ~ .x + 1)) %>%
    rowwise() %>%
    mutate(
      Mean_DMSO = mean(c(`DMSO-1`, `DMSO-2`, `DMSO-3`)),
      Mean_Ferro = mean(c(`Ferroptosis-1`, `Ferroptosis-2`, `Ferroptosis-3`)),
      Log2FC = log2(Mean_Ferro / Mean_DMSO),
      P_value = tryCatch(
        t.test(c(`Ferroptosis-1`, `Ferroptosis-2`, `Ferroptosis-3`), c(`DMSO-1`, `DMSO-2`, `DMSO-3`))$p.value, 
        error = function(e) NA
      )
    ) %>%
    ungroup() %>%
    # 💡 改回 2 倍差異 (abs(Log2FC) > 1)
    filter(P_value < 0.05 & abs(Log2FC) > 1)
  
  if(nrow(stats_df) == 0) {
    next
  }
  
  # 處理重複名稱防呆
  final_df <- stats_df %>%
    group_by(Short_Name) %>%
    mutate(Unique_Name = if(n() > 1) paste0(Short_Name, "_", row_number()) else Short_Name) %>%
    ungroup()
  
  # 步驟 C: 繪製 Heatmap
  mat <- final_df %>%
    select(Unique_Name, `DMSO-1`, `DMSO-2`, `DMSO-3`, `Ferroptosis-1`, `Ferroptosis-2`, `Ferroptosis-3`) %>%
    column_to_rownames("Unique_Name") %>%
    as.matrix() %>%
    log2()
  
  # 輸出 PNG (套用你抓好的完美參數：字體 6, 高度 15)
  pheatmap(
    mat,
    scale = "row",
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    show_colnames = FALSE,       
    border_color = NA,           
    fontsize_row = 6,            # 完美參數：字體 6
    color = colorRampPalette(c("blue", "white", "red"))(100), 
    main = paste0("Significant Metabolites (", mode_name, "): Ferroptosis vs DMSO (>2 Fold)"),
    filename = file.path(output_path, paste0("1_Heatmap_", mode_name, "_Ferroptosis_vs_DMSO_2Fold.png")),
    width = 8,
    height = 15                  # 完美參數：高度 15
  )
  
  # 步驟 D: 匯出資料表
  write_csv(
    final_df, 
    file.path(output_path, paste0("1_Significant_Metabolites_", mode_name, "_Simplified_2Fold.csv"))
  )
}